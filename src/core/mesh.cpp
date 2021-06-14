//
// Created by Mike Smith on 2021/6/12.
//

#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/scene.h>
#include <embree3/rtcore.h>

#include <core/mesh.h>

class Accel {

private:
    RTCDevice _device{nullptr};
    RTCScene _scene{nullptr};
    std::span<const glm::vec3> _vertices;
    std::span<const glm::uvec3> _indices;

private:
    Accel(std::span<const glm::vec3> vertices, std::span<const glm::uvec3> indices) noexcept {

        _device = rtcNewDevice(nullptr);
        _scene = rtcNewScene(_device);
        rtcSetSceneBuildQuality(_scene, RTC_BUILD_QUALITY_HIGH);

        auto geometry = rtcNewGeometry(_device, RTC_GEOMETRY_TYPE_TRIANGLE);
        auto vb = static_cast<glm::vec3 *>(rtcSetNewGeometryBuffer(
            geometry, RTC_BUFFER_TYPE_VERTEX, 0,
            RTC_FORMAT_FLOAT3, sizeof(glm::vec3), vertices.size()));
        std::copy_n(vertices.data(), vertices.size(), vb);
        _vertices = {vb, vertices.size()};

        auto ib = static_cast<glm::uvec3 *>(rtcSetNewGeometryBuffer(
            geometry, RTC_BUFFER_TYPE_INDEX, 0,
            RTC_FORMAT_UINT3, sizeof(glm::uvec3), indices.size()));
        std::copy_n(indices.data(), indices.size(), ib);
        _indices = {ib, indices.size()};

        auto t0 = std::chrono::high_resolution_clock::now();
        rtcCommitGeometry(geometry);
        rtcAttachGeometry(_scene, geometry);
        rtcReleaseGeometry(geometry);
        rtcCommitScene(_scene);
        auto t1 = std::chrono::high_resolution_clock::now();

        using namespace std::chrono_literals;
        std::cout << "Build acceleration structure for "
                  << vertices.size() << " vertices and "
                  << indices.size() << " faces in "
                  << (t1 - t0) / 1ns * 1e-6 << "ms"
                  << std::endl;
    }

    [[nodiscard]] static auto _trace_context(bool coherent) noexcept {
        static thread_local auto coherent_ctx = [] {
            RTCIntersectContext ctx{};
            rtcInitIntersectContext(&ctx);
            ctx.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;
            return ctx;
        }();
        static thread_local auto incoherent_ctx = [] {
            RTCIntersectContext ctx{};
            rtcInitIntersectContext(&ctx);
            ctx.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;
            return ctx;
        }();
        return coherent ? &coherent_ctx : &incoherent_ctx;
    }

public:
    [[nodiscard]] static auto build(std::span<const glm::vec3> vertices, std::span<const glm::uvec3> indices) noexcept {
        return std::unique_ptr<Accel>{new Accel{vertices, indices}};
    }

    ~Accel() noexcept {
        rtcReleaseScene(_scene);
        rtcReleaseDevice(_device);
    }

    [[nodiscard]] auto vertices() const noexcept { return _vertices; }
    [[nodiscard]] auto indices() const noexcept { return _indices; }

    void trace_closest(RayHit &ray_hit) const noexcept {
        rtcIntersect1(_scene, _trace_context(true), reinterpret_cast<RTCRayHit *>(&ray_hit));
    }

    void trace_any(Ray &ray) const noexcept {
        rtcOccluded1(_scene, _trace_context(true), reinterpret_cast<RTCRay *>(&ray));
    }

    void trace_closest(std::span<RayHit> rays, bool coherent) const noexcept {
        rtcIntersect1M(_scene, _trace_context(coherent), reinterpret_cast<RTCRayHit *>(rays.data()), rays.size(), sizeof(RayHit));
    }

    void trace_any(std::span<Ray> rays, bool coherent) const noexcept {
        rtcOccluded1M(_scene, _trace_context(coherent), reinterpret_cast<RTCRay *>(rays.data()), rays.size(), sizeof(Ray));
    }
};

std::unique_ptr<Mesh> Mesh::load(const std::filesystem::path &path) noexcept {

    auto absolute_path = std::filesystem::canonical(path).string();

    Assimp::Importer importer;
    importer.SetPropertyInteger(AI_CONFIG_PP_RVC_FLAGS,
                                aiComponent_ANIMATIONS
                                    | aiComponent_BONEWEIGHTS
                                    | aiComponent_CAMERAS
                                    | aiComponent_COLORS
                                    | aiComponent_LIGHTS
                                    | aiComponent_MATERIALS
                                    | aiComponent_NORMALS
                                    | aiComponent_TANGENTS_AND_BITANGENTS
                                    | aiComponent_TEXCOORDS
                                    | aiComponent_TEXTURES);
    auto model = importer.ReadFile(absolute_path.c_str(),
                                   aiProcess_Triangulate
                                       | aiProcess_JoinIdenticalVertices
                                       | aiProcess_RemoveComponent
                                       | aiProcess_ImproveCacheLocality
                                       | aiProcess_OptimizeMeshes
                                       | aiProcess_OptimizeGraph
                                       | aiProcess_DropNormals);

    if (model == nullptr
        || (model->mFlags & AI_SCENE_FLAGS_INCOMPLETE)
        || model->mRootNode == nullptr
        || model->mRootNode->mNumMeshes == 0) {
        std::cerr << "Failed to load scene." << std::endl;
        exit(-1);
    }

    auto mesh = model->mMeshes[0];
    std::vector<glm::vec3> vertices(mesh->mNumVertices);
    std::memcpy(vertices.data(), mesh->mVertices, vertices.size() * sizeof(glm::vec3));

    // normalize vertices coordinates
    auto bbox = BoundingBox::of(vertices);
    auto centroid = bbox.centroid();
    auto radius = std::accumulate(vertices.cbegin(), vertices.cend(), 0.0f, [centroid](auto r, auto p) noexcept {
        return std::max(r, glm::distance(p, centroid));
    });
    auto scale = 1.0f / radius;
    for (auto &&v : vertices) { v = (v - centroid) * scale; }

    std::vector<glm::uvec3> indices(mesh->mNumFaces);
    std::transform(mesh->mFaces, mesh->mFaces + mesh->mNumFaces, indices.begin(), [](auto &&face) noexcept {
        auto tri = face.mIndices;
        return glm::uvec3{tri[0], tri[1], tri[2]};
    });

    return std::unique_ptr<Mesh>{new Mesh{Accel::build(vertices, indices)}};
}

Mesh::Mesh(std::unique_ptr<Accel> accel) noexcept : _accel{std::move(accel)} {}
Mesh::~Mesh() noexcept = default;

std::span<const glm::vec3> Mesh::vertices() const noexcept { return _accel->vertices(); }
std::span<const glm::uvec3> Mesh::indices() const noexcept { return _accel->indices(); }

void Mesh::trace_closest(RayHit *ray) const noexcept { _accel->trace_closest(*ray); }
void Mesh::trace_any(Ray *ray) const noexcept { _accel->trace_any(*ray); }
void Mesh::trace_closest(std::span<RayHit> rays, bool coherent) const noexcept { _accel->trace_closest(rays, coherent); }
void Mesh::trace_any(std::span<Ray> rays, bool coherent) const noexcept { _accel->trace_any(rays, coherent); }

std::unique_ptr<Mesh> Mesh::build(std::span<const glm::vec3> vertices, std::span<const glm::uvec3> indices) noexcept {
    return std::unique_ptr<Mesh>{new Mesh{Accel::build(vertices, indices)}};
}
