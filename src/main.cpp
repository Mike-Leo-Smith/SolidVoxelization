#include <iostream>
#include <vector>

#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/scene.h>
#include <embree3/rtcore.h>
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/string_cast.hpp>
#include <opencv2/opencv.hpp>

struct BoundingBox {

    glm::vec3 min;
    glm::vec3 max;

    BoundingBox(glm::vec3 min, glm::vec3 max) noexcept
        : min{min}, max{max} {}

    [[nodiscard]] auto valid() const noexcept { return glm::all(glm::lessThanEqual(min, max)); }

    void update(glm::vec3 p) noexcept {
        min = glm::min(min, p);
        max = glm::max(max, p);
    }

    [[nodiscard]] auto centroid() const noexcept { return 0.5f * (min + max); }
    [[nodiscard]] auto extent() const noexcept { return max - min; }
};

struct alignas(16) Ray {
    glm::vec3 o;
    float t_min;
    glm::vec3 d;
    float time;
    float t_max;
    unsigned int mask;
    unsigned int id;
    unsigned int flags;
};

struct alignas(16) Hit {
    glm::vec3 ng;
    glm::vec2 uv;
    unsigned int prim_id;
    unsigned int geom_id;
    unsigned int inst_id;
};

/* Combined ray/hit structure for a single ray */
struct RayHit {
    Ray ray;
    Hit hit;
};

int main() {

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
    auto model = importer.ReadFile("data/bunny.ply",
                                   aiProcess_Triangulate
                                       | aiProcess_JoinIdenticalVertices
                                       | aiProcess_RemoveComponent
                                       | aiProcess_PreTransformVertices
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

    RTCDevice device = rtcNewDevice(nullptr);
    RTCScene scene = rtcNewScene(device);
    rtcSetSceneBuildQuality(scene, RTC_BUILD_QUALITY_HIGH);

    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    auto vb = static_cast<glm::vec3 *>(rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_VERTEX, 0,
        RTC_FORMAT_FLOAT3, sizeof(glm::vec3), mesh->mNumVertices));
    std::memcpy(vb, mesh->mVertices, sizeof(glm::vec3) * mesh->mNumVertices);
    auto ib = static_cast<glm::uvec3 *>(rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_INDEX, 0,
        RTC_FORMAT_UINT3, sizeof(glm::uvec3), mesh->mNumFaces));
    for (auto i = 0u; i < mesh->mNumFaces; i++) {
        auto indices = mesh->mFaces[i].mIndices;
        ib[i] = {indices[0], indices[1], indices[2]};
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    BoundingBox bbox{vb[0], vb[0]};
    for (auto i = 1u; i < mesh->mNumVertices; i++) { bbox.update(vb[i]); }

    std::cout << "Vertices: " << mesh->mNumVertices << "\n"
              << "Faces:    " << mesh->mNumFaces << "\n"
              << "BBox:     " << glm::to_string(bbox.min) << " -> " << glm::to_string(bbox.max) << std::endl;

    auto max_extent = glm::compMax(bbox.extent());
    auto centroid = bbox.centroid();

    std::vector<RayHit> rays;
    rays.resize(1024 * 1024);
    for (auto y = 0u; y < 1024u; y++) {
        for (auto x = 0u; x < 1024u; x++) {
            auto index = y * 1024u + x;
            auto &&r = rays[index];
            auto dx = (static_cast<float>(x) / 1024.0f - 0.5f) * max_extent;
            auto dy = (0.5f - static_cast<float>(y) / 1024.0f) * max_extent;
            auto dz = max_extent;
            r.ray.o = centroid + glm::vec3{dx, dy, dz};
            r.ray.t_min = 0.0f;
            r.ray.t_max = std::numeric_limits<float>::max();
            r.ray.d = {0.0f, 0.0f, -1.0f};
            r.hit.geom_id = RTC_INVALID_GEOMETRY_ID;
        }
    }

    RTCIntersectContext context{};
    rtcInitIntersectContext(&context);
    context.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;

    auto t0 = std::chrono::high_resolution_clock::now();
    rtcIntersect1M(scene, &context, reinterpret_cast<RTCRayHit *>(rays.data()), 1024u * 1024u, sizeof(RayHit));
    auto t1 = std::chrono::high_resolution_clock::now();

    using namespace std::chrono_literals;
    std::cout << "Rendering TIme: " << (t1 - t0) / 1ns * 1e-6 << " ms" << std::endl;

    cv::Mat image{1024, 1024, CV_32FC3, cv::Scalar::all(0)};
    for (auto i = 0u; i < 1024U * 1024u; i++) {
        auto &&pixel = reinterpret_cast<glm::vec3 *>(image.data)[i];
        auto hit = rays[i].hit;
        if (hit.geom_id != RTC_INVALID_GEOMETRY_ID) {
            auto ng = glm::normalize(hit.ng);
            pixel = ng * 0.5f + 0.5f;
        }
    }
    
    cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
    cv::imwrite("normal.exr", image);

    rtcReleaseScene(scene);
    rtcReleaseDevice(device);

    return 0;
}
