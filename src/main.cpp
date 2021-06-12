#include <iostream>
#include <vector>
#include <numeric>

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
    [[nodiscard]] auto diagonal() const noexcept { return glm::length(extent()); }
    [[nodiscard]] auto radius() const noexcept { return 0.5f * diagonal(); }
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

struct Segment {
    int x;
    int y;
    int z_min;
    int z_max;
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
    auto model = importer.ReadFile("data/bunny.obj",
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

    static constexpr auto diameter = 1024;

    BoundingBox bbox{vb[0], vb[0]};
    for (auto i = 1u; i < mesh->mNumVertices; i++) { bbox.update(vb[i]); }
    auto centroid = bbox.centroid();
    auto radius = std::accumulate(vb, vb + mesh->mNumVertices, 0.0f, [centroid](auto r, auto p) noexcept {
        return std::max(r, glm::distance(p, centroid));
    });
    auto scale = static_cast<float>(diameter) * 0.5f / radius;
    for (auto i = 0u; i < mesh->mNumVertices; i++) {
        vb[i] = (vb[i] - centroid) * scale;
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    std::cout << "Vertices: " << mesh->mNumVertices << "\n"
              << "Faces:    " << mesh->mNumFaces << "\n"
              << "BBox:     " << glm::to_string(bbox.min) << " -> " << glm::to_string(bbox.max) << std::endl;

    RTCIntersectContext context{};
    rtcInitIntersectContext(&context);
    context.flags = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;

    cv::Mat image{diameter, diameter, CV_32FC3, cv::Scalar::all(0)};
    std::vector<RayHit> rays;
    rays.resize(diameter * diameter);
    glm::mat4 R{1.0f};
    for (;;) {
        auto omega = glm::mat3{R} * glm::vec3{0.0f, 0.0f, -1.0f};
        for (auto y = 0; y < diameter; y++) {
            for (auto x = 0; x < diameter; x++) {
                auto index = y * diameter + x;
                auto &&r = rays[index];
                auto dx = static_cast<float>(x - diameter / 2) + 0.5f;
                auto dy = static_cast<float>(diameter / 2 - y) - 0.5f;
                auto dz = static_cast<float>(diameter);
                auto d = glm::vec3{R * glm::vec4{dx, dy, dz, 1.0f}};
                r.ray.o = centroid + d;
                r.ray.t_min = 0.0f;
                r.ray.t_max = std::numeric_limits<float>::max();
                r.ray.d = omega;
                r.hit.geom_id = RTC_INVALID_GEOMETRY_ID;
            }
        }

        auto t0 = std::chrono::high_resolution_clock::now();
        rtcIntersect1M(scene, &context, reinterpret_cast<RTCRayHit *>(rays.data()), diameter * diameter, sizeof(RayHit));
        auto t1 = std::chrono::high_resolution_clock::now();

        using namespace std::chrono_literals;
        std::cout << "Rendering TIme: " << (t1 - t0) / 1ns * 1e-6 << " ms" << std::endl;
        for (auto i = 0u; i < 1024U * 1024u; i++) {
            auto &&pixel = reinterpret_cast<glm::vec3 *>(image.data)[i];
            auto hit = rays[i].hit;
            if (hit.geom_id != RTC_INVALID_GEOMETRY_ID) {
                auto ng = glm::normalize(hit.ng);
                ng = glm::vec3{ng.z, ng.y, ng.x};
                pixel = ng * 0.5f + 0.5f;
            } else {
                pixel = {};
            }
        }
        cv::imshow("Display", image);
        auto key = cv::waitKey(1);
        if (key == 'q') { break; }
        if (key == 'a') {
            R = glm::rotate(R, glm::radians(5.0f), glm::vec3{0.0f, 1.0f, 0.0f});
        } else if (key == 'd') {
            R = glm::rotate(R, glm::radians(-5.0f), glm::vec3{0.0f, 1.0f, 0.0f});
        } else if (key == 'w') {
            R = glm::rotate(R, glm::radians(5.0f), glm::vec3{1.0f, 0.0f, 0.0f});
        } else if (key == 's') {
            R = glm::rotate(R, glm::radians(-5.0f), glm::vec3{1.0f, 0.0f, 0.0f});
        }
    }
    cv::imwrite("normal.exr", image);

    rtcReleaseScene(scene);
    rtcReleaseDevice(device);

    return 0;
}
