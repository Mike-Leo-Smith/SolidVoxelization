#include <iostream>

#include <assimp/scene.h>
#include <glm/glm.hpp>
#include <embree3/rtcore.h>

int main() {

    RTCDevice device = rtcNewDevice(nullptr);
    RTCScene scene = rtcNewScene(device);
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    auto vb = static_cast<glm::vec3 *>(rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_VERTEX, 0,
        RTC_FORMAT_FLOAT3, sizeof(glm::vec3), 3));
    vb[0] = {0.0f, 0.0f, 0.0f};
    vb[1] = {1.0f, 0.0f, 0.0f};
    vb[2] = {0.0f, 1.0f, 0.0f};

    auto ib = static_cast<glm::uvec3 *>(rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_INDEX, 0,
        RTC_FORMAT_UINT3, sizeof(glm::uvec3), 1));
    ib[0] = {0u, 1u, 2u};

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    RTCRayHit rayhit{};
    rayhit.ray.org_x = 0.f;
    rayhit.ray.org_y = 0.f;
    rayhit.ray.org_z = -1.f;
    rayhit.ray.dir_x = 0.f;
    rayhit.ray.dir_y = 0.f;
    rayhit.ray.dir_z = 1.f;
    rayhit.ray.tnear = 0.f;
    rayhit.ray.tfar = std::numeric_limits<float>::infinity();
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    RTCIntersectContext context{};
    rtcInitIntersectContext(&context);

    rtcIntersect1(scene, &context, &rayhit);

    if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
        std::cout << "Intersection at t = " << rayhit.ray.tfar << "\n"
                  << "Normal = (" << rayhit.hit.Ng_x << ", " << rayhit.hit.Ng_y << ", " << rayhit.hit.Ng_z << ")"
                  << std::endl;
    } else {
        std::cout << "No Intersection" << std::endl;
    }

    rtcReleaseScene(scene);
    rtcReleaseDevice(device);

    return 0;
}
