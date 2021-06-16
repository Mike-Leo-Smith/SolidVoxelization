#include <gpu/cuda_camera.h>

template<uint32_t times>
[[nodiscard]] SV_XPU auto tea(glm::uvec2 xy) noexcept {
    auto v0 = xy.x;
    auto v1 = xy.y;
    auto s0 = 0u;
    for (auto n = 0u; n < times; n++) {
        s0 += 0x9e3779b9u;
        v0 += ((v1 << 4u) + 0xa341316cu) ^ (v1 + s0) ^ ((v1 >> 5u) + 0xc8013ea4u);
        v1 += ((v0 << 4u) + 0xad90777du) ^ (v0 + s0) ^ ((v0 >> 5u) + 0x7e95761eu);
    }
    return v0;
};

template<uint32_t base>
[[nodiscard]] SV_XPU auto halton(uint32_t i) noexcept {
    auto f = 1.0f;
    auto inv_base = 1.0f / base;
    auto r = 0.0f;
    while (i > 0u) {
        f = f * inv_base;
        r = r + f * (i % base);
        i /= base;
    }
    return r;
}

[[nodiscard]] SV_XPU auto pixel_sample(glm::uvec2 p, uint32_t frame) noexcept {
    auto index = tea<4u>(p) + frame;
    auto x = halton<2>(index);
    auto y = halton<3>(index);
    return glm::vec2{x, y};
}

__global__ void generate_rays(Ray *rays, uint32_t w, uint32_t h, uint32_t frame, glm::vec3 o, glm::mat3 m, float pixel_scale, float z_plane) {
    auto x = threadIdx.x + blockIdx.x * blockDim.x;
    auto y = threadIdx.y + blockIdx.y * blockDim.y;
    if (x >= w || y >= h) { return; }
    auto tid = y * w + x;
    auto xy = glm::uvec2{x, y};
    auto p = pixel_sample(xy, frame) + glm::vec2{xy};
    auto d = glm::normalize(m * glm::vec3{p.x * pixel_scale - 1.0f, 1.0f - p.y * pixel_scale, -z_plane});
    rays[tid] = {o, 0.0f, d, std::numeric_limits<float>::max()};
}

void cuda_camera_generate_rays(
    Camera cam, glm::mat4 world_to_object,
    Ray *rays, uint32_t w, uint32_t h, uint32_t frame) noexcept {
    auto o = glm::vec3{world_to_object * glm::vec4{cam.position(), 1.0f}};
    auto m = glm::mat3{world_to_object} * cam.rotation_to_world();
    auto pixel_scale = cam.pixel_scale();
    auto z_plane = cam.z_plane();
    static constexpr auto block_size = 16u;
    auto bx = (w + block_size - 1u) / block_size;
    auto by = (h + block_size - 1u) / block_size;
    generate_rays<<<dim3(bx, by), dim3(block_size, block_size)>>>(rays, w, h, frame, o, m, pixel_scale, z_plane);
}