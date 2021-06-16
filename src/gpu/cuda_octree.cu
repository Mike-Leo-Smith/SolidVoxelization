//
// Created by Mike on 6/16/2021.
//

#include <glm/gtx/component_wise.hpp>
#include <gpu/cuda_octree.h>

[[nodiscard]] SV_XPU inline auto _intersect_box(Ray ray, glm::vec3 bbox_min, float bbox_r) noexcept {
    auto bbox_max = bbox_min + bbox_r;
    auto t_min = (bbox_min - ray.o) / ray.d;
    auto t_max = (bbox_max - ray.o) / ray.d;
    auto o_mat = glm::mat3{ray.o, ray.o, ray.o};
    glm::vec3 p_min[]{ray.o + t_min.x * ray.d, ray.o + t_min.y * ray.d, ray.o + t_min.z * ray.d};
    glm::vec3 p_max[]{ray.o + t_max.x * ray.d, ray.o + t_max.y * ray.d, ray.o + t_max.z * ray.d};
    auto valid = [ray, bbox_min, bbox_max](auto t, auto p) noexcept {
        return glm::not_(glm::isnan(t))
               && glm::greaterThanEqual(t, glm::vec3{ray.t_min})
               && glm::lessThanEqual(t, glm::vec3{ray.t_max})
               && glm::bvec3{p[0].y >= bbox_min.y && p[0].y <= bbox_max.y
                                 && p[0].z >= bbox_min.z && p[0].z <= bbox_max.z,
                             p[1].z >= bbox_min.z && p[1].z <= bbox_max.z
                                 && p[1].x >= bbox_min.x && p[1].x <= bbox_max.x,
                             p[2].x >= bbox_min.x && p[2].x <= bbox_max.x
                                 && p[2].y >= bbox_min.y && p[2].y <= bbox_max.y};
    };
    auto t_invalid = glm::vec3{std::numeric_limits<float>::max()};
    auto valid_min = valid(t_min, p_min);
    auto valid_max = valid(t_max, p_max);
    t_min = glm::mix(t_invalid, t_min, valid_min);
    t_max = glm::mix(t_invalid, t_max, valid_max);
    auto t = glm::min(t_min, t_max);
    auto is_min = glm::lessThanEqual(t_min, t_max);
    Hit hit{};
    hit.t = std::numeric_limits<float>::max();
    hit.valid = false;
    if (t.x < hit.t) {
        hit.p = is_min.x ? p_min[0] : p_max[0];
        hit.t = t.x;
        hit.ng = {is_min.x ? -1.0f : 1.0f, 0.0f, 0.0f};
        hit.valid = valid_min.x || valid_max.x;
    }
    if (t.y < hit.t) {
        hit.p = is_min.y ? p_min[1] : p_max[1];
        hit.t = t.y;
        hit.ng = {0.0f, is_min.y ? -1.0f : 1.0f, 0.0f};
        hit.valid = valid_min.y || valid_max.y;
    }
    if (t.z < hit.t) {
        hit.p = is_min.z ? p_min[2] : p_max[2];
        hit.t = t.z;
        hit.ng = {0.0f, 0.0f, is_min.z ? -1.0f : 1.0f};
        hit.valid = valid_min.z || valid_max.z;
    }
    return hit;
}

static constexpr auto m000 = 0b00000001u;
static constexpr auto m001 = 0b00000010u;
static constexpr auto m011 = 0b00000100u;
static constexpr auto m010 = 0b00001000u;
static constexpr auto m110 = 0b00010000u;
static constexpr auto m111 = 0b00100000u;
static constexpr auto m101 = 0b01000000u;
static constexpr auto m100 = 0b10000000u;
__constant__ uint32_t _node_m[]{m000, m001, m011, m010, m110, m111, m101, m100};

static constexpr auto d000 = float3{0.0f, 0.0f, 0.0f};
static constexpr auto d001 = float3{0.0f, 0.0f, 1.0f};
static constexpr auto d010 = float3{0.0f, 1.0f, 0.0f};
static constexpr auto d011 = float3{0.0f, 1.0f, 1.0f};
static constexpr auto d100 = float3{1.0f, 0.0f, 0.0f};
static constexpr auto d101 = float3{1.0f, 0.0f, 1.0f};
static constexpr auto d110 = float3{1.0f, 1.0f, 0.0f};
static constexpr auto d111 = float3{1.0f, 1.0f, 1.0f};
__constant__ float3 _node_d[]{d000, d001, d011, d010, d110, d111, d101, d100};

__global__ void octree_trace_closest(const CUDAOctree::Node *nodes, uint32_t resolution, const Ray *rays, Hit *hits, uint32_t w, uint32_t h) {

    auto x = threadIdx.x + blockIdx.x * blockDim.x;
    auto y = threadIdx.y + blockIdx.y * blockDim.y;
    if (x >= w || y >= h) { return; }

    auto tid = x + y * w;

    Hit closest{};
    closest.t = std::numeric_limits<float>::max();
    closest.valid = false;
    if (nodes[0].empty()) {
        hits[tid] = closest;
        return;
    }
    
    struct alignas(16) TraceContext {
        glm::vec3 o;
        float r;
        uint32_t index;
        float t;
    };

    static constexpr auto stack_size = 32u;
    TraceContext stack[stack_size];
    auto sp = 0u;

    auto ray = rays[tid];

    auto add_node = [&sp, ray, &stack, &closest, nodes](auto index, auto o, auto r) noexcept {
        if (auto hit = _intersect_box(ray, o, r); hit.valid && hit.t < closest.t) {
            auto node = nodes[index];
            if (node.full()) {
                closest = hit;
            } else if (r == 2.0f) {
                #pragma unroll
                for (auto i = 0u; i < 8u; i++) {
                    if ((node.child_masks() & _node_m[i])) {
                        auto d = _node_d[i];
                        if (auto child_hit = _intersect_box(ray, o + glm::vec3{d.x, d.y, d.z}, 1.0f);
                            child_hit.valid && child_hit.t < closest.t) {
                            closest = child_hit;
                        }
                    }
                }
            } else {
                if (sp == stack_size) { printf("warning: stack overflows\n"); }
                stack[sp++] = {o, r, index, hit.t};
            }
        }
    };

    add_node(0u, glm::vec3{}, static_cast<float>(resolution));
    while (sp != 0u) {
        auto ctx = stack[--sp];
        if (ctx.t >= closest.t) { continue; }
        auto node = nodes[ctx.index];
        auto half_r = ctx.r * 0.5f;
        #pragma unroll
        for (auto i = 0u; i < 8u; i++) {
            if (node.child_masks() & _node_m[i]) {
                auto d = _node_d[i];
                add_node(ctx.index + node.child_offset() + i, ctx.o + glm::vec3{d.x, d.y, d.z} * half_r, half_r);
            }
        }
    }

    hits[tid] = closest;
}

void CUDAOctree::trace_closest(const Ray *rays, Hit *hits, uint32_t width, uint32_t height) const noexcept {
    static constexpr auto block_size = 16u;
    auto blocks_x = (width + block_size - 1u) / block_size;
    auto blocks_y = (height + block_size - 1u) / block_size;
    octree_trace_closest<<<dim3(blocks_x, blocks_y), dim3(block_size, block_size)>>>(_nodes->data(), _resolution, rays, hits, width, height);
}
