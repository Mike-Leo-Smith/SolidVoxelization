//
// Created by Mike Smith on 2021/6/12.
//

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>
#include <mutex>

#include <core/volume.h>

class Octree {
};

[[nodiscard]] inline auto is_power_of_two(uint64_t x) noexcept {
    return x != 0u && (x & (x - 1u)) == 0u;
}

[[nodiscard]] auto compute_segments(const Mesh &mesh, glm::mat4 M, uint32_t res) noexcept {

    std::vector<Segment> segments;
    segments.reserve(std::min(1024u * 1024u * 4u, res * res * 4u));

    std::mutex mutex;
    std::atomic_uint curr_block{0u};
    auto num_workers = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;
    workers.reserve(num_workers);
    for (auto i = 0u; i < num_workers; i++) {
        workers.emplace_back([&mesh, M, res, &mutex, &curr_block, &segments] {
            auto ray_dir = glm::mat3{M} * glm::vec3{0.0f, 0.0f, -1.0f};
            auto res_f = static_cast<float>(res);
            auto scale = 2.0f / res_f;
            auto eps = std::max(scale / 64.0f, 1e-5f);

            auto block_size = std::min(res, 16u);
            auto blocks_per_edge = (res + block_size - 1u) / block_size;
            auto num_blocks = blocks_per_edge * blocks_per_edge;

            std::vector<Segment> segment_cache;
            segment_cache.reserve(block_size * block_size * 16u);
            for (auto block = curr_block++; block < num_blocks; block = curr_block++) {
                segment_cache.clear();
                auto x_offset = block % blocks_per_edge * block_size;
                auto y_offset = block / blocks_per_edge * block_size;
                for (auto y = y_offset; y < y_offset + block_size; y++) {
                    for (auto x = x_offset; x < x_offset + block_size; x++) {
                        RayHit ray_hit{};
                        auto dx = (static_cast<float>(x) + 0.5f) * scale - 1.0f;
                        auto dy = (static_cast<float>(y) + 0.5f) * scale - 1.0f;
                        auto dz = 2.0f;
                        ray_hit.ray.o = glm::vec3{M * glm::vec4{dx, dy, dz, 1.0f}};
                        ray_hit.ray.d = ray_dir;

                        // find hits along the ray...
                        float t_start;
                        auto counter = 0;
                        for (;;) {
                            ray_hit.ray.t_min = ray_hit.ray.t_max + eps;
                            ray_hit.ray.t_max = std::numeric_limits<float>::infinity();
                            ray_hit.hit.geom_id = Hit::invalid;
                            mesh.trace_closest(&ray_hit);
                            if (ray_hit.hit.geom_id == Hit::invalid) { break; }
                            // process the segment
                            if (auto c = glm::dot(ray_dir, ray_hit.hit.ng); c < 0.0f) {// front face
                                if (++counter == 1) { t_start = ray_hit.ray.t_max; }
                            } else if (c > 0.0f) {
                                if (--counter == 0) {// found segment
                                    auto z_range = (1.5f - glm::vec2{ray_hit.ray.t_max, t_start} * 0.5f) * res_f;
                                    if (auto z = glm::uvec2{glm::clamp(z_range + glm::vec2{0.5f, -0.5f}, 0.0f, res_f - 1.0f)}; z.x <= z.y) {
                                        Segment segment{static_cast<uint16_t>(x),
                                                        static_cast<uint16_t>(y),
                                                        static_cast<uint16_t>(z.x),
                                                        static_cast<uint16_t>(z.y)};
                                        segment_cache.emplace_back(segment);
                                    }
                                } else if (counter < 0) {// bad case, reset
                                    counter = 0;
                                }
                            }
                        }
                    }
                }

                // flush...
                if (!segment_cache.empty()) {
                    std::scoped_lock lock{mutex};
                    std::copy(segment_cache.cbegin(), segment_cache.cend(), std::back_inserter(segments));
                }
            }
        });
    }
    for (auto &&w : workers) { w.join(); }
    return segments;
}

std::unique_ptr<Volume> Volume::from(const Mesh &mesh, size_t resolution, glm::mat4 camera_to_world) noexcept {

    if (resolution > 1024u || !is_power_of_two(resolution)) {
        std::cerr << "Invalid volume resolution: " << resolution << "." << std::endl;
        exit(-1);
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    auto segments = compute_segments(mesh, camera_to_world, resolution);
    auto t1 = std::chrono::high_resolution_clock::now();

    using namespace std::chrono_literals;
    std::cout << "Computed " << segments.size() << " segments in " << (t1 - t0) / 1ns * 1e-6 << " ms" << std::endl;

    // save test data
    {
        std::ofstream file{"test-" + std::to_string(resolution) + ".txt"};
        file << resolution << "\n"
             << segments.size() << "\n";
        for (auto s : segments) {
            file << s.x << " " << s.y << " " << s.z_min << " " << s.z_max << "\n";
        }
    }

    // TODO: construct octree...
    const std::array cube_vertices{
        glm::vec3{0.0f, 0.0f, 0.0f},
        glm::vec3{1.0f, 0.0f, 0.0f},
        glm::vec3{1.0f, 1.0f, 0.0f},
        glm::vec3{0.0f, 1.0f, 0.0f},
        glm::vec3{0.0f, 1.0f, 1.0f},
        glm::vec3{1.0f, 1.0f, 1.0f},
        glm::vec3{1.0f, 0.0f, 1.0f},
        glm::vec3{0.0f, 0.0f, 1.0f}};

    const std::array cube_indices{
        glm::uvec3{0u, 2u, 1u},
        glm::uvec3{0u, 3u, 2u},
        glm::uvec3{0u, 1u, 6u},
        glm::uvec3{0u, 6u, 7u},
        glm::uvec3{1u, 2u, 5u},
        glm::uvec3{1u, 5u, 6u},
        glm::uvec3{3u, 4u, 5u},
        glm::uvec3{3u, 5u, 2u},
        glm::uvec3{0u, 4u, 3u},
        glm::uvec3{0u, 7u, 4u},
        glm::uvec3{4u, 7u, 6u},
        glm::uvec3{4u, 6u, 5u}};

    // save segments as .obj
    std::vector<glm::vec3> vertices;
    std::vector<glm::uvec3> indices;
    std::for_each(segments.cbegin(), segments.cend(), [&](Segment s) noexcept {
        auto e = glm::vec3{1.0f, 1.0f, s.z_max - s.z_min + 1.0f};
        auto t = glm::vec3{s.x, s.y, s.z_min};
        auto dy = glm::vec3{0.0f, 1.0f, 0.0f};
//        auto scale = 2.0f / static_cast<float>(resolution);
        auto index_base = static_cast<uint32_t>(vertices.size());
//        for (auto v : cube_vertices) { vertices.emplace_back((v * e + t) * scale - 1.0f + dy); }
        for (auto v : cube_vertices) { vertices.emplace_back(v * e + t); }
        for (auto f : cube_indices) { indices.emplace_back(f + index_base); }
    });

//    std::ofstream file{"test.obj"};
//    file << std::setprecision(10);
//    for (auto v : vertices) { file << "v " << v.x << " " << v.y << " " << v.z << "\n"; }
//    for (auto i : indices) { file << "f " << i.x + 1u << " " << i.y + 1u << " " << i.z + 1u << "\n"; }

    auto volume = std::unique_ptr<Volume>{new Volume{nullptr /* TODO */, resolution}};
    volume->_mesh = Mesh::build(vertices, indices);
    return volume;
}

Volume::Volume(std::unique_ptr<Octree> octree, size_t resolution) noexcept
    : _octree{std::move(octree)}, _resolution{resolution} {}

Volume::~Volume() noexcept = default;
