//
// Created by Mike Smith on 2021/6/12.
//

#include <array>
#include <bitset>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include <core/volume.h>

[[nodiscard]] inline auto is_power_of_two(uint64_t x) noexcept {
    return x != 0u && (x & (x - 1u)) == 0u;
}

void check_resolution(uint32_t resolution) noexcept {
    if (resolution < 2 || resolution > 1024u || !is_power_of_two(resolution)) {
        std::cerr << "Invalid volume resolution: " << resolution << "." << std::endl;
        exit(-1);
    }
}

static const std::array cube_vertices{
    glm::vec3{0.0f, 0.0f, 0.0f},
    glm::vec3{1.0f, 0.0f, 0.0f},
    glm::vec3{1.0f, 1.0f, 0.0f},
    glm::vec3{0.0f, 1.0f, 0.0f},
    glm::vec3{0.0f, 1.0f, 1.0f},
    glm::vec3{1.0f, 1.0f, 1.0f},
    glm::vec3{1.0f, 0.0f, 1.0f},
    glm::vec3{0.0f, 0.0f, 1.0f}};

static const std::array cube_indices{
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

class Octree {

    class alignas(8) Node {

    public:
        static constexpr auto m000 = 0b00000001u;
        static constexpr auto m001 = 0b00000010u;
        static constexpr auto m011 = 0b00000100u;
        static constexpr auto m010 = 0b00001000u;
        static constexpr auto m110 = 0b00010000u;
        static constexpr auto m111 = 0b00100000u;
        static constexpr auto m101 = 0b01000000u;
        static constexpr auto m100 = 0b10000000u;

        inline static const glm::uvec3 d000{0u, 0u, 0u};
        inline static const glm::uvec3 d001{0u, 0u, 1u};
        inline static const glm::uvec3 d010{0u, 1u, 0u};
        inline static const glm::uvec3 d011{0u, 1u, 1u};
        inline static const glm::uvec3 d100{1u, 0u, 0u};
        inline static const glm::uvec3 d101{1u, 0u, 1u};
        inline static const glm::uvec3 d110{1u, 1u, 0u};
        inline static const glm::uvec3 d111{1u, 1u, 1u};

        static constexpr std::array m{m000, m001, m011, m010, m110, m111, m101, m100};
        inline static const std::array d{d000, d001, d011, d010, d110, d111, d101, d100};

    private:
        uint32_t _child_offset{};
        uint8_t _child_masks{};
        bool _full{};

    public:
        Node() noexcept = default;
        Node(bool full, uint32_t child_offset, uint8_t child_masks) noexcept
            : _full{full}, _child_offset{child_offset}, _child_masks{child_masks} {}
        [[nodiscard]] auto child_masks() const noexcept { return _child_masks; }
        [[nodiscard]] auto child_offset() const noexcept { return _child_offset; }
        [[nodiscard]] auto empty() const noexcept { return child_masks() == 0u; }
        [[nodiscard]] auto full() const noexcept { return _full; }
    };

private:
    std::vector<Node> _nodes;
    uint32_t _resolution;

private:
    explicit Octree(uint32_t resolution) noexcept
        : _resolution{resolution} { _nodes.reserve(65536u); }

    [[nodiscard]] auto _allocate_children() noexcept {
        auto index = static_cast<uint32_t>(_nodes.size());
        _nodes.emplace_back();
        return index;
    }

    [[nodiscard]] std::pair<bool, bool> /* (empty, full) */ _build(uint32_t node_index, std::span<Segment> segments, glm::uvec3 origin, uint32_t res) noexcept {

        // empty node
        if (segments.empty()) { return std::make_pair(true, false); }

        auto half_res = res / 2u;
        auto mid = origin + half_res;

        // split
        auto split_x = std::partition(segments.data(), segments.data() + segments.size(), [m = mid.x](auto s) noexcept { return s.x < m; });
        auto split_0y = std::partition(segments.data(), split_x, [m = mid.y](auto s) noexcept { return s.y < m; });
        auto split_1y = std::partition(split_x, segments.data() + segments.size(), [m = mid.y](auto s) noexcept { return s.y < m; });
        std::span<Segment> s00{segments.data(), split_0y};
        std::span<Segment> s01{split_0y, split_x};
        std::span<Segment> s10{split_x, split_1y};
        std::span<Segment> s11{split_1y, segments.data() + segments.size()};

        // leaf
        if (res == 2) {
            auto masks = 0u;
            auto p0 = [z = origin.z](auto s) noexcept { return std::any_of(s.begin(), s.end(), [z](auto s) { return z >= s.z_min && z <= s.z_max; }); };
            auto p1 = [z = origin.z + 1u](auto s) noexcept { return std::any_of(s.begin(), s.end(), [z](auto s) { return z >= s.z_min && z <= s.z_max; }); };
            if (p0(s00)) { masks |= Node::m000; }
            if (p1(s00)) { masks |= Node::m001; }
            if (p0(s01)) { masks |= Node::m010; }
            if (p1(s01)) { masks |= Node::m011; }
            if (p0(s10)) { masks |= Node::m100; }
            if (p1(s10)) { masks |= Node::m101; }
            if (p0(s11)) { masks |= Node::m110; }
            if (p1(s11)) { masks |= Node::m111; }
            auto full = masks == 0xffu;
            _nodes[node_index] = {full, 0u, static_cast<uint8_t>(masks)};
            return std::make_pair(false, full);
        }

        // inner node, index encoded in xyz order
        auto n000 = _allocate_children();
        auto n001 = _allocate_children();
        auto n011 = _allocate_children();
        auto n010 = _allocate_children();
        auto n110 = _allocate_children();
        auto n111 = _allocate_children();
        auto n101 = _allocate_children();
        auto n100 = _allocate_children();

        auto masks = 0u;
        auto full = true;

        auto p0 = [m = mid.z](auto s) noexcept {
            auto split = std::partition(s.data(), s.data() + s.size(), [m](auto s) noexcept { return s.z_min < m; });
            return std::span{s.data(), split};
        };

        auto p1 = [m = mid.z](auto s) noexcept {
            auto split = std::partition(s.data(), s.data() + s.size(), [m](auto s) noexcept { return s.z_max < m; });
            return std::span{split, s.data() + s.size()};
        };

        auto empty_full_000 = _build(n000, p0(s00), origin + Node::d000 * half_res, half_res);
        if (!empty_full_000.first) { masks |= Node::m000; }
        full &= empty_full_000.second;
        auto empty_full_001 = _build(n001, p1(s00), origin + Node::d001 * half_res, half_res);
        if (!empty_full_001.first) { masks |= Node::m001; }
        full &= empty_full_001.second;
        auto empty_full_011 = _build(n011, p1(s01), origin + Node::d011 * half_res, half_res);
        if (!empty_full_011.first) { masks |= Node::m011; }
        full &= empty_full_011.second;
        auto empty_full_010 = _build(n010, p0(s01), origin + Node::d010 * half_res, half_res);
        if (!empty_full_010.first) { masks |= Node::m010; }
        full &= empty_full_010.second;
        auto empty_full_110 = _build(n110, p0(s11), origin + Node::d110 * half_res, half_res);
        if (!empty_full_110.first) { masks |= Node::m110; }
        full &= empty_full_110.second;
        auto empty_full_111 = _build(n111, p1(s11), origin + Node::d111 * half_res, half_res);
        if (!empty_full_111.first) { masks |= Node::m111; }
        full &= empty_full_111.second;
        auto empty_full_101 = _build(n101, p1(s10), origin + Node::d101 * half_res, half_res);
        if (!empty_full_101.first) { masks |= Node::m101; }
        full &= empty_full_101.second;
        auto empty_full_100 = _build(n100, p0(s10), origin + Node::d100 * half_res, half_res);
        if (!empty_full_100.first) { masks |= Node::m100; }
        full &= empty_full_100.second;
        _nodes[node_index] = {full, n000 - node_index, static_cast<uint8_t>(masks)};
        return std::make_pair(masks == 0u, full);
    }

    void _to_mesh(uint32_t node_index, std::vector<glm::vec3> &vertices, std::vector<glm::uvec3> &indices, glm::uvec3 origin, uint32_t res) const noexcept {

        auto add_cube = [&vertices, &indices](glm::uvec3 origin, uint32_t res) noexcept {
            auto base_index = static_cast<uint32_t>(vertices.size());
            for (auto v : cube_vertices) { vertices.emplace_back(v * static_cast<float>(res) + glm::vec3{origin}); }
            for (auto i : cube_indices) { indices.emplace_back(i + base_index); }
        };

        auto node = _nodes[node_index];
        if (node.empty()) { return; }
        if (node.full()) { return add_cube(origin, res); }
        
        // leaf
        if (auto half_res = res / 2u; half_res == 1) {
            for (auto i = 0u; i < 8u; i++) {
                if (node.child_masks() & Node::m[i]) { add_cube(origin + Node::d[i] * half_res, half_res); }
            }
        } else {// inner nodes
            for (auto i = 0u; i < 8u; i++) {
                if (node.child_masks() & Node::m[i]) {
                    _to_mesh(node_index + node.child_offset() + i, vertices, indices, origin + Node::d[i] * half_res, half_res);
                }
            }
        }
    }
    
    [[nodiscard]] auto _trace_any(glm::vec3 origin, glm::vec3 direction, float t_max) const noexcept {
    
    }

public:
    [[nodiscard]] static auto build(std::vector<Segment> segments, uint32_t resolution) noexcept {
        check_resolution(resolution);
        std::unique_ptr<Octree> tree{new Octree{resolution}};
        auto root = tree->_allocate_children();
        auto t0 = std::chrono::high_resolution_clock::now();
        static_cast<void>(tree->_build(root, segments, glm::uvec3{}, resolution));
        auto t1 = std::chrono::high_resolution_clock::now();
        using namespace std::chrono_literals;
        std::cout << "Build octree with " << tree->_nodes.size() << " nodes in " << (t1 - t0) / 1ns * 1e-6 << " ms" << std::endl;
        return tree;
    }

    [[nodiscard]] auto to_mesh() const noexcept {
        std::vector<glm::vec3> vertices;
        std::vector<glm::uvec3> indices;
        _to_mesh(0u, vertices, indices, glm::uvec3{}, _resolution);
        return Mesh::build(vertices, indices);
    }

    [[nodiscard]] auto trace_closest(RayHit *ray) const noexcept {
    }

    [[nodiscard]] void trace_any(Ray &ray) const noexcept {
    }
};

[[nodiscard]] auto compute_segments(const Mesh &mesh, glm::mat4 M, uint32_t res) noexcept {

    std::vector<Segment> segments;
    segments.reserve(std::min(1024u * 1024u, res * res * 4u));

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
                            // ray exits the scene
                            if (ray_hit.hit.geom_id == Hit::invalid) {
                                //                                if (counter > 0) {// bad case, add surface hit & exit
                                //                                    auto z = glm::clamp((1.5f - t_start * 0.5f) * res_f - 0.5f, 0.0f, res_f - 1.0f);
                                //                                    Segment segment{static_cast<uint16_t>(x),
                                //                                                    static_cast<uint16_t>(y),
                                //                                                    static_cast<uint16_t>(z),
                                //                                                    static_cast<uint16_t>(z)};
                                //                                    segment_cache.emplace_back(segment);
                                //                                }
                                break;
                            }
                            // process the segment
                            if (auto c = glm::dot(ray_dir, ray_hit.hit.ng); c < 0.0f) {// front face
                                if (++counter == 1) { t_start = ray_hit.ray.t_max; }
                            } else if (c > 0.0f) {   // back face
                                if (--counter == 0) {// found segment
                                    auto z_range = (1.5f - glm::vec2{ray_hit.ray.t_max, t_start} * 0.5f) * res_f;
                                    if (auto z = glm::uvec2{glm::clamp(z_range + glm::vec2{0.5f, -0.5f}, 0.0f, res_f - 1.0f)}; z.x <= z.y) {
                                        Segment segment{static_cast<uint16_t>(x),
                                                        static_cast<uint16_t>(y),
                                                        static_cast<uint16_t>(z.x),
                                                        static_cast<uint16_t>(z.y)};
                                        segment_cache.emplace_back(segment);
                                    }
                                } else if (counter < 0) {// bad case, add surface voxel & reset
                                                         //                                    auto z = glm::clamp((1.5f - ray_hit.ray.t_max * 0.5f) * res_f + 0.5f, 0.0f, res_f - 1.0f);
                                                         //                                    Segment segment{static_cast<uint16_t>(x),
                                                         //                                                    static_cast<uint16_t>(y),
                                                         //                                                    static_cast<uint16_t>(z),
                                                         //                                                    static_cast<uint16_t>(z)};
                                                         //                                    segment_cache.emplace_back(segment);
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

    check_resolution(resolution);

    auto t0 = std::chrono::high_resolution_clock::now();
    auto segments = compute_segments(mesh, camera_to_world, resolution);
    auto t1 = std::chrono::high_resolution_clock::now();
    using namespace std::chrono_literals;
    std::cout << "Found " << segments.size() << " segments in " << (t1 - t0) / 1ns * 1e-6 << " ms" << std::endl;

    auto volume = std::unique_ptr<Volume>{new Volume{Octree::build(segments, resolution), resolution}};
    volume->_mesh = volume->_octree->to_mesh();
    return volume;
}

Volume::Volume(std::unique_ptr<Octree> octree, size_t resolution) noexcept
    : _octree{std::move(octree)}, _resolution{resolution} {}

Volume::~Volume() noexcept = default;
