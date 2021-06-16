//
// Created by Mike on 6/16/2021.
//

#pragma once

#include <memory>
#include <gpu/common.h>

class Volume;

class CUDAOctree {

    friend class Volume;

public:
    class alignas(8) Node {

    private:
        uint32_t _child_offset{};
        uint8_t _child_masks{};
        bool _full{};

    public:
        [[nodiscard]] SV_XPU auto child_masks() const noexcept { return _child_masks; }
        [[nodiscard]] SV_XPU auto child_offset() const noexcept { return _child_offset; }
        [[nodiscard]] SV_XPU auto empty() const noexcept { return child_masks() == 0u; }
        [[nodiscard]] SV_XPU auto full() const noexcept { return _full; }
    };

private:
    std::unique_ptr<CUDABuffer<Node>> _nodes{nullptr};
    uint32_t _resolution{0u};

    [[nodiscard]] static auto create(const void *nodes, size_t count, uint32_t res) noexcept {
        auto tree = std::make_unique<CUDAOctree>();
        tree->_nodes = CUDABuffer<Node>::create(count);
        tree->_nodes->upload(reinterpret_cast<const Node *>(nodes));
        tree->_resolution = res;
        return tree;
    }

public:
    void trace_closest(const Ray *rays, Hit *hits, uint32_t w, uint32_t h) const noexcept;
};
