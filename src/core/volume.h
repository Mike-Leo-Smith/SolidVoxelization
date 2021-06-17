//
// Created by Mike Smith on 2021/6/12.
//

#pragma once

#include <memory>
#include <core/mesh.h>
#include <core/segment.h>

#ifdef SV_CUDA_AVAILABLE
#include <gpu/cuda_octree.h>
#endif

class Octree;

class Volume {

private:
    std::unique_ptr<Octree> _octree;
    std::unique_ptr<Mesh> _mesh;
    size_t _resolution;

private:
    Volume(std::unique_ptr<Octree> octree, size_t resolution) noexcept;

public:
    ~Volume() noexcept;
    [[nodiscard]] static std::unique_ptr<Volume> from(
        const Mesh &mesh, size_t resolution,
        glm::mat4 camera_to_world = glm::identity<glm::mat4>()) noexcept;
    [[nodiscard]] auto resolution() const noexcept { return _resolution; }
    [[nodiscard]] bool trace_any(Ray ray) const noexcept;
    [[nodiscard]] Hit trace_closest(Ray ray) const noexcept;
    void dump(const std::filesystem::path &path) const noexcept;

#ifdef SV_CUDA_AVAILABLE
    [[nodiscard]] std::unique_ptr<CUDAOctree> cuda() const noexcept;
#endif

};
