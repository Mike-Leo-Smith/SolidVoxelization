//
// Created by Mike Smith on 2021/6/12.
//

#pragma once

#include <core/mesh.h>
#include <core/segment.h>

class Octree;

class Volume {

private:
    std::unique_ptr<Octree> _octree;
    size_t _resolution;

private:
    Volume(std::unique_ptr<Octree> octree, size_t resolution) noexcept;

public:
    ~Volume() noexcept;
    [[nodiscard]] static std::unique_ptr<Volume> from(
        const Mesh &mesh, size_t resolution,
        glm::mat4 camera_to_world = glm::identity<glm::mat4>()) noexcept;
    [[nodiscard]] auto resolution() const noexcept { return _resolution; }
};
