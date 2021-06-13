//
// Created by Mike Smith on 2021/6/12.
//

#pragma once

#include <filesystem>
#include <memory>

#include <core/bbox.h>
#include <core/ray.h>

class Accel;

class Mesh {

private:
    std::unique_ptr<Accel> _accel;
    explicit Mesh(std::unique_ptr<Accel> accel) noexcept;

public:
    ~Mesh() noexcept;
    [[nodiscard]] static std::unique_ptr<Mesh> load(const std::filesystem::path &path) noexcept;
    [[nodiscard]] static std::unique_ptr<Mesh> build(std::span<const glm::vec3> vertices, std::span<const glm::uvec3> indices) noexcept;
    [[nodiscard]] std::span<const glm::vec3> vertices() const noexcept;
    [[nodiscard]] std::span<const glm::uvec3> indices() const noexcept;
    void trace_closest(RayHit *ray) const noexcept;
    void trace_any(Ray *ray) const noexcept;
    void trace_closest(std::span<RayHit> rays, bool coherent = true) const noexcept;
    void trace_any(std::span<Ray> rays, bool coherent = true) const noexcept;
};
