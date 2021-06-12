//
// Created by Mike Smith on 2021/6/12.
//

#pragma once

#include <span>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

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
    
    [[nodiscard]] static auto of(std::span<const glm::vec3> vertices) noexcept {
        BoundingBox bbox{vertices.front(), vertices.front()};
        for (auto v : vertices.subspan(1)) { bbox.update(v); }
        return bbox;
    }
};
