//
// Created by Mike Smith on 2021/6/13.
//

#pragma once

#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

class Camera {

private:
    glm::mat3 _rotation_to_world{1.0f};
    float _distance{1.0f};
    float _pixel_scale{1.0f};
    float _z_plane{1.0f};
    
public:
    explicit Camera(glm::uvec2 resolution, float fov = 35.0f, float distance = 3.0f) noexcept;
    [[nodiscard]] glm::vec3 direction(glm::vec2 pixel) const noexcept;
    [[nodiscard]] glm::vec3 position() const noexcept;
    void set_resolution(glm::uvec2 res) noexcept;
    void set_fov(float fov) noexcept;
    void set_distance(float d) noexcept;
    void rotate_x(float angle) noexcept;
    void rotate_y(float angle) noexcept;
    void rotate_z(float angle) noexcept;
    [[nodiscard]] auto rotation_to_world() const noexcept { return _rotation_to_world; }
};
