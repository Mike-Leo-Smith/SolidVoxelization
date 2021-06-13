//
// Created by Mike Smith on 2021/6/13.
//

#include <core/camera.h>

glm::vec3 Camera::direction(glm::vec2 pixel) const noexcept {
    glm::vec3 w{pixel.x * _pixel_scale - 1.0f, 1.0f - pixel.y * _pixel_scale, -_z_plane};
    return glm::normalize(_rotation_to_world * w);
}

Camera::Camera(glm::uvec2 resolution, float fov, float distance) noexcept {
    set_resolution(resolution);
    set_fov(fov);
    set_distance(distance);
}

glm::vec3 Camera::position() const noexcept {
    return _rotation_to_world * glm::vec3{0.0f, 0.0f, _distance};
}

void Camera::set_resolution(glm::uvec2 res) noexcept {
    _pixel_scale = 2.0f / static_cast<float>(res.y);
}

void Camera::set_fov(float fov) noexcept {
    _z_plane = 1.0f / std::tan(0.5f * glm::radians(fov));
}

void Camera::set_distance(float d) noexcept { _distance = d; }

void Camera::rotate_x(float angle) noexcept {
    _rotation_to_world = glm::mat3{glm::rotate(glm::mat4{_rotation_to_world}, glm::radians(angle), glm::vec3{1.0f, 0.0f, 0.0f})};
}

void Camera::rotate_y(float angle) noexcept {
    _rotation_to_world = glm::mat3{glm::rotate(glm::mat4{_rotation_to_world}, glm::radians(angle), glm::vec3{0.0f, 1.0f, 0.0f})};
}

void Camera::rotate_z(float angle) noexcept {
    _rotation_to_world = glm::mat3{glm::rotate(glm::mat4{_rotation_to_world}, glm::radians(angle), glm::vec3{0.0f, 0.0f, 1.0f})};
}
