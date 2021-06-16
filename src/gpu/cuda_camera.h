#pragma once

#include <gpu/common.h>
#include <core/ray.h>
#include <core/camera.h>

void cuda_camera_generate_rays(
    Camera cam, glm::mat4 world_to_object,
    Ray *rays, uint32_t w, uint32_t h, uint32_t frame) noexcept;
