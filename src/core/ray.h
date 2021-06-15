//
// Created by Mike Smith on 2021/6/12.
//

#pragma once

#include <limits>
#include <glm/glm.hpp>

struct alignas(16) EmbreeRay {
    glm::vec3 o;
    float t_min;
    glm::vec3 d;
    float time;
    float t_max;
    unsigned int mask;
    unsigned int id;
    unsigned int flags;
};

struct alignas(16) EmbreeHit {
    
    static constexpr auto invalid = std::numeric_limits<uint32_t>::max();

    glm::vec3 ng;
    glm::vec2 uv;
    unsigned int prim_id;
    unsigned int geom_id;
    unsigned int inst_id;
};

struct EmbreeRayHit {
    EmbreeRay ray;
    EmbreeHit hit;
};

struct alignas(16) Ray {
    glm::vec3 o;
    float t_min;
    glm::vec3 d;
    float t_max;
};

struct alignas(16) Hit {
    glm::vec3 p;
    float t;
    glm::vec3 ng;
    int valid;
};
