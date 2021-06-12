//
// Created by Mike Smith on 2021/6/12.
//

#pragma once

#include <cstdint>

struct alignas(8) Segment {
    uint16_t x;
    uint16_t y;
    uint16_t z_min;
    uint16_t z_max;
};


