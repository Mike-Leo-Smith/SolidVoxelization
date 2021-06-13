//
// Created by Mike Smith on 2021/6/12.
//

#pragma once

#include <thread>
#include <vector>
#include <queue>
#include <memory>
#include <functional>

#include <glm/glm.hpp>

class ThreadPool;

class Stream {

public:
    using Task1D = std::function<void(uint32_t)>;
    using Task2D = std::function<void(glm::uvec2)>;

private:
    std::unique_ptr<ThreadPool> _thread_pool;

public:
    explicit Stream(size_t num_workers = 0u) noexcept;
    ~Stream() noexcept;
    void dispatch_1d(uint32_t dispatch_size, uint32_t block_size, Task1D task) noexcept;
    void dispatch_2d(glm::uvec2 dispatch_size, glm::uvec2 block_size, Task2D task) noexcept;
    void dispatch_1d(uint32_t dispatch_size, Task1D task) noexcept { dispatch_1d(dispatch_size, 256u, std::move(task)); }
    void dispatch_2d(glm::uvec2 dispatch_size, Task2D task) noexcept { dispatch_2d(dispatch_size, {16u, 16u}, std::move(task)); }
};
