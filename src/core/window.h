//
// Created by Mike Smith on 2021/6/13.
//

#pragma once

#include <type_traits>

#include <glm/glm.hpp>
#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <imgui/imgui.h>

class Window {

public:
    using Handle = GLFWwindow *;

private:
    Handle _handle{nullptr};

private:
    static void _begin_frame(Handle) noexcept;
    static void _end_frame(Handle) noexcept;

public:
    Window(const char *title, uint32_t width, uint32_t height) noexcept;
    Window(const char *title, glm::uvec2 size = {1280, 720}) noexcept : Window{title, size.x, size.y} {}
    ~Window() noexcept;

    [[nodiscard]] auto handle() const noexcept { return _handle; }
    [[nodiscard]] bool should_close() const noexcept;
    
    void notify_close() noexcept;
    
    template<typename Frame, std::enable_if_t<std::is_invocable_v<Frame>, int> = 0>
    void with_frame(Frame &&frame) {
        _begin_frame(_handle);
        frame();
        _end_frame(_handle);
    }

    template<typename Frame, std::enable_if_t<std::is_invocable_v<Frame>, int> = 0>
    void run(Frame &&frame) {
        while (!should_close()) { with_frame(std::forward<Frame>(frame)); }
    }
};
