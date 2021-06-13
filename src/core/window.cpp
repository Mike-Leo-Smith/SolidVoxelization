//
// Created by Mike Smith on 2021/6/13.
//

#include <iostream>

#include <glad/glad.h>
#include <imgui/backends/imgui_impl_glfw.h>
#include <imgui/backends/imgui_impl_opengl3.h>

#include <core/window.h>

static void glfw_error_callback(int error, const char *description) {
    std::cerr << "Glfw Error " << error << ": " << description << std::endl;
}

Window::Window(const char *title, uint32_t width, uint32_t height) noexcept {

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    if (_handle = glfwCreateWindow(static_cast<int>(width), static_cast<int>(height), title, nullptr, nullptr); _handle == nullptr) {
        std::cerr << "Failed to create window" << std::endl;
        exit(-1);
    }
    glfwMakeContextCurrent(_handle);
    glfwSwapInterval(1);

    if (gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress)) == 0) {
        std::cerr << "Failed to initialize OpenGL loader" << std::endl;
        exit(-1);
    }

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsLight();
    ImGui_ImplGlfw_InitForOpenGL(_handle, true);
    ImGui_ImplOpenGL3_Init("#version 150");
}

void Window::_begin_frame(Window::Handle) noexcept {
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void Window::_end_frame(Window::Handle handle) noexcept {
    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(handle, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(handle);
}

Window::~Window() noexcept {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(_handle);
    glfwTerminate();
}

bool Window::should_close() const noexcept { return glfwWindowShouldClose(_handle); }
void Window::notify_close() noexcept { glfwSetWindowShouldClose(_handle, true); }
