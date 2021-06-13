// Dear ImGui: standalone example application for GLFW + OpenGL 3, using programmable pipeline
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)
// If you are new to Dear ImGui, read documentation from the docs/ folder + read the top of imgui.cpp.
// Read online: https://github.com/ocornut/imgui/tree/master/docs

#include <iostream>
#include <random>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/color_space.hpp>
#include <stb/stb_image_write.h>

#include <core/camera.h>
#include <core/framerate.h>
#include <core/mesh.h>
#include <core/stream.h>
#include <core/texture.h>
#include <core/volume.h>
#include <core/window.h>

template<typename F>
void with_imgui_window(const char *name, bool *is_open, ImGuiWindowFlags flags, F &&f) {
    if (ImGui::Begin(name, is_open, flags)) { f(); }
    ImGui::End();
}

int main(int, char **) {

    auto mesh = Mesh::load("data/bunny.obj");
    Stream stream;

    static constexpr auto resolution = 512;

    std::vector<RayHit> rays;
    rays.resize(resolution * resolution);

    static constexpr auto rand = [] {
        static thread_local std::mt19937 random{std::random_device{}()};
        return std::uniform_real_distribution<float>{0.0f, 1.0f}(random);
    };

    Window window{"Solid Voxelization"};
    Texture<glm::u8vec4> display_buffer{resolution, resolution};
    std::vector<glm::vec3> accum_buffer(resolution * resolution);

    glm::vec3 light_position{-3.0f, 2.0f, 5.0f};
    glm::vec3 light_emission{20.0f};
    glm::vec3 albedo{1.0f, 1.0f, 1.0f};

    float camera_fov = 35.0f;
    float camera_distance = 3.0f;
    Camera camera{glm::uvec2{resolution}, camera_fov, camera_distance};

    Framerate framerate;
    window.run([&] {
        auto ray_origin = camera.position();
        stream.dispatch_2d({resolution, resolution}, [&](glm::uvec2 xy) noexcept {
            auto x = xy.x;
            auto y = xy.y;
            auto index = y * resolution + x;
            auto &&r = rays[index];
            auto dx = rand();
            auto dy = rand();
            r.ray.o = ray_origin;
            r.ray.t_min = 0.0f;
            r.ray.t_max = std::numeric_limits<float>::infinity();
            r.ray.d = camera.direction(glm::vec2{xy} + glm::vec2{dx, dy});
            r.hit.geom_id = Hit::invalid;
        });

        stream.dispatch_1d(resolution, 16u, [&mesh, &rays](uint32_t tid) noexcept {
            mesh->trace_closest(std::span{rays}.subspan(tid * resolution, resolution));
        });

        stream.dispatch_1d(resolution * resolution, [&, vb = mesh->vertices(), ib = mesh->indices()](uint32_t tid) {
            auto hit = rays[tid].hit;
            auto radiance = [&] {
                if (hit.geom_id == Hit::invalid) { return glm::vec3{}; }
                auto ng = glm::normalize(hit.ng);
                auto tri = ib[hit.prim_id];
                auto p0 = vb[tri.x];
                auto p1 = vb[tri.y];
                auto p2 = vb[tri.z];
                auto p = (1.0f - hit.uv.x - hit.uv.y) * p0 + hit.uv.x * p1 + hit.uv.y * p2;
                auto L = light_position - p;
                auto inv_dd = 1.0f / glm::dot(L, L);
                L = glm::normalize(L);
                auto ray = rays[tid].ray;
                ray.o = p + 1e-4f * ng;
                ray.d = L;
                ray.t_min = 0.0f;
                ray.t_max = std::numeric_limits<float>::max();
                mesh->trace_any(&ray);
                auto radiance = ray.t_max < 0.0f
                                    ? glm::vec3{}
                                    : glm::max(glm::dot(L, ng), 0.0f) * light_emission * inv_dd * albedo;
                return radiance + 0.003f * albedo;
            }();
            auto &&accum = accum_buffer[tid];
            auto t = 1.0f / static_cast<float>(framerate.count() + 1u);
            accum = glm::mix(accum, radiance, t);
        });

        display_buffer.with_pixels_uploading([&](auto pixels) noexcept {
            stream.dispatch_1d(resolution * resolution, [&](auto tid) noexcept {
                pixels[tid] = {glm::clamp(glm::convertLinearToSRGB(accum_buffer[tid]) * 255.0f, 0.0f, 255.0f), 255};
            });
        });

        auto dt = framerate.tick();
        auto fps = framerate.fps();
        auto spp = framerate.count();

        with_imgui_window("Console", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize, [&] {
            ImGui::Text("Samples: %llu", spp);
            ImGui::SameLine();
            ImGui::Spacing();
            ImGui::SameLine();
            ImGui::Text("Time: %.2lf ms", dt);
            ImGui::SameLine();
            ImGui::Spacing();
            ImGui::SameLine();
            ImGui::Text("FPS: %.2lf", fps);
            ImGui::Text("Camera Distance: %.2f", camera_distance);
            ImGui::SameLine();
            ImGui::Spacing();
            ImGui::SameLine();
            ImGui::Text("FoV: %.2f", camera_fov);
            if (ImGui::SliderFloat3("Light Position", &light_position.x, -10.0f, 10.0f, "%.1f")) { framerate.clear(); }
            if (ImGui::SliderFloat3("Light Emission", &light_emission.x, 0.0f, 100.0f, "%.1f")) { framerate.clear(); }
            if (ImGui::ColorEdit3("Mesh Albedo", &albedo.x, ImGuiColorEditFlags_Float)) { framerate.clear(); }
        });

        if (glfwGetKey(window.handle(), GLFW_KEY_Q) == GLFW_PRESS) {
            display_buffer.with_pixels_downloaded([](auto pixels) noexcept {
                stbi_write_png("test.png", resolution, resolution, 4, pixels.data(), 0);
            });
            window.notify_close();
        }

        if (glfwGetKey(window.handle(), GLFW_KEY_LEFT) == GLFW_PRESS) {
            camera.rotate_y(-3.0f);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_RIGHT) == GLFW_PRESS) {
            camera.rotate_y(3.0f);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_UP) == GLFW_PRESS) {
            camera.rotate_x(-3.0f);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_DOWN) == GLFW_PRESS) {
            camera.rotate_x(3.0f);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_COMMA) == GLFW_PRESS) {
            camera.rotate_z(-3.0f);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_PERIOD) == GLFW_PRESS) {
            camera.rotate_z(3.0f);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_MINUS) == GLFW_PRESS) {
            camera_distance = std::max(camera_distance / 1.02f, 1.0f);
            camera.set_distance(camera_distance);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_EQUAL) == GLFW_PRESS) {
            camera_distance = camera_distance * 1.02f;
            camera.set_distance(camera_distance);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_Z) == GLFW_PRESS) {
            camera_fov = std::max(camera_fov / 1.05f, 1.0f);
            camera.set_fov(camera_fov);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_X) == GLFW_PRESS) {
            camera_fov = std::min(camera_fov * 1.05f, 150.0f);
            camera.set_fov(camera_fov);
            framerate.clear();
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_SPACE) == GLFW_PRESS) {
            auto volume = Volume::from(*mesh, 256u, camera.rotation_to_world());// TODO...
        }

        with_imgui_window("Mesh", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize, [&] {
            ImGui::Image(display_buffer.handle(), {resolution, resolution});
        });
    });
}
