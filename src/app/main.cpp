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
    Stream stream{16};

    static constexpr auto resolution = 512;

    std::vector<RayHit> rays;
    rays.resize(resolution * resolution);
    glm::mat4 R{1.0f};

    static constexpr auto rand = [] {
        static thread_local std::mt19937 random{std::random_device{}()};
        return std::uniform_real_distribution<float>{0.0f, 1.0f}(random);
    };

    Window window{"Solid Voxelization"};
    Texture<glm::u8vec4> display_buffer{resolution, resolution};
    std::vector<glm::vec3> accum_buffer(resolution * resolution);

    auto spp = 0u;
    window.run([&] {
        spp++;

        auto t0 = std::chrono::high_resolution_clock::now();
        auto omega = glm::mat3{R} * glm::vec3{0.0f, 0.0f, -1.0f};
        stream.dispatch_2d({resolution, resolution}, [omega, R, &rays](glm::uvec2 xy) noexcept {
            auto x = xy.x;
            auto y = xy.y;
            auto index = y * resolution + x;
            auto &&r = rays[index];
            constexpr auto scale = 2.0f / resolution;
            auto dx = (x + rand()) * scale - 1.0f;
            auto dy = 1.0f - (y + rand()) * scale;
            auto dz = 2.0f;
            r.ray.o = glm::vec3{R * glm::vec4{dx, dy, dz, 1.0f}};
            r.ray.t_min = 0.0f;
            r.ray.t_max = std::numeric_limits<float>::infinity();
            r.ray.d = omega;
            r.hit.geom_id = Hit::invalid;
        });

        stream.dispatch_1d(resolution, 16u, [&mesh, &rays](uint32_t tid) noexcept {
            mesh->trace_closest(std::span{rays}.subspan(tid * resolution, resolution));
        });

        stream.dispatch_1d(resolution * resolution, [&, vb = mesh->vertices(), ib = mesh->indices()](uint32_t tid) {
            auto hit = rays[tid].hit;
            glm::vec3 light_position{-3.0f, 2.0f, 5.0f};
            glm::vec3 light_emission{20.0f};
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
                                    : glm::max(glm::dot(L, ng), 0.0f) * light_emission * inv_dd;
                return radiance + 0.01f;
            }();
            auto &&accum = accum_buffer[tid];
            auto t = 1.0f / static_cast<float>(spp);
            accum = glm::mix(accum, radiance, t);
        });

        display_buffer.with_pixels_uploading([&](auto pixels) noexcept {
            stream.dispatch_1d(resolution * resolution, [&](auto tid) noexcept {
                pixels[tid] = {glm::clamp(glm::convertLinearToSRGB(accum_buffer[tid]) * 255.0f, 0.0f, 255.0f), 255};
            });
        });
        auto t1 = std::chrono::high_resolution_clock::now();

        using namespace std::chrono_literals;
        std::cout << "Rendering Time: " << (t1 - t0) / 1ns * 1e-6 << " ms" << std::endl;

        if (glfwGetKey(window.handle(), GLFW_KEY_Q) == GLFW_PRESS) {
            display_buffer.with_pixels_downloaded([](auto pixels) noexcept {
                stbi_write_png("test.png", resolution, resolution, 4, pixels.data(), 0);
            });
            window.notify_close();
        }

        if (glfwGetKey(window.handle(), GLFW_KEY_LEFT) == GLFW_PRESS) {
            R = glm::rotate(R, glm::radians(-5.0f), glm::vec3{0.0f, 1.0f, 0.0f});
            spp = 0u;
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_RIGHT) == GLFW_PRESS) {
            R = glm::rotate(R, glm::radians(5.0f), glm::vec3{0.0f, 1.0f, 0.0f});
            spp = 0u;
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_UP) == GLFW_PRESS) {
            R = glm::rotate(R, glm::radians(-5.0f), glm::vec3{1.0f, 0.0f, 0.0f});
            spp = 0u;
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_DOWN) == GLFW_PRESS) {
            R = glm::rotate(R, glm::radians(5.0f), glm::vec3{1.0f, 0.0f, 0.0f});
            spp = 0u;
        }
        if (glfwGetKey(window.handle(), GLFW_KEY_SPACE) == GLFW_PRESS) {
            auto volume = Volume::from(*mesh, 256u, R);// TODO...
        }

        with_imgui_window("Mesh", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoResize, [&] {
            ImGui::Image(display_buffer.handle(), {resolution, resolution});
        });
    });
}
