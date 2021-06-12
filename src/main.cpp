#include <opencv2/opencv.hpp>

#include <mesh.h>
#include <stream.h>
#include <volume.h>

#include <condition_variable>

int main() {

    auto mesh = Mesh::load("test.obj");
    Stream stream;

    static constexpr auto resolution = 512;

    cv::Mat image{resolution, resolution, CV_32FC3, cv::Scalar::all(0)};
    std::vector<RayHit> rays;
    rays.resize(resolution * resolution);
    glm::mat4 R{1.0f};
    for (;;) {
        auto omega = glm::mat3{R} * glm::vec3{0.0f, 0.0f, -1.0f};
        stream.dispatch_2d({resolution, resolution}, [omega, R, &rays](glm::uvec2 xy) noexcept {
            auto x = xy.x;
            auto y = xy.y;
            auto index = y * resolution + x;
            auto &&r = rays[index];
            constexpr auto scale = 2.0f / resolution;
            auto dx = (x + 0.5f) * scale - 1.0f;
            auto dy = 1.0f - (y + 0.5f) * scale;
            auto dz = 2.0f;
            r.ray.o = glm::vec3{R * glm::vec4{dx, dy, dz, 1.0f}};
            r.ray.t_min = 0.0f;
            r.ray.t_max = std::numeric_limits<float>::infinity();
            r.ray.d = omega;
            r.hit.geom_id = Hit::invalid;
        });

        auto t0 = std::chrono::high_resolution_clock::now();
        stream.dispatch_1d(resolution * resolution, [&mesh, &rays](uint32_t tid) noexcept {
            mesh->trace_closest(&rays[tid]);
        });
        auto t1 = std::chrono::high_resolution_clock::now();

        using namespace std::chrono_literals;
        std::cout << "Rendering Time: " << (t1 - t0) / 1ns * 1e-6 << " ms" << std::endl;

        stream.dispatch_1d(resolution * resolution, [&image, &rays](uint32_t tid) {
            auto &&pixel = reinterpret_cast<glm::vec3 *>(image.data)[tid];
            auto hit = rays[tid].hit;
            if (hit.geom_id != Hit::invalid) {
                auto ng = glm::normalize(hit.ng);
                ng = glm::vec3{ng.z, ng.y, ng.x};
                pixel = ng * 0.5f + 0.5f;
            } else {
                pixel = {};
            }
        });
        cv::imshow("Display", image);
        auto key = cv::waitKey(1);
        if (key == 'q') { break; }
        if (key == 'a') {
            R = glm::rotate(R, glm::radians(5.0f), glm::vec3{0.0f, 1.0f, 0.0f});
        } else if (key == 'd') {
            R = glm::rotate(R, glm::radians(-5.0f), glm::vec3{0.0f, 1.0f, 0.0f});
        } else if (key == 'w') {
            R = glm::rotate(R, glm::radians(5.0f), glm::vec3{1.0f, 0.0f, 0.0f});
        } else if (key == 's') {
            R = glm::rotate(R, glm::radians(-5.0f), glm::vec3{1.0f, 0.0f, 0.0f});
        } else if (key == ' ') {
            auto volume = Volume::from(*mesh, 256u, R);
        }
    }
    cv::imwrite("normal.exr", image);

    return 0;
}
