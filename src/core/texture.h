//
// Created by Mike Smith on 2021/6/13.
//

#pragma once

#include <iostream>
#include <memory>
#include <span>
#include <type_traits>

#include <glad/glad.h>
#include <glm/glm.hpp>

template<typename T>
struct always_false : std::false_type {};

template<typename T>
constexpr auto always_false_v = always_false<T>::value;

template<typename T>
class Texture {

public:
    static constexpr auto pixel_format = []() -> uint32_t {
        if constexpr (std::is_same_v<T, uint8_t> || std::is_same_v<T, float>) {
            return GL_RED;
        } else if constexpr (std::is_same_v<T, glm::u8vec2> || std::is_same_v<T, glm::vec2>) {
            return GL_RG;
        } else if constexpr (std::is_same_v<T, glm::u8vec3> || std::is_same_v<T, glm::vec3>) {
            return GL_RGB;
        } else if constexpr (std::is_same_v<T, glm::u8vec4> || std::is_same_v<T, glm::vec4>) {
            return GL_RGBA;
        } else {
            static_assert(always_false_v<T>);
        }
    }();

    static constexpr auto data_format = []() -> uint32_t {
        if constexpr (std::is_same_v<T, uint8_t> || std::is_same_v<T, glm::u8vec2> || std::is_same_v<T, glm::u8vec3> || std::is_same_v<T, glm::u8vec4>) {
            return GL_UNSIGNED_BYTE;
        } else if constexpr (std::is_same_v<T, float> || std::is_same_v<T, glm::vec2> || std::is_same_v<T, glm::vec3> || std::is_same_v<T, glm::vec4>) {
            return GL_FLOAT;
        } else {
            static_assert(always_false_v<T>);
        }
    }();

private:
    glm::uvec2 _size;
    uint32_t _handle{0u};
    mutable std::unique_ptr<T[]> _pixels;

private:
    [[nodiscard]] auto _pixel_cache() const noexcept {
        if (_pixels == nullptr) { _pixels = std::make_unique<T[]>(pixel_count()); }
        return std::span{_pixels.get(), pixel_count()};
    }

public:
    explicit Texture(glm::uvec2 size) noexcept : _size{size} {
        glGenTextures(1, &_handle);
        glBindTexture(GL_TEXTURE_2D, _handle);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexImage2D(GL_TEXTURE_2D, 0, pixel_format, _size.x, _size.y, 0, pixel_format, data_format, nullptr);
    }

    explicit Texture(uint32_t width, uint32_t height) noexcept : Texture{glm::uvec2{width, height}} {}

    ~Texture() noexcept {
        if (_handle != 0u) { glDeleteTextures(1, &_handle); }
    }

    [[nodiscard]] auto size() const noexcept { return _size; }
    [[nodiscard]] auto handle() const noexcept { return _handle; }
    [[nodiscard]] auto pixel_count() const noexcept { return static_cast<size_t>(_size.x * _size.y); }

    void upload(std::span<const T> pixels) noexcept {
        if (pixels.size() < _size.x * _size.y) {
            std::cerr << "Invalid pixel count " << pixels.size() << std::endl;
            exit(-1);
        }
        glBindTexture(GL_TEXTURE_2D, _handle);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, _size.x, _size.y, pixel_format, data_format, pixels.data());
    }

    void download(std::span<T> pixels) const noexcept {
        if (pixels.size() < _size.x * _size.y) {
            std::cerr << "Invalid pixel count " << pixels.size() << std::endl;
            exit(-1);
        }
        glBindTexture(GL_TEXTURE_2D, _handle);
        glGetTexImage(GL_TEXTURE_2D, 0, pixel_format, data_format, pixels.data());
    }

    template<typename Update, std::enable_if_t<std::is_invocable_v<Update, std::span<T>>, int> = 0>
    void with_pixels_uploading(Update &&update) {
        update(_pixel_cache());
        upload(_pixel_cache());
    }

    template<typename Update, std::enable_if_t<std::is_invocable_v<Update, std::span<const T>>, int> = 0>
    void with_pixels_downloaded(Update &&update) const {
        download(_pixel_cache());
        update(_pixel_cache());
    }
};
