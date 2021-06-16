#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>

#include <cuda_runtime.h>
#include <core/ray.h>

#ifdef __CUDACC__
#define SV_GPU_CODE
#endif

#ifdef SV_GPU_CODE
#define SV_XPU __host__ __device__
#else
#define SV_XPU
#endif

namespace detail {
    inline void check_cuda(cudaError_t error, const char *file, int line) noexcept {
        if (error != cudaSuccess) {
            fprintf(stderr, "CUDA Error: %s [%s:$d]", cudaGetErrorString(error), file, line);
            exit(-1);
        }
    }
}

#define CHECK_CUDA(...) detail::check_cuda(__VA_ARGS__, __FILE__, __LINE__)

template<typename T>
class CUDABuffer {

private:
    T *_data{nullptr};
    size_t _size;

    explicit CUDABuffer(size_t size) noexcept
        : _size{size} { CHECK_CUDA(cudaMalloc(&_data, sizeof(T) * size)); }

public:
    CUDABuffer(CUDABuffer &&) noexcept = delete;
    CUDABuffer(const CUDABuffer &) noexcept = delete;
    CUDABuffer &operator=(CUDABuffer &&) noexcept = delete;
    CUDABuffer &operator=(const CUDABuffer &) noexcept = delete;

    ~CUDABuffer() noexcept { cudaFree(_data); }
    [[nodiscard]] static auto create(size_t size) noexcept { return std::unique_ptr<CUDABuffer>{new CUDABuffer{size}}; }
    [[nodiscard]] auto size() const noexcept { return _size; }
    [[nodiscard]] auto size_bytes() const noexcept { return _size * sizeof(T); }
    [[nodiscard]] auto data() noexcept { return _data; }
    [[nodiscard]] auto data() const noexcept { return const_cast<const T *>(_data); }
    void upload(const T *data) noexcept { CHECK_CUDA(cudaMemcpyAsync(_data, data, size_bytes(), cudaMemcpyHostToDevice)); }
    void download(T *data) noexcept { CHECK_CUDA(cudaMemcpyAsync(data, _data, size_bytes(), cudaMemcpyDeviceToHost)); }
};
