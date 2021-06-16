include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

include(CheckLanguage)
check_language(CUDA)

if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "Compile with CUDA")
    set(SV_ENABLE_CUDA ON)
endif()

if (SV_ENABLE_CUDA)
    add_subdirectory(gpu)
endif ()

add_subdirectory(core)
add_subdirectory(app)
