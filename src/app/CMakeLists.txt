add_executable(SolidVoxelization main.cpp)
target_link_libraries(SolidVoxelization PRIVATE sv-core-cpu)
target_compile_definitions(SolidVoxelization PRIVATE PROJECT_BASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../..")

if (SV_ENABLE_CUDA)
    add_executable(SolidVoxelizationGPU main-gpu.cpp)
    target_link_libraries(SolidVoxelizationGPU PRIVATE sv-core-gpu)
    target_compile_definitions(SolidVoxelizationGPU PRIVATE PROJECT_BASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../..")
endif ()
