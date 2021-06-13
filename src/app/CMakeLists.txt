add_executable(SolidVoxelization main.cpp)
target_link_libraries(SolidVoxelization PRIVATE sv-core)
target_compile_definitions(SolidVoxelization PRIVATE PROJECT_BASE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../..")
