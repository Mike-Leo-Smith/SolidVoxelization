include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(SolidVoxelization main.cpp)
target_link_libraries(SolidVoxelization PRIVATE sv-ext)
