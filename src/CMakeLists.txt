include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(SolidVoxelization main.cpp mesh.cpp mesh.h ray.h bbox.h segment.h volume.cpp volume.h)
target_link_libraries(SolidVoxelization PRIVATE sv-ext)
