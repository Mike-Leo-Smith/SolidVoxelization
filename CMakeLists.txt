cmake_minimum_required(VERSION 3.20)
project(SolidVoxelization)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_subdirectory(ext)
add_subdirectory(src)
