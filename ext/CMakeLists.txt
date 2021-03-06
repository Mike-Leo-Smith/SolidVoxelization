add_library(sv-ext INTERFACE)

if (MSVC)
    add_compile_options("/wd5105")
endif ()

set(ASSIMP_NO_EXPORT ON CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ASSIMP_INSTALL OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ALL_IMPORTERS_BY_DEFAULT OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_OBJ_IMPORTER ON CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_PLY_IMPORTER ON CACHE BOOL "" FORCE)
add_subdirectory(assimp)
target_link_libraries(sv-ext INTERFACE assimp::assimp)

add_subdirectory(glm)
target_link_libraries(sv-ext INTERFACE glm::glm)

find_package(embree 3.0 CONFIG)
if (NOT embree_FOUND)
    set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "" FORCE)
    set(EMBREE_GEOMETRY_QUAD OFF CACHE BOOL "" FORCE)
    set(EMBREE_GEOMETRY_CURVE OFF CACHE BOOL "" FORCE)
    set(EMBREE_GEOMETRY_SUBDIVISION OFF CACHE BOOL "" FORCE)
    set(EMBREE_GEOMETRY_USER OFF CACHE BOOL "" FORCE)
    set(EMBREE_GEOMETRY_INSTANCE OFF CACHE BOOL "" FORCE)
    set(EMBREE_GEOMETRY_GRID OFF CACHE BOOL "" FORCE)
    set(EMBREE_GEOMETRY_POINT OFF CACHE BOOL "" FORCE)
    set(EMBREE_TASKING_SYSTEM "INTERNAL" CACHE STRING "" FORCE)
    set(EMBREE_TUTORIALS OFF CACHE BOOL "" FORCE)
    if (APPLE)
        set(EMBREE_STATIC_LIB OFF CACHE BOOL "" FORCE)
    else ()
        set(EMBREE_STATIC_LIB ON CACHE BOOL "" FORCE)
    endif ()
    add_subdirectory(embree)
endif ()
target_link_libraries(sv-ext INTERFACE embree)

find_package(OpenGL REQUIRED)
target_link_libraries(sv-ext INTERFACE OpenGL::GL)

find_package(GLFW3 CONFIG)
if (NOT GLFW3_FOUND)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF)
    add_subdirectory(glfw)
endif ()
target_link_libraries(sv-ext INTERFACE glfw)

add_subdirectory(glad)
target_link_libraries(sv-ext INTERFACE glad)

add_subdirectory(imgui)
target_link_libraries(sv-ext INTERFACE imgui)

add_subdirectory(stb)
target_link_libraries(sv-ext INTERFACE stb)
