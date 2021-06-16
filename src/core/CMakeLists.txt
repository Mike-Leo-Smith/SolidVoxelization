set(SV_SOURCES
    mesh.cpp mesh.h
    ray.h
    bbox.h
    segment.h
    volume.cpp volume.h
    stream.cpp stream.h
    texture.h
    window.cpp window.h
    camera.cpp camera.h
    framerate.cpp framerate.h)

add_library(sv-core-cpu STATIC ${SV_SOURCES})
target_link_libraries(sv-core-cpu PUBLIC sv-ext)

if (SV_ENABLE_CUDA)
    add_library(sv-core-gpu STATIC ${SV_SOURCES})
    target_link_libraries(sv-core-gpu PUBLIC sv-gpu)
endif ()