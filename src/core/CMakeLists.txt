add_library(sv-core STATIC
            mesh.cpp mesh.h
            ray.h
            bbox.h
            segment.h
            volume.cpp volume.h
            stream.cpp stream.h
            texture.h window.cpp window.h)
target_link_libraries(sv-core PUBLIC sv-ext)
