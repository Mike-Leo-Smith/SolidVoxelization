add_library(imgui STATIC
            imgui/imgui.cpp imgui/imgui.h
            imgui/imgui_demo.cpp
            imgui/imgui_draw.cpp
            imgui/imgui_tables.cpp
            imgui/imgui_widgets.cpp
            imgui/backends/imgui_impl_glfw.cpp imgui/backends/imgui_impl_glfw.h
            imgui/backends/imgui_impl_opengl3.cpp imgui/backends/imgui_impl_opengl3.h)
target_include_directories(imgui PUBLIC imgui "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(imgui PUBLIC glad glfw)
target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLAD ImTextureID=uint32_t)
