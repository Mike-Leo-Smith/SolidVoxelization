add_library(glad STATIC glad/glad.c glad/glad.h glad/khrplatform.h)
target_link_libraries(glad PUBLIC OpenGL::GL)
target_include_directories(glad PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
