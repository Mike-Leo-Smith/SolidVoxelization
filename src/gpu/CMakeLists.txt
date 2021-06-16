add_library(sv-gpu STATIC common.h cuda_octree.cu cuda_octree.h)
target_link_libraries(sv-gpu PUBLIC sv-ext CUDA::cudart)
target_compile_definitions(sv-gpu PUBLIC SV_CUDA_AVAILABLE=1)
