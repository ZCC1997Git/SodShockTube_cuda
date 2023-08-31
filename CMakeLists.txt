cmake_minimum_required(VERSION 3.5)
project(SodShockWave LANGUAGES C CXX CUDA)
SET(CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_language(CUDA)
find_package(CUDA REQUIRED)
set(CMAKE_CUDA_ARCHITECTURES off)
# set(CUDA_NVCC_FLAGS "-rdc=true")

file(GLOB_RECURSE CPU_sources CONFIGURE_DEPENDS ./src/CPU/*.cpp)
file(GLOB_RECURSE CPU_header CONFIGURE_DEPENDS ./src/CPU/include/*.h)
file(GLOB_RECURSE GPU_sources CONFIGURE_DEPENDS ./src/GPU/*.cu)
file(GLOB_RECURSE GPU_header CONFIGURE_DEPENDS ./src/GPU/include/*.cuh)

add_executable(SodWave_CPU ${CPU_sources})
target_include_directories(SodWave_CPU PUBLIC ./src/CPU/include)
add_executable(SodWave_GPU ${GPU_sources})
target_include_directories(SodWave_GPU PUBLIC ./src/GPU/include)
target_compile_options(SodWave_GPU PUBLIC -rdc=true)
set_target_properties(SodWave_GPU PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON
                             RUNTIME_OUTPUT_DIRECTORY
                                 "${CMAKE_SOURCE_DIR}/build")

# find python
find_package(Python3 COMPONENTS Interpreter)