cmake_minimum_required(VERSION 3.16)
project(CMPMC VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Threads REQUIRED)

add_library(cmpmc INTERFACE)
target_include_directories(cmpmc INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cmpmc INTERFACE Threads::Threads)

add_executable(phonon_benchmark examples/phonon_benchmark.cpp)
target_link_libraries(phonon_benchmark PRIVATE cmpmc)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(phonon_benchmark PRIVATE -O3 -march=native -Wall -Wextra)
endif()
