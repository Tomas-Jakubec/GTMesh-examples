cmake_minimum_required(VERSION 3.5)

project(GTMesh-example LANGUAGES CXX)
project(GTMesh-example-omp_parallel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_SOURCE_DIR}/libs/GTMesh/src)

add_executable(GTMesh-example main.cpp)
target_sources(GTMesh-example PRIVATE RKMSolver.hpp)


find_package(OpenMP REQUIRED)
add_executable(GTMesh-example-omp_parallel heat_conduction_omp_graph.cpp)
target_sources(GTMesh-example-omp_parallel PRIVATE RKMSolver_omp.hpp)
target_link_libraries(GTMesh-example-omp_parallel PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(GTMesh-example-omp_parallel PRIVATE "${OpenMP_CXX_FLAGS}")
