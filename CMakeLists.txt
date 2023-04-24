cmake_minimum_required(VERSION 3.5)

project(GTMesh-example LANGUAGES CXX)
project(GTMesh-example-omp_parallel_graph LANGUAGES CXX)
project(GTMesh-example-omp_parallel_struct LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_SOURCE_DIR}/libs/GTMesh/src)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(GTMesh-example src/problems/heat_conduction.cpp)
target_sources(GTMesh-example PRIVATE src/solvers/RKMSolver.hpp)


find_package(OpenMP REQUIRED)
add_executable(GTMesh-example-omp_parallel_graph src/problems/heat_conduction_omp_graph.cpp)
target_sources(GTMesh-example-omp_parallel_graph PRIVATE src/solvers/RKMSolver_omp.hpp)
target_link_libraries(GTMesh-example-omp_parallel_graph PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(GTMesh-example-omp_parallel_graph PRIVATE "${OpenMP_CXX_FLAGS}")


add_executable(GTMesh-example-omp_parallel_struct src/problems/heat_conduction_omp_structure.cpp)
target_sources(GTMesh-example-omp_parallel_struct PRIVATE src/solvers/RKMSolver_omp.hpp)
target_link_libraries(GTMesh-example-omp_parallel_struct PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(GTMesh-example-omp_parallel_struct PRIVATE "${OpenMP_CXX_FLAGS}")
