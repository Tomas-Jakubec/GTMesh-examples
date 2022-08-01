cmake_minimum_required(VERSION 3.5)

project(GTMesh-example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(GTMesh-example main.cpp)
include_directories(${CMAKE_SOURCE_DIR}/GTMesh/src)
