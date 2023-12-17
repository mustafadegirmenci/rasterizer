cmake_minimum_required(VERSION 3.25)
project(ray_tracer)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "include/*.h")

add_executable(rasterizer ${SOURCES} ${HEADERS} include/DummyLine.h src/DummyLine.cpp include/DummyTriangle.h src/DummyTriangle.cpp)
