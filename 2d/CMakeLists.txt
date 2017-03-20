# Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

# Open terminal
# cd path_to_the_folder/2d/bin/
# cmake -DCMAKE_BUILD_TYPE=Debug .. or cmake -DCMAKE_BUILD_TYPE=Release ..
# make
# ./blmv

cmake_minimum_required(VERSION 3.3)

project(blmv)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -std=c++11 ")

set(SRC 
        extra/BarycentricCoordinatesR2.hpp
        extra/SegmentCoordinatesR2.hpp
        extra/TriangleCoordinatesR2.hpp
        extra/VertexExpressionsR2.hpp
        extra/VertexR2.hpp
        extra/Halfedge.hpp
        extra/Face.hpp
        extra/MeshR2.hpp
        extra/TriangulatorR2.hpp
        extra/libs/triangle/triangle.hpp
        extra/libs/triangle/triangle.cpp
        coords/BlendedR2.hpp
        main.cpp
)

add_executable(blmv ${SRC})
