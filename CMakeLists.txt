cmake_minimum_required(VERSION 3.8)
project(SDF)

# Set C++ version
set(CMAKE_CXX_STANDARD 17)
# Set verbose makefile
set(CMAKE_VERBOSE_MAKEFILE ON)

set(SOURCE_FILES main.cpp common.hpp sdf.hpp lodepng.cpp lodepng.hpp plotting.hpp grid.hpp)
add_executable(SDF ${SOURCE_FILES})