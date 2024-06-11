cmake_minimum_required(VERSION 3.22.1)
project(prog3 CXX)

set(CMAKE_CXX_STANDARD 23)

add_executable(prog3 prog3.cpp
        Consts.hpp)

target_link_libraries(prog3 m gsl)