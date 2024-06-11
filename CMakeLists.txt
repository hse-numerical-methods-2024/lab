cmake_minimum_required(VERSION 3.22.1)
project(prog2 CXX)

set(CMAKE_CXX_STANDARD 23)

add_executable(prog2 prog2.cpp
        Consts.hpp)

target_link_libraries(prog2 m gsl)