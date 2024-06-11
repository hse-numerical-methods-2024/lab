cmake_minimum_required(VERSION 3.22.1)
project(prog1 CXX)

set(CMAKE_CXX_STANDARD 23)

add_executable(prog1 prog1.cpp
        Consts.hpp
        Exp.hpp)

target_link_libraries(prog1 m)