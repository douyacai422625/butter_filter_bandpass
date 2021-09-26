cmake_minimum_required(VERSION 3.20)
project(ButterFilter)

set(CMAKE_CXX_STANDARD 14)

add_executable(ButterFilter Signal_Filter.cpp Signal_Filter.h main.cpp)

target_link_libraries(ButterFilter)