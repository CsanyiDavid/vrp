cmake_minimum_required(VERSION 3.10)
project(vrp)

add_executable(vrp main.cc)
target_link_libraries(vrp emon)