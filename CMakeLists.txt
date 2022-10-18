cmake_minimum_required(VERSION 3.6)
project(4dchrom_Lucas)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(/home/janni/Code/boost_1_79_0)
include_directories(/home/janni/Code/eigen-3.4.0)
#include_directories(/home/joris/Downloads/eigen-eigen-323c052e1731)
#include_directories(/home/joris/Downloads/boost_1_76_0)
set(SOURCE_FILES main.cpp initialization.h moves.h energy_changes.h global.h RandomGenerator.h)
add_executable(Inverse_3D_new ${SOURCE_FILES})
