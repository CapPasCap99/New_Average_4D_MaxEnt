cmake_minimum_required(VERSION 3.6)
project(4dchrom_Lucas)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(/home/capucine/Documents/test/4D/AverageChromMaxEnt4D/eigen-3.4.0)

set(SOURCE_FILES main.cpp initialization.h moves.h energy_changes.h global.h RandomGenerator.h)
add_executable(Inverse_4D_new ${SOURCE_FILES})
