# Standalone project that consumes an installed SeQuant via find_package;
# not part of the SeQuant build, see the consumer step in .github/workflows/cmake.yml
cmake_minimum_required(VERSION 3.28)
project(sequant-consumer CXX)

find_package(SeQuant CONFIG REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE SeQuant::SeQuant)
