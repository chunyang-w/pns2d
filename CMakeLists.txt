# set cmake version
cmake_minimum_required(VERSION 3.10)

# project name
project(NS_Solver)

# add include path
include_directories(inc)

# add the src path
add_subdirectory(src)

# add test
include(CTest)
enable_testing()
add_subdirectory(test)