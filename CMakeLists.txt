# set cmake version
cmake_minimum_required(VERSION 3.10)

if (APPLE)
    set(CMAKE_CXX_COMPILER "mpicxx")
endif()

# Find C++ compiler
enable_language(CXX)

# Find and enable OpenMPI support
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

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