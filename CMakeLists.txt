cmake_minimum_required (VERSION 2.8)
project (IF64)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_definitions("-std=c++11 -mlzcnt")
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -mlzcnt")

add_executable(benchmark main.cc if64Sum.cc if64AllReduceSum.cc)
add_executable(testlzcnt testlzcnt.c)

