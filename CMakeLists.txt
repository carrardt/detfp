cmake_minimum_required (VERSION 2.8)
project (IF64)

find_package(OpenMP)
find_package(MPI)

set(IS64_SOURCES f64Sum.cc if64Sum.cc)

if (MPI_FOUND)
	set(IS64_SOURCES ${IS64_SOURCES} if64AllReduceSum.cc)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp-simd")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp-simd")
endif()

add_definitions("-mlzcnt")
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -mlzcnt")

add_library(if64 ${IS64_SOURCES})

add_executable(benchmark benchmark.cc)
target_link_libraries(benchmark if64)

add_executable(testlzcnt testlzcnt.c)

