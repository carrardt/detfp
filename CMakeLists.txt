cmake_minimum_required (VERSION 3.0)
project (IF64)

set(ENABLE_OPENMP "true" CACHE BOOL "Enable OpenMP")
set(ENABLE_MPI "true" CACHE BOOL "Enable MPI")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

set(IS64_SOURCES src/f64Sum.cc src/if64Sum.cc)

if(ENABLE_OPENMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)    
	    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()

if(ENABLE_MPI)
	find_package(MPI)
	if (MPI_FOUND)
		set(IS64_SOURCES ${IS64_SOURCES} src/if64AllReduceSum.cc)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
		include_directories(${MPI_INCLUDE_PATH})
	endif()
endif()

add_definitions("-mlzcnt")
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -mlzcnt")

add_library(if64 ${IS64_SOURCES})

add_executable(benchmark src/benchmark.cc)
target_link_libraries(benchmark if64)

add_executable(debugadd src/debugadd.cc)
target_link_libraries(debugadd if64)

add_executable(testlzcnt src/testlzcnt.c)

if (MPI_FOUND)
	add_executable(mpibenchmark src/mpibenchmark.cc)
	target_link_libraries(mpibenchmark if64 ${MPI_LIBRARIES})
endif()

