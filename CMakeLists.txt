cmake_minimum_required(VERSION 2.8)

project(masbcpp)

# Figure if we can enable OpenMP
set(LINK_LIBS)
set(COMPILE_OPTIONS "-funroll-loops -ffast-math")

if(EXISTS "/usr/local/bin/clang-omp++")
	set(CMAKE_CXX_COMPILER "/usr/local/bin/clang-omp++" )
	add_definitions(-DWITH_OPENMP)
	add_definitions(-DCLANG_OMP)

	set(LINK_LIBS ${LINK_LIBS} -liomp5)
	set(COMPILE_OPTIONS "${COMPILE_OPTIONS} -fopenmp")
else()
	find_package(OpenMP)
	if(OPENMP_FOUND)
		add_definitions(-DWITH_OPENMP)
		set(COMPILE_OPTIONS "${COMPILE_OPTIONS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()

# global
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COMPILE_OPTIONS}")
set(CMAKE_BUILD_TYPE Release)
set_property(GLOBAL PROPERTY LINKER_LANGUAGE CXX)

# include directories
include_directories(${CMAKE_SOURCE_DIR}/thirdparty ${CMAKE_SOURCE_DIR}/thirdparty/vrui)

# bundle all thirdparty stuff as a static library
FILE(GLOB_RECURSE THIRDPARTY thirdparty/*.cpp)
add_library(thirdparty STATIC ${THIRDPARTY})
set(LINK_LIBS ${LINK_LIBS} thirdparty)

# set targets
add_executable(compute_ma src/compute_ma.cpp)
add_executable(compute_normals src/compute_normals.cpp)

# link targets
target_link_libraries(compute_ma ${LINK_LIBS})
target_link_libraries(compute_normals ${LINK_LIBS})

install(TARGETS compute_ma DESTINATION bin)