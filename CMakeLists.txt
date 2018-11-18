cmake_minimum_required(VERSION 2.8)

project(masbcpp)

# Figure out if we can enable OpenMP
set(LINK_LIBS)
set(COMPILE_OPTIONS "-funroll-loops -ffast-math")

# disable openmp with xcode, doesn't work
if(NOT CMAKE_GENERATOR MATCHES "Xcode")
if(APPLE)
	if(EXISTS "/usr/local/bin/clang-omp++")
		message(STATUS "Building with OMP support")

		set(CMAKE_CXX_COMPILER "/usr/local/bin/clang-omp++" )
		add_definitions(-DWITH_OPENMP)
		#add_definitions(-DCLANG_OMP)

		set(LINK_LIBS ${LINK_LIBS} -liomp5)
		set(COMPILE_OPTIONS "${COMPILE_OPTIONS} -fopenmp")
		include_directories(/usr/local/include/libiomp)
		link_directories(/usr/local/lib)
	endif()
else()
	find_package(OpenMP)
	if(OPENMP_FOUND)
	
		add_definitions(-DWITH_OPENMP)
		set(COMPILE_OPTIONS "${COMPILE_OPTIONS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()
endif()

# enable verbose printing by default
add_definitions(-DVERBOSEPRINT)

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

# build a library from the masbpcpp processing functions
add_library(masbcpp STATIC src/compute_ma_processing.cpp src/compute_normals_processing.cpp src/simplify_processing.cpp)

# set excutables
add_executable(compute_ma src/compute_ma.cpp)
add_executable(compute_normals src/compute_normals.cpp)
add_executable(simplify src/simplify.cpp)
set_target_properties(
  compute_ma compute_normals simplify thirdparty masbcpp
  PROPERTIES CXX_STANDARD 11
)

# link targets
target_link_libraries(masbcpp ${LINK_LIBS})

target_link_libraries(compute_ma masbcpp)
target_link_libraries(compute_normals masbcpp)
target_link_libraries(simplify masbcpp)

install(TARGETS compute_ma compute_normals simplify DESTINATION bin)
