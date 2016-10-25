cmake_minimum_required(VERSION 2.8)

project(masbcpp)

# C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

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

# Find eigen, should set EIGEN3_INCLUDE_DIR
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Eigen3 REQUIRED)

find_package(PCL 1.8 REQUIRED COMPONENTS common search features)

# enable verbose printing by default
# PCL_NO_PRECOMPILE is needed to make the pcl::NormalEstimationOMP work properly
add_definitions(${PCL_DEFINITIONS} -DVERBOSEPRINT -DPCL_NO_PRECOMPILE)

# global
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COMPILE_OPTIONS}")
set(CMAKE_BUILD_TYPE Release)
set_property(GLOBAL PROPERTY LINKER_LANGUAGE CXX)

# include directories
include_directories(${CMAKE_SOURCE_DIR}/thirdparty ${PCL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
link_directories(${PCL_LIBRARY_DIRS})

# bundle all thirdparty stuff as a static library
FILE(GLOB_RECURSE THIRDPARTY thirdparty/*.cpp)
add_library(thirdparty STATIC ${THIRDPARTY})

set(LINK_LIBS ${LINK_LIBS} thirdparty ${PCL_COMMON_LIBRARIES} ${PCL_SEARCH_LIBRARIES} ${PCL_FEATURES_LIBRARIES})

# build a library from the masbpcpp processing functions
# add_library(masbcpp STATIC src/compute_ma_processing.cpp src/compute_normals_processing.cpp src/simplify_processing.cpp)
add_library(masbcpp STATIC src/io.cpp src/compute_normals_processing.cpp src/compute_ma_processing.cpp src/simplify_processing.cpp)

# set excutables
add_executable(compute_ma src/compute_ma.cpp)
add_executable(compute_normals src/compute_normals.cpp)
add_executable(simplify src/simplify.cpp)

# link targets
target_link_libraries(masbcpp ${LINK_LIBS})

target_link_libraries(compute_ma masbcpp)
target_link_libraries(compute_normals masbcpp)
target_link_libraries(simplify masbcpp)

# install(TARGETS compute_ma compute_normals simplify DESTINATION bin)
