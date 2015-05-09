cmake_minimum_required(VERSION 2.8)

project(masbcpp)

set(CMAKE_BUILD_TYPE Release)

set(LINK_LIBS )
set(COMPILE_OPTIONS "-funroll-loops -ffast-math")

find_package(OpenMP)

if(OPENMP_FOUND)
	add_definitions(-DWITH_OPENMP)
	set(COMPILE_OPTIONS "${COMPILE_OPTIONS} ${OpenMP_CXX_FLAGS}")
endif()

if(EXISTS "/usr/local/bin/clang-omp++")
	set(CMAKE_CXX_COMPILER "/usr/local/bin/clang-omp++")
	add_definitions(-DWITH_OPENMP)
	add_definitions(-DCLANG_OMP)

	set(LINK_LIBS ${LINK_LIBS} -liomp5)
	set(COMPILE_OPTIONS "${COMPILE_OPTIONS} -fopenmp")
endif()

FILE(GLOB_RECURSE THIRDPARTY thirdparty/*.cpp)

add_executable(compute_ma src/compute_ma.cpp ${THIRDPARTY})

set_target_properties(compute_ma PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(compute_ma PRIVATE ${CMAKE_SOURCE_DIR}/thirdparty ${CMAKE_SOURCE_DIR}/thirdparty/vrui)
target_link_libraries(compute_ma ${LINK_LIBS})
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COMPILE_OPTIONS}" )
#target_compile_options(compute_ma PRIVATE "${COMPILE_OPTIONS}")

install(TARGETS compute_ma DESTINATION bin)