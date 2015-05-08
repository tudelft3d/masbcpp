cmake_minimum_required(VERSION 3.2)

project(masbcpp)

set(CMAKE_BUILD_TYPE Release)

set(LINK_LIBS )
set(COMPILE_OPTIONS -funroll-loops -ffast-math)

if(EXISTS "/usr/local/bin/clang-omp++")
	set(CMAKE_CXX_COMPILER "/usr/local/bin/clang-omp++")
	add_definitions(-DWITH_OPENMP)

	set(LINK_LIBS ${LINK_LIBS} -liomp5)
	set(COMPILE_OPTIONS ${COMPILE_OPTIONS} -fopenmp)
endif()

FILE(GLOB_RECURSE THIRDPARTY thirdparty/*.cpp)

add_executable(compute_ma src/compute_ma.cpp ${THIRDPARTY})

set_target_properties(compute_ma PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(compute_ma PRIVATE "${CMAKE_SOURCE_DIR}/thirdparty/vrui")
target_link_libraries(compute_ma ${LINK_LIBS})
target_compile_options(compute_ma PRIVATE "${COMPILE_OPTIONS}")

install(TARGETS compute_ma DESTINATION bin)