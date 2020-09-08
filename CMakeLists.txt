cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(cali LANGUAGES CXX C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SRC "./")

add_subdirectory("./cdnest")

add_executable(cali 
    ${SRC}/main.cpp 
    )

add_library(utilities
    SHARED
    ${SRC}/utilities.hpp
    ${SRC}/utilities.cpp
    ${SRC}/mathfun.h 
    ${SRC}/mathfun.c
)

find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h PATHS /usr/include/lapacke)
find_library(LAPACKE_LIB lapacke)
include_directories(${LAPACKE_INCLUDE_DIR})

find_path(CBLAS_INCLUDE_DIR NAMES cblas.h PATHS /usr/include/cblas)
find_library(CBLAS_LIB cblas)
include_directories(${CBLAS_INCLUDE_DIR})

target_link_libraries(cali utilities dnest ${LAPACKE_LIB} ${CBLAS_LIB})

find_package(pybind11 REQUIRED)

add_library(pycali MODULE wrapper.cpp)

target_link_libraries(pycali PUBLIC utilities dnest ${LAPACKE_LIB} ${CBLAS_LIB} PRIVATE pybind11::module)
set_target_properties(pycali PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}")