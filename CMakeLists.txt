cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(cali LANGUAGES CXX C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SRC "./")

find_path(CDNEST_INCLUDE_DIR NAMES dnestvars.h PATHS /home/liyropt/Projects/GIT/CDNestSingle)
find_library(CDNEST_LIB NAMES dnest PATHS /home/liyropt/Projects/GIT/CDNestSingle)
include_directories(${CDNEST_INCLUDE_DIR})

find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h PATHS /usr/include/lapacke)
find_library(LAPACKE_LIB lapacke)
include_directories(${LAPACKE_INCLUDE_DIR})

find_path(CBLAS_INCLUDE_DIR NAMES cblas.h PATHS /usr/include/cblas)
find_library(CBLAS_LIB cblas)
include_directories(${CBLAS_INCLUDE_DIR})

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

target_link_libraries(cali utilities ${CDNEST_LIB} ${LAPACKE_LIB} ${CBLAS_LIB})