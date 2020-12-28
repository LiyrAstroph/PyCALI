cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(pycali LANGUAGES CXX C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
# when installing directly in CMake, redirect its intall root to dist dir
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/dist" CACHE PATH "default install path" FORCE )
endif()

add_subdirectory("src")

set(SRC "src/pycali/pycali")
set(SRC_CDNEST "src/pycali/cdnest")
add_executable(cali 
    ${SRC}/main.cpp 
    ${SRC}/utilities.hpp
    ${SRC}/utilities.cpp
    ${SRC}/mathfun.h 
    ${SRC}/mathfun.c
    ${SRC_CDNEST}/dnest.c
    ${SRC_CDNEST}/dnestpostprocess.c
    ${SRC_CDNEST}/dnestvars.h
    ${SRC_CDNEST}/dnestvars.c
    )

find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h PATHS /usr/include/lapacke)
find_library(LAPACKE_LIB lapacke)

if(LAPACKE_INCLUDE_DIR AND LAPACKE_LIB)
  include_directories(${LAPACKE_INCLUDE_DIR})
  target_link_libraries(cali PUBLIC ${LAPACKE_LIB})
else(LAPACKE_INCLUDE_DIR AND LAPACKE_LIB) 
  message(FATAL_ERROR "Lapacke library not found!")
endif(LAPACKE_INCLUDE_DIR AND LAPACKE_LIB)

find_path(CBLAS_INCLUDE_DIR NAMES cblas.h PATHS /usr/include/cblas)
find_library(CBLAS_LIB cblas)
if(CBLAS_INCLUDE_DIR AND CBLAS_LIB)
  include_directories(${CBLAS_INCLUDE_DIR})
  target_link_libraries(cali PUBLIC ${CBLAS_LIB})
else(CBLAS_INCLUDE_DIR AND CBLAS_LIB)
message(FATAL_ERROR "CBLAS library not found!")
endif(CBLAS_INCLUDE_DIR AND CBLAS_LIB)

find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(cali PUBLIC ${GSL_LIBRARY})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL library not found")
endif(GSL_FOUND)
