cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(cali LANGUAGES CXX C)

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

set(SRC "./src/")
set(SRC_CDNEST "./src/cdnest")

# create cali
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

# link Lapacke
find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h PATHS /usr/include/lapacke)
find_library(LAPACKE_LIB lapacke)
if(LAPACKE_LIB AND LAPACKE_INCLUDE_DIR)
  include_directories(${LAPACKE_INCLUDE_DIR})
  target_link_libraries(cali ${LAPACKE_LIB})
else(LAPACKE_LIB AND LAPACKE_INCLUDE_DIR)
  message(FATAL_ERROR "Lapacke library not found")
endif(LAPACKE_LIB AND LAPACKE_INCLUDE_DIR)

# link CBLAS
find_path(CBLAS_INCLUDE_DIR NAMES cblas.h PATHS /usr/include/cblas)
find_library(CBLAS_LIB cblas)
if(CBLAS_LIB AND CBLAS_INCLUDE_DIR)
  include_directories(${CBLAS_INCLUDE_DIR})
  target_link_libraries(cali ${CBLAS_LIB})
else(CBLAS_LIB AND CBLAS_INCLUDE_DIR)
  message(FATAL_ERROR "CBLAS library not found")
endif(CBLAS_LIB AND CBLAS_INCLUDE_DIR)

# link GSL
find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(cali ${GSL_LIBRARY})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL library not found")
endif(GSL_FOUND)

install(TARGETS cali)

# create pycali
find_package(pybind11 REQUIRED)
pybind11_add_module(pycali
    ${SRC}/wrapper.cpp 
    ${SRC}/utilities.hpp
    ${SRC}/utilities.cpp
    ${SRC}/mathfun.h 
    ${SRC}/mathfun.c
    ${SRC_CDNEST}/dnest.c
    ${SRC_CDNEST}/dnestpostprocess.c
    ${SRC_CDNEST}/dnestvars.h
    ${SRC_CDNEST}/dnestvars.c
    )

target_link_libraries(pycali PUBLIC ${LAPACKE_LIB} ${CBLAS_LIB} PRIVATE pybind11::module)

find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(pycali PUBLIC ${GSL_LIBRARY})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL library not found")
endif(GSL_FOUND)

set_target_properties(pycali PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}")