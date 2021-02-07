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

# header files for Lapacke
find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h PATHS /usr/include/lapacke)
if(LAPACKE_INCLUDE_DIR)
  include_directories(${LAPACKE_INCLUDE_DIR})
else(LAPACKE_INCLUDE_DIR) 
  message(FATAL_ERROR "Lapacke header file lapacke.h not found!")
endif(LAPACKE_INCLUDE_DIR)
# header files for CBLAS
find_path(CBLAS_INCLUDE_DIR NAMES cblas.h PATHS /usr/include/cblas)
if(CBLAS_INCLUDE_DIR)
  include_directories(${CBLAS_INCLUDE_DIR})
else(CBLAS_INCLUDE_DIR)
message(FATAL_ERROR "CBLAS header file cblas.h not found!")
endif(CBLAS_INCLUDE_DIR)
# header files for GSL
find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL header files not found")
endif(GSL_FOUND)

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

# link libraries of Lapacke
find_library(LAPACKE_LIB lapacke)
if(LAPACKE_LIB)
  target_link_libraries(cali PUBLIC ${LAPACKE_LIB})
else(LAPACKE_LIB) 
  message(FATAL_ERROR "Lapacke library liblapacke not found!")
endif(LAPACKE_LIB)
# link libraries of CBLAS
find_library(CBLAS_LIB cblas)
if(CBLAS_LIB)
  include_directories(${CBLAS_INCLUDE_DIR})
  target_link_libraries(cali PUBLIC ${CBLAS_LIB})
else(CBLAS_LIB)
message(FATAL_ERROR "CBLAS library libcblas not found!")
endif(CBLAS_LIB)
# link libraries of GSL
find_package(GSL REQUIRED)
if(GSL_FOUND)
  target_link_libraries(cali PUBLIC ${GSL_LIBRARY})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL library not found")
endif(GSL_FOUND)

# put this at the end to ensure configurations being passed to Python builder 
add_subdirectory("src")