cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(pycali LANGUAGES CXX C Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O3")

set(CMAKE_Fortran_FLAGS "-frecursive")
set(CMAKE_Fortran_FLAGS_DEBUG "-g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

##if using clang compiler, add c++11/c++14/c++17... support (https://clang.llvm.org/cxx_status.html)
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  message(${CMAKE_CXX_COMPILER_ID})
  set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "-std=c++11 ${CMAKE_C_FLAGS}")
endif()

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
# when installing directly in CMake, redirect its intall root to dist dir
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/dist" CACHE PATH "default install path" FORCE )
endif()

set(SRC "src/pycali/pycali")
set(SRC_CDNEST "src/pycali/cdnest")
set(SRC_LAPACKE "src/pycali/lapacke")
set(SRC_CBLAS "src/pycali/cblas")

# header files for CBLAS
find_path(CBLAS_INCLUDE_DIR NAMES cblas.h PATHS /usr/include/cblas)
if(CBLAS_INCLUDE_DIR)
  include_directories(${CBLAS_INCLUDE_DIR})
else(CBLAS_INCLUDE_DIR)
  message("CBLAS header file cblas.h not found, use the source code!")
  set(CBLAS_INCLUDE_DIR "${SRC_CBLAS}")
  include_directories(${CBLAS_INCLUDE_DIR})
endif(CBLAS_INCLUDE_DIR)
# header files for Lapacke
find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h PATHS /usr/include/lapacke)
if(LAPACKE_INCLUDE_DIR)
  include_directories(${LAPACKE_INCLUDE_DIR})
else(LAPACKE_INCLUDE_DIR) 
  message("Lapacke header file lapacke.h not found, use the source code!")
  set(LAPACKE_INCLUDE_DIR "${SRC_LAPACKE}")
  include_directories(${LAPACKE_INCLUDE_DIR})
endif(LAPACKE_INCLUDE_DIR)
# header files for GSL
find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL header files not found")
endif(GSL_FOUND)

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

# source files for CBLAS
file(READ "${SRC_CBLAS}/cblas_srcs.txt" SOURCES)
STRING(REGEX REPLACE ";" "\\\\;" SOURCES "${SOURCES}")
STRING(REGEX REPLACE "\n" ";" SOURCES "${SOURCES}")
FOREACH(VAR ${SOURCES})
  string(STRIP "${VAR}" VARSTRIP)
  list(APPEND CBLAS_SRCS "${SRC_DIR}/${SRC_CBLAS}/${VARSTRIP}")
ENDFOREACH()

# source files for LAPACKE
file(READ "${SRC_LAPACKE}/lapacke_srcs.txt" SOURCES)
STRING(REGEX REPLACE ";" "\\\\;" SOURCES "${SOURCES}")
STRING(REGEX REPLACE "\n" ";" SOURCES "${SOURCES}")
FOREACH(VAR ${SOURCES})
  string(STRIP "${VAR}" VARSTRIP)
  list(APPEND LAPACKE_SRCS "${SRC_DIR}/${SRC_LAPACKE}/${VARSTRIP}")
ENDFOREACH()

# link libraries of CBLAS
find_library(BLAS_LIB blas)
find_library(CBLAS_LIB cblas)
if(CBLAS_LIB)
  include_directories(${CBLAS_INCLUDE_DIR})
  target_link_libraries(cali PUBLIC ${BLAS_LIB} ${CBLAS_LIB})
else(CBLAS_LIB)
  message(WARNING "CBLAS library libcblas not found!\n"
                  "Compile from source code.")
  target_sources(cali PUBLIC ${CBLAS_SRCS})
endif(CBLAS_LIB)

# link libraries of Lapacke
find_library(LAPACK_LIB lapack)
find_library(LAPACKE_LIB lapacke)
if(LAPACKE_LIB)
  target_link_libraries(cali PUBLIC ${LAPACKE_LIB} ${LAPACK_LIB})
else(LAPACKE_LIB) 
  message(WARNING "Lapack libraries liblapack and liblapacke not found!\n"
                  "Compile from source code.")
  target_sources(cali PUBLIC ${LAPACKE_SRCS})
endif(LAPACKE_LIB)

# link libraries of GSL
find_package(GSL REQUIRED)
if(GSL_FOUND)
  target_link_libraries(cali PUBLIC ${GSL_LIBRARY})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL library not found.")
endif(GSL_FOUND)

# put this at the end to ensure configurations being passed to Python builder 
add_subdirectory("src")