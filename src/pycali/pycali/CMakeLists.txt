# src folder name is used as the target name 
get_filename_component(TARGET ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(SRC ".")

find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h PATHS /usr/include/lapacke)
find_library(LAPACKE_LIB lapacke)
include_directories(${LAPACKE_INCLUDE_DIR})

find_path(CBLAS_INCLUDE_DIR NAMES cblas.h PATHS /usr/include/cblas)
find_library(CBLAS_LIB cblas)
include_directories(${CBLAS_INCLUDE_DIR})

find_package(pybind11 REQUIRED)

# create the target
pybind11_add_module(${TARGET} 
    ${SRC}/wrapper.cpp 
    ${SRC}/utilities.hpp
    ${SRC}/utilities.cpp
    ${SRC}/mathfun.h 
    ${SRC}/mathfun.c
    ../cdnest/dnest.c
    ../cdnest/dnestpostprocess.c
    ../cdnest/dnestvars.h
    ../cdnest/dnestvars.c
    )
#add_library(${TARGET} MODULE ${SRC}/wrapper.cpp)

target_link_libraries(${TARGET} PUBLIC ${LAPACKE_LIB} ${CBLAS_LIB} PRIVATE pybind11::module)
set_target_properties(${TARGET} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}")

find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(${TARGET} PUBLIC ${GSL_LIBRARY})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL library not found")
endif(GSL_FOUND)

# install pyd files one directory structure UP in the destination as 
# CMAKE_CURRENT_SOURCE_DIR is holding the source files of the MODULE
# which matches the name of CMAKE_CURRENT_SOURCE_DIR
get_filename_component(CURRENT_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
if(${SRC_DIR} STREQUAL ${CURRENT_SRC_DIR}) 
    set(DST_DIR ".") # if parent is the base folder
else()
    file(RELATIVE_PATH DST_DIR ${SRC_DIR} ${CURRENT_SRC_DIR})
endif()
install(TARGETS ${TARGET} DESTINATION ${DST_DIR} COMPONENT "EXT")
