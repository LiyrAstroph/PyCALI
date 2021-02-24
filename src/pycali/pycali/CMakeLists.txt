# src folder name is used as the target name 
get_filename_component(TARGET ${CMAKE_CURRENT_SOURCE_DIR} NAME)

include_directories(${LAPACKE_INCLUDE_DIR})
include_directories(${CBLAS_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIRS})

find_package(pybind11 REQUIRED)
# create the target
pybind11_add_module(${TARGET} 
    wrapper.cpp 
    utilities.hpp
    utilities.cpp
    mathfun.h 
    mathfun.c
    ../cdnest/dnest.c
    ../cdnest/dnestpostprocess.c
    ../cdnest/dnestvars.h
    ../cdnest/dnestvars.c
    )
#add_library(${TARGET} MODULE ${SRC}/wrapper.cpp)

target_link_libraries(${TARGET} PUBLIC ${LAPACKE_LIB} ${LAPACK_LIB} ${CBLAS_LIB} ${GSL_LIBRARY} PRIVATE pybind11::module)
set_target_properties(${TARGET} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}")

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
