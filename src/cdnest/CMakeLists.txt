
add_library(dnest
  SHARED
  dnest.c
  dnestpostprocess.c
  dnestvars.h
  dnestvars.c
)

find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(dnest ${GSL_LIBRARY})
else(GSL_FOUND)
  message(FATAL_ERROR "GSL library not found")
endif(GSL_FOUND)
