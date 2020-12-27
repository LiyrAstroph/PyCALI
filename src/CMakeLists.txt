# root source directory should only contain this file (unless creating a root-package Python module)

# variable to hold the path to src folder to compute relative path of its decendents
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# traverse the tree
add_subdirectory("pycali")
