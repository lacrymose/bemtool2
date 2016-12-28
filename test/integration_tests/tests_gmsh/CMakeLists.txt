# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_gmsh processing.cpp tests_gmsh.cpp)
