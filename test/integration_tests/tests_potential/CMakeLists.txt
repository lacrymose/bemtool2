# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_potential processing.cpp
                    tests_potential.cpp tests_potential.hpp)
