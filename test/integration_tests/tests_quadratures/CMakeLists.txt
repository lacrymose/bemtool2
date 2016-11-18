# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_quadratures processing.cpp
                    tests_quadratures.cpp tests_quadratures.hpp "${bemtool_dir}")