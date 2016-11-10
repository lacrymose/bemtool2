add_subdirectory(tests_circle_2D)

# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include")
include_directories(${bemtool_dir})

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
# add_executable(Tests_quadrature tests_quadrature.cpp)