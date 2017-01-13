# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_quadratures tests_quadratures.cpp)

target_include_directories(Tests_quadratures PRIVATE 
								${bemtool_dir})
