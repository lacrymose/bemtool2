#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_quadratures tests_quadratures.cpp)

target_include_directories(Tests_quadratures PRIVATE 
								${bemtool_include_dir})
