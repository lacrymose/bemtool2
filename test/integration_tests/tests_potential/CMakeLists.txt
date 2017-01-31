#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_potential tests_potential.cpp)

target_include_directories(Tests_potential PRIVATE 
								${bemtool_include_dir}
								${Boost_INCLUDE_DIRS})

