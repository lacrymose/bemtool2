#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_htool tests_htool.cpp)

target_include_directories(Tests_htool PRIVATE 
								${bemtool_include_dir}
								${Boost_INCLUDE_DIRS}
								${htool_include_dir}
								${MPI_INCLUDE_PATH}
								${EIGEN3_INCLUDE_DIR})
								
target_link_libraries (Tests_htool ${MPI_LIBRARIES}
								${BLAS_LIBRARIES})
								