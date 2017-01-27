#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_circle_2D processing.cpp
                    tests_circle_2D.cpp)

target_include_directories(Tests_circle_2D PRIVATE
							${bemtool_include_dir}
							${Boost_INCLUDE_DIRS}
							${GMM_INCLUDE_DIR})