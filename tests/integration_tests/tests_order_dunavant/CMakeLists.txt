#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#

add_executable(Tests_order_dunavant processing.cpp tests_sphere_3D.cpp tests_sphere_3D.hpp)
target_link_libraries(Tests_order_dunavant bemtool) 
