#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#

add_executable(Tests_order_dunavant processing.cpp tests_sphere_3D.cpp tests_sphere_3D.hpp)
target_include_directories(Tests_order_dunavant PRIVATE ${GMM_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${bemtool_include_dir}) 
