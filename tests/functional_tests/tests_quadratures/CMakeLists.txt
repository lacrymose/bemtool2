#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_quadratures tests_quadratures.cpp)
target_link_libraries(Tests_quadratures bemtool)


add_executable(Tests_potential_2D tests_potential_2D.cpp)
target_link_libraries(Tests_potential_2D bemtool)

# add_executable(Tests_circle_2D tests_circle_2D.cpp)
# target_link_libraries(Tests_circle_2D bemtool)

add_executable(Tests_overlapp_2D tests_overlapp_2D.cpp)
target_link_libraries(Tests_overlapp_2D bemtool)

add_executable(Tests_overlapp_3D tests_overlapp_3D.cpp)
target_link_libraries(Tests_overlapp_3D bemtool)
