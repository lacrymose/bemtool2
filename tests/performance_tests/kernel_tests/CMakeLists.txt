#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_SLP_PH_2D tests_SLP_PH_2D.cpp)
target_link_libraries(Tests_SLP_PH_2D bemtool)

add_executable(Tests_HSP_PH_2D tests_HSP_PH_2D.cpp)
target_link_libraries(Tests_HSP_PH_2D bemtool)

add_executable(Tests_TDLP_PH_2D tests_TDLP_PH_2D.cpp)
target_link_libraries(Tests_TDLP_PH_2D bemtool)


add_executable(Tests_DLP_PH_2D tests_DLP_PH_2D.cpp)
target_link_libraries(Tests_DLP_PH_2D bemtool)
