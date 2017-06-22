#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_htool tests_htool.cpp)
target_link_libraries (Tests_htool bemtool htool)

add_executable(Tests_htool_gui tests_htool_gui.cpp)
target_link_libraries (Tests_htool_gui bemtool htool htool_visu)
target_compile_options(Tests_htool_gui PUBLIC -DHTOOL_WITH_GUI)