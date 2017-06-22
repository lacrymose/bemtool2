#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_hpddm tests_hpddm.cpp)

target_compile_options(Tests_hpddm PRIVATE -DHPDDM_MPI=0)

target_link_libraries(Tests_hpddm bemtool hpddm)

