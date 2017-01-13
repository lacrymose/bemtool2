# GMM package - linear algebra
find_package(GMM REQUIRED)
if(NOT GMM_FOUND)
    message("GMM not found, linear algebra not possible")
endif()

# Boost package - special functions
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message("Boost not found, no access to special functions")
endif()

# BLAS
find_package(BLAS REQUIRED)

# LAPACK
find_package(LAPACK REQUIRED)

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_hpddm tests_hpddm.cpp)

target_compile_options(Tests_hpddm PRIVATE -DHPDDM_MPI=0)

target_include_directories(Tests_hpddm PRIVATE 
								${bemtool_include_dir} 
								${GMM_INCLUDE_DIR}
								${Boost_INCLUDE_DIRS}
								${hpddm_include_dir})

target_link_libraries(Tests_hpddm 
						${BLAS_LIBRARIES}
						${LAPACK_LIBRARIES})

