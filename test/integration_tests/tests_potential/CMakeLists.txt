#=============================================================================#
#=========================== Libraries =======================================#
#=============================================================================#
# Boost package - special functions
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message("Boost not found, no access to special functions")
endif()

# Eigen package
find_package (Eigen3 REQUIRED)

# MPI
find_package (MPI REQUIRED)

# BLAS
find_package(BLAS REQUIRED)

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_potential processing.cpp
                     tests_potential.hpp)

target_include_directories(Tests_potential PRIVATE 
								${bemtool_include_dir}
								${Boost_INCLUDE_DIRS}
								${htool_include_dir}
								${MPI_INCLUDE_PATH}
								${EIGEN3_INCLUDE_DIR})
							
							
target_link_libraries (Tests_potential ${MPI_LIBRARIES}
								${BLAS_LIBRARIES})
