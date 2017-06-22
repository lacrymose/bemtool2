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
# include_directories(${EIGEN3_INCLUDE_DIR})

find_package (MPI REQUIRED)
# include_directories(${MPI_INCLUDE_PATH})

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_potential processing.cpp
                     tests_potential.hpp)

target_include_directories(Tests_potential PRIVATE 
								${bemtool_include_dir}
								${Boost_INCLUDE_DIRS})
							
                     
# target_include_directories(Tests_potential PRIVATE 
# 								${bemtool_dir} 
# 								${htool_dir}
# 								${EIGEN3_INCLUDE_DIR}
# 								${MPI_INCLUDE_PATH})
# 								
# target_link_libraries (Tests_potential ${MPI_LIBRARIES})
