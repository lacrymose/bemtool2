#=============================================================================#
#=========================== Libraries =======================================#
#=============================================================================#
# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
# include_directories(${bemtool_dir})

# Boost package - special functions
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message("Boost not found, no access to special functions")
endif()

# HTOOL package
set(htool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../../htool/include")
# include_directories(${htool_dir})

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
								${bemtool_dir}
								${Boost_INCLUDE_DIRS})
							
                     
# target_include_directories(Tests_potential PRIVATE 
# 								${bemtool_dir} 
# 								${htool_dir}
# 								${EIGEN3_INCLUDE_DIR}
# 								${MPI_INCLUDE_PATH})
# 								
# target_link_libraries (Tests_potential ${MPI_LIBRARIES})
