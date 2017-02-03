#=============================================================================#
#============================== Project ======================================#
#=============================================================================#
cmake_minimum_required(VERSION 3.0)
Project(Bemtool)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release)
endif()

# Usefule variable
set(bemtool_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/include")

# C++ 11 required
if ("${CMAKE_VERSION}" VERSION_LESS 3.1)
	add_compile_options(-std=c++11)
else ()
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

# Files to do find_package for some module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")

#=============================================================================#
#========================== External Libraries ===============================#
#=============================================================================#
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

# Eigen package
# find_package (Eigen3)


#=== MPI
find_package(MPI)
if(NOT MPI_FOUND)
    message("MPI not found. Needed in htool and hpddm")
else()
	add_library(mpi INTERFACE)
	target_include_directories(mpi INTERFACE ${MPI_INCLUDE_PATH})
endif()
#===

#=== BLAS
find_package(BLAS)
if(NOT BLAS_FOUND)
    message("Blas not found. Needed in htool and hpddm")
else()
	add_library(blas INTERFACE)
	target_link_libraries(blas INTERFACE ${BLAS_LIBRARIES})
endif()
#===

#=== LAPACK
find_package(LAPACK)
if(NOT LAPACK_FOUND)
    message("Lapack not found. Needed in htool and hpddm")
else()
	add_library(lapack INTERFACE)
	target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES})
endif()
#===

#=============================================================================#
#=========================== Libraries =======================================#
#=============================================================================#
#=== BEMTOOL
if(Boost_FOUND AND GMM_FOUND)
	message("-- Bemtool and dependencies found")
	add_library(bemtool INTERFACE)
	target_include_directories(bemtool INTERFACE 
						${bemtool_include_dir}
						${Boost_INCLUDE_DIRS}
						${GMM_INCLUDE_DIR})
else()
	message("Bemtool and dependencies not found")
endif()
#===

#=== HPDDM
find_path(hpddm_include_dir "HPDDM.hpp" "${CMAKE_CURRENT_SOURCE_DIR}/external/hpddm/include")
if(hpddm_include_dir AND MPI_FOUND AND BLAS_FOUND AND LAPACK_FOUND)
	message("-- HPDDM and dependencies found")
elseif(NOT hpddm_include_dir)
	execute_process(COMMAND git submodule update --init --recursive
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set(hpddm_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/external/hpddm/include")
endif()
add_library(hpddm INTERFACE)
target_include_directories(hpddm INTERFACE "${hpddm_include_dir}")
target_link_libraries(hpddm INTERFACE
					blas
					lapack)
#===



#=== Htool and dependencies
message("-- Htool and its dependencies will be built as submodules")
execute_process(COMMAND git submodule update --init --recursive
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
                    
    
# Configuration options for htool
set(HTOOL_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)

# Add the configurations from htool
add_subdirectory(external/htool)
#===

#=============================================================================#
#========================== Repertories ======================================#
#=============================================================================#

# Add documentation
# add_subdirectory(doc)

# Add test files
add_subdirectory(test)

# add_subdirectory(externals)
