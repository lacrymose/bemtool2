#=============================================================================#
#============================== Project ======================================#
#=============================================================================#
cmake_minimum_required(VERSION 3.0)
Project(Bemtool)

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

#=== HTOOL
# find_path(htool_dir "${CMAKE_CURRENT_SOURCE_DIR}/externals/htool")
# if(htool_dir)
# 	add_library(htool INTERFACE)
# 	target_include_directories(htool INTERFACE ${htool_include_dir})
# 	add_subdirectory(externals)
# else()
# 	... # read on
# endif()

#=== HPDDM
find_path(hpddm_include_dir "HPDDM.hpp" "${CMAKE_CURRENT_SOURCE_DIR}/externals/hpddm/include")
if(hpddm_include_dir AND MPI_FOUND AND BLAS_FOUND AND LAPACK_FOUND)
	message("-- HPDDM and dependencies found")
elseif(NOT hpddm_include_dir)
	execute_process(COMMAND git submodule init
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	execute_process(COMMAND git submodule update
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set(hpddm_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/external/hpddm/include")
endif()
add_library(hpddm INTERFACE)
target_include_directories(hpddm INTERFACE "${hpddm_include_dir}")
target_link_libraries(hpddm INTERFACE
					blas
					lapack)
#===
					
					
					
# get_property(test TARGET hpddm PROPERTY INTERFACE_LINK_LIBRARIES)
# message(${test})

# #OpenGL - htool/hpddm
# find_package(OpenGL)
# if(NOT OPENGL_FOUND)
#     message(ERROR " OPENGL not found!")
# endif(NOT OPENGL_FOUND)
# 
# #FreeGlut - htool/hpddm
# find_package(GLUT)
# if(NOT GLUT_FOUND)
#     message(ERROR " GLUT not found!")
# endif(NOT GLUT_FOUND)
# 
# configure_file(${htool_include_dir}/config.h.in ${CMAKE_BINARY_DIR}/generated/config.h )
# include_directories( ${CMAKE_BINARY_DIR}/generated/ )


#=============================================================================#
#========================== Repertories ======================================#
#=============================================================================#

# Add documentation
# add_subdirectory(doc)

# Add test files
add_subdirectory(test)

# add_subdirectory(externals)
