#=============================================================================#
#============================== Project ======================================#
#=============================================================================#
cmake_minimum_required(VERSION 3.0)
Project(Bemtool)

# Usefule variable
set(bemtool_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(htool_include_dir   "${CMAKE_CURRENT_SOURCE_DIR}/libs/htool/include")
set(hpddm_include_dir   "${CMAKE_CURRENT_SOURCE_DIR}/libs/hpddm/include")

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
#=========================== Libraries =======================================#
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
find_package (Eigen3)

# MPI - htool/hpddm
find_package (MPI)
if(NOT MPI_FOUND)
    message("MPI not found. Needed in htool and hpddm")
endif()

# BLAS - htool/hpddm
find_package(BLAS)
if(NOT BlAS_FOUND)
    message("Blas not found. Needed in htool and hpddm")
endif()

# LAPACK - hpddm
find_package(LAPACK REQUIRED)
if(NOT LAPACK_FOUND)
    message("Lapack not found. Needed in htool and hpddm")
endif()


#OpenGL - htool/hpddm
find_package(OpenGL)
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

#FreeGlut - htool/hpddm
find_package(GLUT)
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

configure_file(${htool_include_dir}/config.h.in ${CMAKE_BINARY_DIR}/generated/config.h )
include_directories( ${CMAKE_BINARY_DIR}/generated/ )


#=============================================================================#
#========================== Repertories ======================================#
#=============================================================================#

# Add documentation
# add_subdirectory(doc)

# Add test files
add_subdirectory(test)
