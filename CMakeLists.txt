#=============================================================================#
#============================== Project ======================================#
#=============================================================================#
cmake_minimum_required(VERSION 3.0)
Project(Bemtool)


# To force c++11
if (${CMAKE_VERSION} VERSION_LESS 3.1)
	add_compile_options(-std=c++11)
elseif (${CMAKE_VERSION} VERSION_LESS 3.6.3 AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
	add_compile_options(-std=c++11)
else()
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release)
endif()

# Usefule variable
set(bemtool_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Files to do find_package for some module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")

# Options
option(BEMTOOL_WITH_EXAMPLES "Build htool examples ?" ON)
option(BEMTOOL_WITH_TESTS    "Build htool tests ?" ON)


#=============================================================================#
#========================== External Libraries ===============================#
#=============================================================================#
# GMM package - linear algebra
# find_package(GMM REQUIRED)
# if(NOT GMM_FOUND)
#   message("GMM not found, linear algebra not possible")
# endif()
# set (GMM_INCLUDE_DIR "external" )
# set (GMM_FOUND 1)


# EIGEN3
find_package(EIGEN3 REQUIRED)

# Boost package - special functions
find_package(Boost REQUIRED)

# Htool
find_package(Htool QUIET)

# HPDDM
find_package(HPDDM QUIET)

# MPI
find_package (MPI REQUIRED QUIET)
message("-- MPI libraries found in " "${MPI_LIBRARIES}")
message("-- MPI include files found in " "${MPI_INCLUDE_PATH}")

# BLAS
find_package(BLAS REQUIRED QUIET)
message("-- Found Blas implementation:" "${BLAS_LIBRARIES}")

# LAPACK
find_package(LAPACK REQUIRED QUIET)
message("-- Found Lapack:" "${LAPACK_LIBRARIES}")

#MUMPS
find_package(MUMPS)


#=============================================================================#
#=========================== Libraries =======================================#
#=============================================================================#
#=== BEMTOOL
add_library(bemtool INTERFACE)
target_include_directories(bemtool INTERFACE
					${bemtool_include_dir}
					${Boost_INCLUDE_DIRS}
					${HTOOL_INCLUDE_DIRS}
					${MPI_INCLUDE_PATH}
					${HPDDM_INCLUDE_DIRS}
					${MUMPS_INCLUDE_DIRS}
				)
target_link_libraries(bemtool INTERFACE
					${MPI_LIBRARIES}
					${BLAS_LIBRARIES}
					${LAPACK_LIBRARIES}
					${MUMPS_LIBRARIES}
				)
#===



#=============================================================================#
#========================== Repertories ======================================#
#=============================================================================#

if (BEMTOOL_WITH_EXAMPLES)
	enable_testing()
	add_custom_target(build-examples)
	# add_subdirectory(examples EXCLUDE_FROM_ALL)
endif()
if (BEMTOOL_WITH_TESTS)
	enable_testing()
	add_custom_target(build-tests)
	add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()



# Add documentation
# add_subdirectory(doc)

# Add test files
# add_subdirectory(test)

# add_subdirectory(externals)
