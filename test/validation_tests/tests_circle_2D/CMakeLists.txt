
#=============================================================================#
#============================= Packages ======================================#
#=============================================================================#
# GMM package - linear algebra
find_package(GMM REQUIRED)
if(NOT GMM_FOUND)
    message("GMM not found, linear algebra not possible")
else()
	include_directories(${GMM_INCLUDE_DIR})
endif()

# Boost package - special functions
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message("Boost not found, no access to special functions")
else()
include_directories(${Boost_INCLUDE_DIRS})
endif()

# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_circle_2D processing.cpp
                    tests_circle_2D.cpp)
