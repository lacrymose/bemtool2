
#=============================================================================#
#============================= Packages ======================================#
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

# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")


#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_circle_2D processing.cpp
                    tests_circle_2D.cpp)

target_include_directories(Tests_circle_2D PRIVATE
							${bemtool_dir}
							${Boost_INCLUDE_DIRS}
							${GMM_INCLUDE_DIR})