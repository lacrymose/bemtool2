
#=============================================================================#
#============================= Packages ======================================#
#=============================================================================#
# GMM package - linear algebra
find_package(GMM REQUIRED)
if(NOT GMM_FOUND)
  message("GMM not found, linear algebra not possible")
endif()

# Boost package - special functions
# find_path(BOOST math/special_functions/bessel.hpp PATH_SUFFIXES boost)
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message("Boost not found, no access to special functions")
endif()

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#

add_executable(Tests_order_dunavant processing.cpp tests_sphere_3D.cpp tests_sphere_3D.hpp)
target_include_directories(Tests_order_dunavant PRIVATE ${GMM_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${bemtool_include_dir}) 
