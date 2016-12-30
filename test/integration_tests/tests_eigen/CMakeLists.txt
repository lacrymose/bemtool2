# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

# Eigen package
find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
message ("${EIGEN3_INCLUDE_DIR}")
#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_eigen processing.cpp tests_eigen.cpp)
