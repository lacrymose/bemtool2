# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

# Boost package - special functions
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
message("Boost not found, no access to special functions")
else()
include_directories(${Boost_INCLUDE_DIRS})
endif()

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_gmsh processing.cpp tests_gmsh.cpp)
