



# GMM package - linear algebra
find_package(GMM REQUIRED)
include_directories(${GMM_INCLUDE_DIR})

# Boost package - special functions
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include") 

add_library(bemtool "triangle_dunavant_rule.cpp")
target_include_directories(	bemtool PUBLIC "${bemtool_dir}"
					PRIVATE GMM Boost)
