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

# Add documentation
# add_subdirectory(doc)

# Add test files
add_subdirectory(test)
