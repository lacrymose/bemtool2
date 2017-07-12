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
# Boost package - special functions
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message("Boost not found, no access to special functions")
endif()


# Htool
message("-- Htool and its dependencies will be built as submodules")
execute_process(COMMAND git submodule update --init --remote external/htool
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(external/htool)

#=============================================================================#
#=========================== Libraries =======================================#
#=============================================================================#
#=== BEMTOOL
if(Boost_FOUND)
	message("-- Bemtool and dependencies found")
	add_library(bemtool INTERFACE)
	target_include_directories(bemtool INTERFACE
						${bemtool_include_dir}
						${Boost_INCLUDE_DIRS})
  target_link_libraries(bemtool INTERFACE htool)

else()
	message("Bemtool and dependencies not found")
endif()
#===



#=============================================================================#
#========================== Repertories ======================================#
#=============================================================================#

# Add documentation
# add_subdirectory(doc)

# Add test files
# add_subdirectory(test)

# add_subdirectory(externals)
