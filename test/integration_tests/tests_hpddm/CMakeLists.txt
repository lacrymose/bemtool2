# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

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

# HPDDM
set(hpddm_include_dir /Users/pierremarchand/Partage/phd/lib/hpddm/include)
include_directories(${hpddm_include_dir})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHPDDM_MPI=0")

# BLAS
find_package(BLAS REQUIRED)

# LAPACK
find_package(LAPACK REQUIRED)

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
add_executable(Tests_hpddm processing.cpp tests_hpddm.cpp)


target_link_libraries(Tests_hpddm ${BLAS_LIBRARIES})
target_link_libraries(Tests_hpddm ${LAPACK_LIBRARIES})
# set_property(TARGET Tests_hpddm PROPERTY CXX_STANDARD 11)
# set_property(TARGET Tests_hpddm PROPERTY CXX_STANDARD_REQUIRED ON)
