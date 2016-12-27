# BEMTOOL package
set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../../include")
include_directories(${bemtool_dir})

# HPDDM
set(hpddm_include_dir /Volumes/Data/Repositories/hpddm/include)
include_directories(${hpddm_include_dir})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHPDDM_MPI=0")

# BLAS
find_package(BLAS REQUIRED)

# LAPACK
find_package(LAPACK REQUIRED)

#=============================================================================#
#=========================== Executables =====================================#
#=============================================================================#
# add_executable(Tests_hpddm tests_hpddm.cpp
#                     "${bemtool_dir}")
#
#
# target_link_libraries(Tests_circle_2D ${BLAS_LIBRARIES})
# target_link_libraries(Tests_circle_2D ${LAPACK_LIBRARIES})
