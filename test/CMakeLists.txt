set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include") 
include_directories(${bemtool_dir})

# Intial test files
# add_executable(Prog prog.cxx)

# add_executable(Prog2 prog2.cxx)

# add_executable(Prog3 prog3.cxx)





# add_executable(Tests_circle processing.cpp tests_circle_2D.cpp)
# target_link_libraries(Tests_circle bemtool)

add_executable(Test plane_wave_harmonics_2D.cpp)










# BEMTOOL package
# set(bemtool_dir "${CMAKE_CURRENT_SOURCE_DIR}/../include") 
# include_directories(${bemtool_dir})

# Tests
# add_library(Test_circle tests_circle_2D.cpp tests_circle_2D.cpp)
# # Python embedding
# set(Python_ADDITIONAL_VERSIONS 3.4)
# find_package( PythonInterp 3 REQUIRED )
# find_package( PythonLibs 3 REQUIRED )
# include_directories ( ${PYTHON_INCLUDE_DIRS} )
# message("test")
# message("${Python_ADDITIONAL_VERSIONS}")
# Creation du test Test_loading_mesh a partir de loading_mesh
# add_executable(Test_loading_mesh loading_mesh.cpp )
# target_link_libraries (Test_loading_mesh  dunavant)

# add_executable(Test processing.cpp)
# target_link_libraries (Test  dunavant Test_circle)
# Executable pour first kind dirichlet 2D
# add_executable(Tests_circle_2D tests_circle_2D.cpp first_kind_dirichlet_2D.cpp tests_circle_2D.hpp)
# target_link_libraries (Tests_circle_2D  dunavant)

# Executable pour first kind neumann 2D
# add_executable(Test_first_kind_neumann_2D first_kind_neumann_2D.cpp )

# 
# # Executable pour second kind dirichlet 2D
# add_executable(Test_second_kind_dirichlet_2D second_kind_dirichlet_2D.cpp )
# target_link_libraries (Test_second_kind_dirichlet_2D  dunavant)
# 
# # Executable pour second kind neumann 2D
# add_executable(Test_second_kind_neumann_2D second_kind_neumann_2D.cpp )
# target_link_libraries (Test_second_kind_neumann_2D  dunavant)

# Executable pour les valeurs propres des operateurs integraux 2D
# add_executable(Test_fourier_harmonics_2D fourier_harmonics_2D.cpp )
# target_link_libraries (Test_fourier_harmonics_2D  dunavant)

# # Post-processing
# add_executable(Test_post_processing post_processing.cpp )
# target_link_libraries(Test_post_processing ${PYTHON_LIBRARIES})


# add_executable(Test_first_kind_dirichlet_3D first_kind_dirichlet_3D.cpp )
# target_link_libraries (Test_first_kind_dirichlet_3D  "${bemtool_dir}/quad2D.o")