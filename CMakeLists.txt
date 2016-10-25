cmake_minimum_required(VERSION 2.6)
Project(Bemtool)

# Cmake files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")

# Source files
add_subdirectory(src)


# # List the headers and add them into the project
# FILE(GLOB_RECURSE MY_HEADERS "include/*.hpp")
# include_directories(include)

# Add documentation
add_subdirectory(doc)

# Add test files
add_subdirectory(test)

# Testing project
# enable_testing()
# add_test(NAME loading_mesh_test COMMAND Test_loading_mesh)

# add_test(NAME tests_circle_2D COMMAND Tests_circle_2D)
# add_test(NAME first_kind_neumann_2D_test COMMAND Test_first_kind_neumann_2D)
# add_test(NAME second_kind_dirichlet_2D_test COMMAND Test_second_kind_dirichlet_2D)
# add_test(NAME second_kind_neumann_2D_test COMMAND Test_second_kind_neumann_2D)

# add_test(NAME fourier_harmonics_2D_test COMMAND Test_fourier_harmonics_2D)

# add_test(NAME post_processing_test COMMAND Test_post_processing)
# add_test(NAME first_kind_dirichlet_3D_test COMMAND Test_first_kind_dirichlet_3D)

