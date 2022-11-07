add_executable(2dtest_manual 2d_distribute.cpp)
target_link_libraries(2dtest_manual PUBLIC MPI::MPI_CXX lsb)