add_executable(2dtest 2d_distribute.cpp)
target_link_libraries(2dtest PUBLIC MPI::MPI_CXX)