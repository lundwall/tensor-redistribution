add_executable(2dtest 2d_distribute.cpp)
target_link_libraries(2dtest PUBLIC MPI::MPI_CXX)

add_executable(pgtest pg_distribute.cpp)
target_link_libraries(pgtest PUBLIC MPI::MPI_CXX)