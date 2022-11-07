add_executable(2dtest 2d_distribute.cpp)
target_link_libraries(2dtest PUBLIC MPI::MPI_CXX)

add_executable(pgtest pg_distribute.cpp)
target_link_libraries(pgtest PUBLIC MPI::MPI_CXX)

add_executable(2dtransmit single_block_transmit/2d_trasmit.cpp)
target_link_libraries(2dtransmit PUBLIC MPI::MPI_CXX lsb)

add_executable(3dtransmit single_block_transmit/3d_trasmit.cpp)
target_link_libraries(3dtransmit PUBLIC MPI::MPI_CXX lsb)

add_executable(4dtransmit single_block_transmit/4d_trasmit.cpp)
target_link_libraries(4dtransmit PUBLIC MPI::MPI_CXX lsb)