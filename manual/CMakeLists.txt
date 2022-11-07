add_executable(2dtest_manual 2d_distribute.cpp)
target_link_libraries(2dtest_manual PUBLIC MPI::MPI_CXX lsb)
add_executable(2dtransmit_manual single_block_transmit/2d_trasmit.cpp)
target_link_libraries(2dtransmit_manual PUBLIC MPI::MPI_CXX lsb)