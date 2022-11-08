add_executable(2dtest_manual 2d_distribute.cpp)
target_link_libraries(2dtest_manual PUBLIC MPI::MPI_CXX lsb)

add_executable(2dtransmit_manual single_block_transmit/2d_trasmit.cpp)
target_link_libraries(2dtransmit_manual PUBLIC MPI::MPI_CXX lsb)

add_executable(3dtransmit_manual single_block_transmit/3d_trasmit.cpp)
target_link_libraries(3dtransmit_manual PUBLIC MPI::MPI_CXX lsb)

add_executable(4dtransmit_manual single_block_transmit/4d_trasmit.cpp)
target_link_libraries(4dtransmit_manual PUBLIC MPI::MPI_CXX lsb)

add_executable(5dtransmit_manual single_block_transmit/5d_trasmit.cpp)
target_link_libraries(5dtransmit_manual PUBLIC MPI::MPI_CXX lsb)