add_executable(bandwidth "bandwidth.cpp")
target_link_libraries(bandwidth PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX lsb)
