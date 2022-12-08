add_executable(redistribute "redistribute.cpp")
target_include_directories(redistribute PUBLIC "../auto/include")
target_link_libraries(redistribute PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX lsb)