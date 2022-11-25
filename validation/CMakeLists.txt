add_executable(validation_demo 2d_trasmit.cpp)
target_link_libraries(validation_demo PUBLIC MPI::MPI_CXX lsb)