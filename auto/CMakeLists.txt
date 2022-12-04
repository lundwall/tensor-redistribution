function(CREATE_TEST_AUTO NUM_D)
    add_executable(${NUM_D}d_test_auto "single_block_transmit/${NUM_D}d_transmit.cpp")
    target_include_directories(${NUM_D}d_test_auto PUBLIC "include")
    target_link_libraries(${NUM_D}d_test_auto PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX lsb)
endfunction()

CREATE_TEST_AUTO(2)
CREATE_TEST_AUTO(3)