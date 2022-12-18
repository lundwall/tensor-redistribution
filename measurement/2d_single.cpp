#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <omp.h>
#include <string>
#include <tuple>
#include <cassert>
#include <send_recv.hpp>
#include <validation.hpp>

int main(int argc, char** argv){
    constexpr size_t RUNS = 10;
    constexpr size_t N = 2;
    using T = int;
    int size, rank, received_threads;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &received_threads);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // size should be 2!
    if (size != 2) {
        throw std::runtime_error("When testing only 1 block transmission, only 2 processors are needed");
    }

    
    constexpr int NI = 20;
    constexpr int NJ = 30;
    constexpr int NI_NEW = 30;
    constexpr int NJ_NEW = 20;
    constexpr int SUB_NI = 3;
    constexpr int SUB_NJ = 6;
    constexpr int FROM_I = 5;
    constexpr int FROM_J = 10;
    constexpr int FROM_I_NEW = 5;
    constexpr int FROM_J_NEW = 10;
    constexpr int TO_I = FROM_I + SUB_NI;
    constexpr int TO_J = FROM_J + SUB_NJ;
    constexpr int TO_I_NEW = FROM_I_NEW + SUB_NI;
    constexpr int TO_J_NEW = FROM_J_NEW + SUB_NJ;

    constexpr int CHUNK_I = 1;
    constexpr int CHUNK_J = 1;
    constexpr int NUM_CHUNKS = 1;
    constexpr int CHUNK_SIZE = SUB_NI*SUB_NJ; // send all at once (same as no chunk)

    T* current_array;
    T* new_array;
    T* send_array;
    T* recv_array;

    NdIndices<N> chunk_num = {CHUNK_I, CHUNK_J};
    NdIndices<N> current_size = {NI, NJ};
    NdIndices<N> new_size = {NI_NEW, NJ_NEW};
    NdIndices<N> from = {FROM_I, FROM_J};
    NdIndices<N> to = {TO_I, TO_J};
    NdIndices<N> from_recv = {FROM_I_NEW, FROM_J_NEW};
    NdIndices<N> to_recv = {TO_I_NEW, TO_J_NEW};

    
    MPI_Datatype send_type, recv_type;
    int subarray_size[2] = {SUB_NI,SUB_NJ};
    int send_array_size[2] = {NI,NJ};
    int send_start[2] = {FROM_I,FROM_J};
    int recv_array_size[2] = {NI_NEW,NJ_NEW};
    int recv_start[2] = {FROM_I_NEW,FROM_J_NEW};
    MPI_Request* sendreq = new MPI_Request[1];
    MPI_Request* recvreq = new MPI_Request[1];

    char processor_name[256];
    int len_processor_name = 0;
    MPI_Get_processor_name(processor_name, &len_processor_name);
    std::cout << processor_name << std::endl;

    size_t current_total = get_product<N>(current_size);
    size_t new_total = get_product<N>(new_size);

    current_array = new T[current_total];
    new_array = new T[new_total];
    send_array = new T[SUB_NI*SUB_NJ];
    recv_array = new T[SUB_NI*SUB_NJ];

    // START METHOD 1
    std::string file_name = std::to_string(N) + std::string("d_transmit_with_API") + std::string(std::getenv("OMP_NUM_THREADS"));
    LSB_Init(file_name.c_str(), 0);
    LSB_Set_Rparam_int("rank", rank);
    set_lsb_chunk_size<N>(chunk_num);

    if (rank == 0) {
        for (int k = 0; k < RUNS; ++k) {
            LSB_Res();
            send<T, N>(current_array, 1, current_size, from, to, chunk_num);
            LSB_Rec(k);
        }
    }

    if(rank == 1){
        for (int k = 0; k < RUNS; ++k){
            LSB_Res();
            recv<T, N>(new_array, 0, new_size, from, to, chunk_num);
            LSB_Rec(k);
        }
    }
    LSB_Finalize();
    LSB_chunk_dim_cstr_free_all<N>();
    // END METHOD 1

    // START METHOD 2
    file_name = std::to_string(N) + std::string("d_transmit_custom_datatype") + std::string(std::getenv("OMP_NUM_THREADS"));
    LSB_Init(file_name.c_str(), 0);
    LSB_Set_Rparam_int("rank", rank);

    MPI_Type_create_subarray(2, send_array_size, subarray_size, send_start, MPI_ORDER_C, MPI_INT, &send_type);
    MPI_Type_commit(&send_type);

    MPI_Type_create_subarray(2, recv_array_size, subarray_size, recv_start, MPI_ORDER_C, MPI_INT, &recv_type);
    MPI_Type_commit(&recv_type);

    for (int k = 0; k < RUNS; ++k) {
        int count = 0;

        LSB_Res();

        if (rank == 0)
        {
            MPI_Isend(&(current_array[0]), 1, send_type, 1, 0, MPI_COMM_WORLD, &sendreq[0]);
            MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
        }

        if (rank == 1)
        {

            MPI_Recv(&(new_array[0]), 1, recv_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        LSB_Rec(k);
    }
    LSB_Finalize();

    // END METHOD 2


    MPI_Finalize();
    delete[] current_array; 
    current_array = nullptr;
    delete[] new_array;
    new_array = nullptr;
}
