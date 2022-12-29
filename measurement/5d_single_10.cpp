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
#include <send_recv_5d.hpp>
#include "oneside_helper.h"

int main(int argc, char** argv){
    constexpr size_t RUNS = 100;
    constexpr size_t N = 5;
    using T = int;
    int size, rank, received_threads;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &received_threads);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // size should be 2!
    if (size != 2) {
        throw std::runtime_error("When testing only 1 block transmission, only 2 processors are needed");
    }

    int thread_num;
    for (thread_num = 1; thread_num <= 8; thread_num++)
    {    
        omp_set_num_threads(thread_num);
        constexpr int NI = 70;
        constexpr int NJ = 70;
        constexpr int NK = 70;
        constexpr int NL = 70;
        constexpr int NM = 70;

        constexpr int NI_NEW = 70;
        constexpr int NJ_NEW = 70;
        constexpr int NK_NEW = 70;
        constexpr int NL_NEW = 70;
        constexpr int NM_NEW = 70;

        constexpr int SUB_NI = 10;
        constexpr int SUB_NJ = 10;
        constexpr int SUB_NK = 10;
        constexpr int SUB_NL = 10;
        constexpr int SUB_NM = 10;

        constexpr int FROM_I = 5;
        constexpr int FROM_J = 5;
        constexpr int FROM_K = 5;
        constexpr int FROM_L = 5;
        constexpr int FROM_M = 5;

        constexpr int FROM_I_NEW = 5;
        constexpr int FROM_J_NEW = 5;
        constexpr int FROM_K_NEW = 5;
        constexpr int FROM_L_NEW = 5;
        constexpr int FROM_M_NEW = 5;

        constexpr int TO_I = FROM_I + SUB_NI;
        constexpr int TO_J = FROM_J + SUB_NJ;
        constexpr int TO_K = FROM_K + SUB_NK;
        constexpr int TO_L = FROM_L + SUB_NL;
        constexpr int TO_M = FROM_M + SUB_NM;

        constexpr int TO_I_NEW = FROM_I_NEW + SUB_NI;
        constexpr int TO_J_NEW = FROM_J_NEW + SUB_NJ;
        constexpr int TO_K_NEW = FROM_K_NEW + SUB_NK;
        constexpr int TO_L_NEW = FROM_L_NEW + SUB_NL;
        constexpr int TO_M_NEW = FROM_M_NEW + SUB_NM;

        constexpr int CHUNK_I = 1;
        constexpr int CHUNK_J = 1;
        constexpr int CHUNK_K = 1;
        constexpr int CHUNK_L = 1;
        constexpr int CHUNK_M = 1;
        constexpr int NUM_CHUNKS = 1;
        constexpr int CHUNK_SIZE = SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM; // send all at once (same as no chunk)

        T* current_array;
        T* new_array;
        T* send_array;
        T* recv_array;

        NdIndices<N> chunk_num = {CHUNK_I, CHUNK_J, CHUNK_J, CHUNK_L, CHUNK_M};
        NdIndices<N> current_size = {NI, NJ, NK, NL, NM};
        NdIndices<N> new_size = {NI_NEW, NJ_NEW, NK_NEW, NL_NEW, NM_NEW};
        NdIndices<N> from = {FROM_I, FROM_J, FROM_K, FROM_L, FROM_M};
        NdIndices<N> to = {TO_I, TO_J, TO_K, TO_L, TO_M};
        NdIndices<N> from_recv = {FROM_I_NEW, FROM_J_NEW, FROM_K_NEW, FROM_L_NEW, FROM_M_NEW};
        NdIndices<N> to_recv = {TO_I_NEW, TO_J_NEW, TO_K_NEW, TO_L_NEW, TO_M_NEW};

        NdIndices<N> range = to - from;
        size_t sending_total = get_product<N>(range);
        T* send_buffer = new T[sending_total];

        int from_int[N] = {FROM_I, FROM_J, FROM_K, FROM_L, FROM_M};
        int to_int[N] = {TO_I, TO_J, TO_K, TO_L, TO_M};
        int from_rec_int[N] = {FROM_I_NEW, FROM_J_NEW, FROM_K_NEW, FROM_L_NEW, FROM_M_NEW};
        int to_rec_int[N] = {TO_I_NEW, TO_J_NEW, TO_K_NEW, TO_L_NEW, TO_M_NEW};
        int current_size_int[N] = {NI, NJ, NK, NL, NM};
        int new_size_int[N] = {NI_NEW, NJ_NEW, NK_NEW, NL_NEW, NM_NEW};
        
        MPI_Datatype send_type, recv_type;
        int subarray_size[N] = {SUB_NI,SUB_NJ, SUB_NK, SUB_NL, SUB_NM};
        int send_array_size[N] = {NI,NJ,NK,NL,NM};
        int send_start[N] = {FROM_I,FROM_J,FROM_K,FROM_L,FROM_M};
        int recv_array_size[N] = {NI_NEW,NJ_NEW,NK_NEW,NL_NEW,NM_NEW};
        int recv_start[N] = {FROM_I_NEW,FROM_J_NEW,FROM_K_NEW,FROM_L_NEW,FROM_M_NEW};
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
        send_array = new T[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];
        recv_array = new T[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];

        // MPI_Barrier(MPI_COMM_WORLD);

        // START METHOD 1
        std::string file_name = std::to_string(N) + std::string("d_transmit_with_API_t") + std::to_string(omp_get_max_threads()) + std::string("_") + std::to_string(SUB_NI);
        LSB_Init(file_name.c_str(), 0);
        LSB_Set_Rparam_int("rank", rank);
        set_lsb_chunk_size<N>(chunk_num);
        for (int k = 0; k < RUNS; ++k)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            LSB_Res();
            if (rank == 0) 
            {
                send<T, N>(current_array, 1, current_size, from, to, chunk_num, &sendreq[0], send_buffer);
                MPI_Waitall(1, &sendreq[0], MPI_STATUS_IGNORE);
            }

            if(rank == 1)
            {
                recv<T, N>(new_array, 0, new_size, from, to, chunk_num);
            }
            LSB_Rec(k);
        }
        LSB_Finalize();
        LSB_chunk_dim_cstr_free_all<N>();
        // END METHOD 1

        // MPI_Barrier(MPI_COMM_WORLD);

        // START METHOD 2
        file_name = std::to_string(N) + std::string("d_transmit_custom_datatype_t") + std::to_string(omp_get_max_threads()) + std::string("_") + std::to_string(SUB_NI);
        LSB_Init(file_name.c_str(), 0);    LSB_Set_Rparam_int("rank", rank);

        MPI_Type_create_subarray(N, send_array_size, subarray_size, send_start, MPI_ORDER_C, MPI_INT, &send_type);
        MPI_Type_commit(&send_type);

        MPI_Type_create_subarray(N, recv_array_size, subarray_size, recv_start, MPI_ORDER_C, MPI_INT, &recv_type);
        MPI_Type_commit(&recv_type);

        for (int k = 0; k < RUNS; ++k) 
        {
            MPI_Barrier(MPI_COMM_WORLD);
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

        MPI_Barrier(MPI_COMM_WORLD);

        // START METHOD 3
        file_name = std::to_string(N) + std::string("d_transmit_without_API_t") + std::to_string(omp_get_max_threads()) + std::string("_") + std::to_string(SUB_NI);
        LSB_Init(file_name.c_str(), 0);
        LSB_Set_Rparam_int("rank", rank);

        for (int k = 0; k < RUNS; ++k)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            LSB_Res();
            if (rank == 0) 
            {
                send_5d(current_array, 1, current_size_int, from_int, to_int, &sendreq[0], send_buffer);
                MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
            }

            if(rank == 1)
            {
                recv_5d(new_array, 0, new_size_int, from_rec_int, to_rec_int);
            }
            LSB_Rec(k);
        }
        LSB_Finalize();
        // END METHOD 3

        // MPI_Barrier(MPI_COMM_WORLD);

        // START METHOD 4 datatype with one-sided put
        file_name = std::to_string(N) + std::string("d_transmit_custom_datatype_put_t") + std::to_string(omp_get_max_threads()) + std::string("_") + std::to_string(SUB_NI);
        LSB_Init(file_name.c_str(), 0);
        LSB_Set_Rparam_int("rank", rank);

        MPI_Win window1;
        MPI_Win_create(new_array, new_total * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window1);
        MPI_Win_fence(0, window1);

        for (int k = 0; k < RUNS; ++k) 
        {
            MPI_Barrier(MPI_COMM_WORLD);
            LSB_Res();
            if (rank == 0)
            {
                MPI_Put(current_array, 1, send_type, 1, 0, 1, recv_type, window1);
            }
            MPI_Win_fence(0, window1);
            LSB_Rec(k);
        }
        MPI_Win_free(&window1);
        LSB_Finalize();
        // END METHOD 4

        // MPI_Barrier(MPI_COMM_WORLD);

        // START METHOD 5 manual with one-sided put
        file_name = std::to_string(N) + std::string("d_transmit_manual_put_t") + std::to_string(omp_get_max_threads()) + std::string("_") + std::to_string(SUB_NI);
        LSB_Init(file_name.c_str(), 0);
        LSB_Set_Rparam_int("rank", rank);

        int transmit_size = SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM;
        MPI_Win window2;
        MPI_Win_create(recv_array,  transmit_size * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window2);
        MPI_Win_fence(0, window2);

        for (int k = 0; k < RUNS; ++k) {
            MPI_Barrier(MPI_COMM_WORLD);
            LSB_Res();
            if (rank == 0)
            {
                prepare_send_buffer(current_array, current_size_int, from_int, to_int, send_array);
                MPI_Put(send_array, transmit_size, MPI_INT, 1, 0, transmit_size, MPI_INT, window2);
            }
            MPI_Win_fence(0, window2);
            if(rank == 1)
            {
                unpack_recv_buffer(new_array, new_size_int, from_rec_int, to_rec_int, recv_array);
            }
            LSB_Rec(k);
        }
        MPI_Win_free(&window2);
        LSB_Finalize();

        MPI_Finalize();
        delete[] current_array; 
        current_array = nullptr;
        delete[] new_array;
        new_array = nullptr;
        delete[] send_buffer;
        send_buffer = nullptr;
    }
}
