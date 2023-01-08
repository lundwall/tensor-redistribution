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
#include "utils.hpp"
#include <algorithm>

int main(int argc, char** argv){
    size_t RUNS = 1000;
    size_t WARMUP = 10;
    size_t SYNC = 10;
    size_t N = 5;
    using T = int;
    int size, rank, received_threads;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &received_threads);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // size should be 2!
    if (size != 2) {
        throw std::runtime_error("When testing only 1 block transmission, only 2 processors are needed");
    }

    int experiment_size[6][5] = {
        {6,6,6,6,6},{3,3,3,3,24},{12,12,12,12,12},{6,6,6,6,48},{24,24,24,24,24},{12,12,12,12,96}
    };

    for (int exp_it = 0; exp_it < 6; exp_it++)
    {
	printf("start new experiment\n");
        int NI = 50;
        int NJ = 50;
        int NK = 50;
        int NL = 50;
        int NM = 100;

        int NI_NEW = 50;
        int NJ_NEW = 50;
        int NK_NEW = 50;
        int NL_NEW = 50;
        int NM_NEW = 100;

        int SUB_NI = experiment_size[exp_it][0];
        int SUB_NJ = experiment_size[exp_it][1];
        int SUB_NK = experiment_size[exp_it][2];
        int SUB_NL = experiment_size[exp_it][3];
        int SUB_NM = experiment_size[exp_it][4];

        int FROM_I = 0;
        int FROM_J = 0;
        int FROM_K = 0;
        int FROM_L = 0;
        int FROM_M = 0;

        int FROM_I_NEW = 0;
        int FROM_J_NEW = 0;
        int FROM_K_NEW = 0;
        int FROM_L_NEW = 0;
        int FROM_M_NEW = 0;

        int TO_I = FROM_I + SUB_NI;
        int TO_J = FROM_J + SUB_NJ;
        int TO_K = FROM_K + SUB_NK;
        int TO_L = FROM_L + SUB_NL;
        int TO_M = FROM_M + SUB_NM;

        int TO_I_NEW = FROM_I_NEW + SUB_NI;
        int TO_J_NEW = FROM_J_NEW + SUB_NJ;
        int TO_K_NEW = FROM_K_NEW + SUB_NK;
        int TO_L_NEW = FROM_L_NEW + SUB_NL;
        int TO_M_NEW = FROM_M_NEW + SUB_NM;

        T* current_array;
        T* new_array;
        T* send_array;
        T* recv_array;

        int from_int[N] = {FROM_I, FROM_J, FROM_K, FROM_L, FROM_M};
        int to_int[N] = {TO_I, TO_J, TO_K, TO_L, TO_M};
        int from_rec_int[N] = {FROM_I_NEW, FROM_J_NEW, FROM_K_NEW, FROM_L_NEW, FROM_M_NEW};
        int to_rec_int[N] = {TO_I_NEW, TO_J_NEW, TO_K_NEW, TO_L_NEW, TO_M_NEW};
        int current_size_int[N] = {NI, NJ, NK, NL, NM};
        int new_size_int[N] = {NI_NEW, NJ_NEW, NK_NEW, NL_NEW, NM_NEW};
        T* send_buffer = new T[(TO_I-FROM_I) * (TO_J-FROM_J) * (TO_K-FROM_K) * (TO_L-FROM_L) * (TO_M-FROM_M)];

        MPI_Datatype send_type, recv_type;
        int subarray_size[N] = {SUB_NI,SUB_NJ, SUB_NK, SUB_NL, SUB_NM};
        int send_array_size[N] = {NI,NJ,NK,NL,NM};
        int send_start[N] = {FROM_I,FROM_J,FROM_K,FROM_L,FROM_M};
        int recv_array_size[N] = {NI_NEW,NJ_NEW,NK_NEW,NL_NEW,NM_NEW};
        int recv_start[N] = {FROM_I_NEW,FROM_J_NEW,FROM_K_NEW,FROM_L_NEW,FROM_M_NEW};
        MPI_Request* sendreq = new MPI_Request[10];
        MPI_Request* recvreq = new MPI_Request[10];

        char processor_name[256];
        int len_processor_name = 0;
        MPI_Get_processor_name(processor_name, &len_processor_name);
        std::cout << processor_name << std::endl;

        size_t current_total = NI*NJ*NK*NL*NM;
        size_t new_total = NI_NEW*NJ_NEW*NK_NEW*NL_NEW*NM_NEW;

        current_array = new T[current_total];
        new_array = new T[new_total];
        send_array = new T[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];
        recv_array = new T[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];

        double win;
        double* recorded_values;
        bool all_finished;

        int thread_num;
        for (thread_num = 1; thread_num <= 8; thread_num++)
        {    
            omp_set_num_threads(thread_num);
	    MPI_Barrier(MPI_COMM_WORLD);
            // START METHOD 1
            std::string file_name = std::to_string(N) + std::string("d_transmit_custom_datatype_t") + std::to_string(omp_get_max_threads()) + std::string("_") + std::to_string(SUB_NI) + std::string("_") + std::to_string(SUB_NJ) + std::string("_") + std::to_string(SUB_NK) + std::string("_") + std::to_string(SUB_NL) + std::string("_") + std::to_string(SUB_NM);
            LSB_Init(file_name.c_str(), 0);
            LSB_Set_Rparam_int("rank", rank);
            LSB_Set_Rparam_string("mode", "datatype");
            LSB_Set_Rparam_int("threads", omp_get_max_threads());

            MPI_Type_create_subarray(N, send_array_size, subarray_size, send_start, MPI_ORDER_C, MPI_INT, &send_type);
            MPI_Type_commit(&send_type);

            MPI_Type_create_subarray(N, recv_array_size, subarray_size, recv_start, MPI_ORDER_C, MPI_INT, &recv_type);
            MPI_Type_commit(&recv_type);

            LSB_Set_Rparam_string("type", "sync");
            LSB_Set_Rparam_double("err", 0); // meaningless here
            LSB_Rec_disable();
            recorded_values = new double[1000];
            all_finished = false;
            for (int k = 0; k < WARMUP + SYNC + RUNS && !all_finished; ++k)
            {
                configure_LSB_and_sync(k, WARMUP, SYNC, &win);
                LSB_Res();

                if (rank == 0)
                {
                    MPI_Isend(&(current_array[0]), 1, send_type, 1, 0, MPI_COMM_WORLD, &sendreq[0]);
                    MPI_Recv(&(new_array[0]), 1, recv_type, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
                }

                if (rank == 1)
                {
                    MPI_Isend(&(current_array[0]), 1, send_type, 0, 0, MPI_COMM_WORLD, &sendreq[0]);
                    MPI_Recv(&(new_array[0]), 1, recv_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
                }

                int num_recorded_values = k - WARMUP - SYNC + 1;
                LSB_Rec(std::max(num_recorded_values, 0));
                if (num_recorded_values >= 1)
                {
                    aggregate_CIs(num_recorded_values, recorded_values, size, &all_finished);
                }
            }
            delete[] recorded_values;
            LSB_Finalize();

            // END METHOD 1

            // START METHOD 2
	    // int num_chunk = 1;
	    MPI_Barrier(MPI_COMM_WORLD);
            for (int num_chunk = 1; num_chunk < 6; num_chunk++)
            {
                file_name = std::to_string(N) + std::string("d_transmit_without_API_t") + std::to_string(omp_get_max_threads()) + std::string("_") + std::to_string(SUB_NI) + std::string("_") + std::to_string(SUB_NJ) + std::string("_") + std::to_string(SUB_NK) + std::string("_") + std::to_string(SUB_NL) + std::string("_") + std::to_string(SUB_NM);
                LSB_Init(file_name.c_str(), 0);
                LSB_Set_Rparam_int("rank", rank);
                LSB_Set_Rparam_string("mode", "manual");
                LSB_Set_Rparam_int("threads", omp_get_max_threads());
                LSB_Set_Rparam_int("chunk", num_chunk);

                LSB_Set_Rparam_string("type", "sync");
                LSB_Set_Rparam_double("err", 0); // meaningless here
                LSB_Rec_disable();
                recorded_values = new double[1000];
                all_finished = false;
	        MPI_Barrier(MPI_COMM_WORLD);
                for (int k = 0; k < WARMUP + SYNC + RUNS && !all_finished; ++k)
                {
                    configure_LSB_and_sync(k, WARMUP, SYNC, &win);
                    LSB_Res();
                    if (rank == 0) 
                    {
                        send_5d(current_array, 1, current_size_int, from_int, to_int, num_chunk, &sendreq[0], send_buffer);
                        recv_5d(new_array, 1, new_size_int, from_rec_int, to_rec_int, num_chunk);
                        MPI_Waitall(num_chunk, sendreq, MPI_STATUSES_IGNORE);
                    }

                    if(rank == 1)
                    {
                        send_5d(current_array, 0, current_size_int, from_int, to_int, num_chunk, &sendreq[0], send_buffer);
                        recv_5d(new_array, 0, new_size_int, from_rec_int, to_rec_int, num_chunk);
                        MPI_Waitall(num_chunk, sendreq, MPI_STATUSES_IGNORE);
                    }
                    int num_recorded_values = k - WARMUP - SYNC + 1;
                    LSB_Rec(std::max(num_recorded_values, 0));
                    if (num_recorded_values >= 1)
                    {
                        aggregate_CIs(num_recorded_values, recorded_values, size, &all_finished);
                    }
                }
                delete[] recorded_values;
                LSB_Finalize();
            }
            // END METHOD 2

            // START METHOD 3 manual with one-sided put
            file_name = std::to_string(N) + std::string("d_transmit_manual_put_t") + std::to_string(omp_get_max_threads()) + std::string("_") + std::to_string(SUB_NI) + std::string("_") + std::to_string(SUB_NJ) + std::string("_") + std::to_string(SUB_NK) + std::string("_") + std::to_string(SUB_NL) + std::string("_") + std::to_string(SUB_NM);
            LSB_Init(file_name.c_str(), 0);
            LSB_Set_Rparam_int("rank", rank);
            LSB_Set_Rparam_string("mode", "put_manual");
            LSB_Set_Rparam_int("threads", omp_get_max_threads());

            int transmit_size = SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM;
            MPI_Win window2;
            MPI_Win_create(recv_array,  transmit_size * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window2);

            LSB_Set_Rparam_string("type", "sync");
            LSB_Set_Rparam_double("err", 0); // meaningless here
            LSB_Rec_disable();
            recorded_values = new double[1000];
            all_finished = false;
	    MPI_Barrier(MPI_COMM_WORLD);
            for (int k = 0; k < WARMUP + SYNC + RUNS && !all_finished; ++k)
            {
                MPI_Win_fence(0, window2);
                configure_LSB_and_sync(k, WARMUP, SYNC, &win);
                LSB_Res();
                if (rank == 0)
                {
                    prepare_send_buffer(current_array, current_size_int, from_int, to_int, send_array);
                    MPI_Put(send_array, transmit_size, MPI_INT, 1, 0, transmit_size, MPI_INT, window2);
                    MPI_Win_fence(0, window2);
                    unpack_recv_buffer(new_array, new_size_int, from_rec_int, to_rec_int, recv_array);
                }
                if(rank == 1)
                {
                    prepare_send_buffer(current_array, current_size_int, from_int, to_int, send_array);
                    MPI_Put(send_array, transmit_size, MPI_INT, 0, 0, transmit_size, MPI_INT, window2);
                    MPI_Win_fence(0, window2);
                    unpack_recv_buffer(new_array, new_size_int, from_rec_int, to_rec_int, recv_array);
                }
                int num_recorded_values = k - WARMUP - SYNC + 1;
                LSB_Rec(std::max(num_recorded_values, 0));
                if (num_recorded_values >= 1)
                {
                    aggregate_CIs(num_recorded_values, recorded_values, size, &all_finished);
                }
            }
            MPI_Win_free(&window2);
            delete[] recorded_values;
            LSB_Finalize();
            // END METHOD 3

        }
        delete[] current_array; 
        current_array = nullptr;
        delete[] new_array;
        new_array = nullptr;
        delete[] send_buffer;
        send_buffer = nullptr;
        delete[] send_array;
        send_array = nullptr;
        delete[] recv_array;
        recv_array = nullptr;
        delete[] sendreq;
        sendreq = nullptr;
        delete[] recvreq;
        recvreq = nullptr;
    }
    MPI_Finalize();
}
