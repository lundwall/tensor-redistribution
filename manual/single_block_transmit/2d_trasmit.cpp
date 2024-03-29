#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <omp.h>
#include <string>

#define NI 19200
#define NJ 24000

#define NI_NEW 24000
#define NJ_NEW 19200

#define SUB_NI 5400
#define SUB_NJ 4500

#define RUNS 100
// #define CHUNK_SIZE 4050000 // 900*4500, 1/6 of all the data
// #define NUM_CHUNKS SUB_NI*SUB_NJ/CHUNK_SIZE

int main(int argc, char** argv)
{
    int CHUNK_SIZE = atoi(argv[1]);
    int NUM_CHUNKS = SUB_NI*SUB_NJ/CHUNK_SIZE;

    int size, rank, received_threads;
    std::string name_string = "2d_transmit_manual"+std::string(std::getenv("OMP_NUM_THREADS"));
    const char* liblsb_fname = name_string.c_str();

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &received_threads);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    LSB_Init(liblsb_fname, 0);
    LSB_Set_Rparam_int("rank", rank);
    int max_threads = omp_get_max_threads();
    // LSB_Set_Rparam_int("threads", max_threads);
    LSB_Set_Rparam_int("chunk size", CHUNK_SIZE);

    // size should be 2!
    if (size != 2)
    {
        printf("When testing only 1 block transmission, only 2 processors are needed\n");
        return 0;
    }

    char* processor_name = new char[256]; int len_processor_name = 0;
    MPI_Get_processor_name(processor_name, &len_processor_name);
    std::cout << processor_name << std::endl;

    int* originalArray = new int[NI*NJ];
    int* newArray = new int[NI_NEW*NJ_NEW];
    int* sendArray = new int[SUB_NI*SUB_NJ];
    int* recvArray = new int[SUB_NI*SUB_NJ];
    MPI_Request* sendreq = new MPI_Request[1];
    MPI_Request* recvreq = new MPI_Request[1];

    for (int i = 0; i < NI*NJ; i++)
        originalArray[i] = 1;
    for (int i = 0; i < NI_NEW*NJ_NEW; i++)
        newArray[i] = 0;

    for (int k = 0; k < RUNS; ++k)
    {
        int count = 0;

        if (rank == 0)
        {
            LSB_Res();
            for (int chunk = 0; chunk < NUM_CHUNKS; chunk++)
            {
                #pragma omp parallel
                #pragma omp single
                for (int i = 0; i < CHUNK_SIZE / SUB_NJ; i++)
                {
                    #pragma omp task
                    memcpy(sendArray + chunk * CHUNK_SIZE + i * SUB_NJ, originalArray + (chunk * CHUNK_SIZE / SUB_NJ) * NJ + i * NJ, sizeof(int) * SUB_NJ);
                }
                if (chunk > 0) {
                    MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
                }
                MPI_Isend(&(sendArray[chunk*CHUNK_SIZE]), CHUNK_SIZE, MPI_INT, 1, chunk, MPI_COMM_WORLD, sendreq);
            }
            MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
            LSB_Rec(k);
        }

        if (rank == 1)
        {
            LSB_Res();
            MPI_Irecv(&(recvArray[0]), CHUNK_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD, recvreq);
            for (int chunk = 0; chunk < NUM_CHUNKS; chunk++)
            {
                MPI_Waitall(1, recvreq, MPI_STATUSES_IGNORE);
                if (chunk != NUM_CHUNKS - 1)
                {
                    MPI_Irecv(&(recvArray[(chunk+1)*CHUNK_SIZE]), CHUNK_SIZE, MPI_INT, 0, chunk+1, MPI_COMM_WORLD, recvreq);
                }
                #pragma omp parallel
                #pragma omp single
                for (int i = 0; i < CHUNK_SIZE / SUB_NJ; i++)
                {
                    #pragma omp task
                    memcpy(newArray + (chunk * CHUNK_SIZE / SUB_NJ) * NJ_NEW + i * NJ_NEW, recvArray + chunk * CHUNK_SIZE + i * SUB_NJ, sizeof(int) * SUB_NJ);
                }
            }
            LSB_Rec(k);
        }
    }

    LSB_Finalize();
    MPI_Finalize();

    // if (rank == 1)
    // {
    //     for (int i = 0; i < NI_NEW; i++)
    //     {
    //         for (int j = 0; j < NJ_NEW; j++)
    //             std::cout << newArray[i*NJ_NEW+j] << " ";
    //         std::cout << std::endl;
    //     }
    // }
    delete[] originalArray; originalArray = nullptr;
    delete[] newArray; newArray = nullptr;
}
