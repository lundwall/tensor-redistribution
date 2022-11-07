#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <time.h>

#define NI 4
#define NJ 6 
#define NK 8

#define NI_NEW 6
#define NJ_NEW 4
#define NK_NEW 8

#define SUB_NI 4
#define SUB_NJ 2
#define SUB_NK 3

#define RUNS 10

int main(int argc, char** argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    LSB_Init("3d_transmit", 0);
    LSB_Set_Rparam_int("rank", rank);
    LSB_Set_Rparam_int("P", size);

    // size should be 2!
    if (size != 2)
    {
        printf("When testing only 1 block transmission, only 2 processors are needed\n");
        return 0;
    }

    char* processor_name = new char[256]; int len_processor_name = 0;
    MPI_Get_processor_name(processor_name, &len_processor_name);
    std::cout << processor_name << std::endl;

    int* originalArray = new int[NI*NJ*NK];
    int* newArray = new int[NI_NEW*NJ_NEW*NK_NEW];
    int* sendArray = new int[SUB_NI*SUB_NJ*SUB_NK];
    int* recvArray = new int[SUB_NI*SUB_NJ*SUB_NK];

    for (int i = 0; i < NI*NJ*NK; i++)
        originalArray[i] = 1;
    for (int i = 0; i < NI_NEW*NJ_NEW*NK_NEW; i++)
        newArray[i] = 0;
    MPI_Datatype send, recv;
    MPI_Request* sendreq = new MPI_Request[1];

    srand(time(NULL));

    for (int k = 0; k < RUNS; ++k) {
        int count = 0;

        LSB_Res();

        if (rank == 0)
        {
            for (int i = 0; i < SUB_NI; i++)
            {
                for (int j = 0; j < SUB_NJ; j++)
                {
                    memcpy(sendArray+i*SUB_NJ*SUB_NK+j*SUB_NK, originalArray+i*NJ*NK+j*NK, sizeof(int)*SUB_NK);
                }
            }
            MPI_Isend(&(sendArray[0]), SUB_NI*SUB_NJ*SUB_NK, MPI_INT, 1, 0, MPI_COMM_WORLD, &sendreq[0]);
            MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
        }

        if (rank == 1)
        {
            MPI_Recv(&(recvArray[0]), SUB_NI*SUB_NJ*SUB_NK, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < SUB_NI; i++)
            {
                for (int j = 0; j < SUB_NJ; j++)
                {
                    memcpy(newArray+i*NJ_NEW*NK_NEW+j*NK_NEW, recvArray+i*SUB_NJ*SUB_NK+j*SUB_NK, sizeof(int)*SUB_NK);
                }
            }
        }

        LSB_Rec(k);
    }

    LSB_Finalize();
    MPI_Finalize();

    if (rank == 1)
    {
        for (int i = 0; i < NI_NEW; i++)
        {
            for (int j = 0; j < NJ_NEW; j++)
            {
                for (int k = 0; k < NK_NEW; k++)
                    std::cout << newArray[i*NJ_NEW*NK_NEW+j*NK_NEW+k] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl << std::endl;
        }
    }
    delete[] originalArray; originalArray = nullptr;
    delete[] newArray; newArray = nullptr;
}