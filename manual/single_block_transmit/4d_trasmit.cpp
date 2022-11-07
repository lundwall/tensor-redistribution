#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <time.h>

#define NI 20
#define NJ 40
#define NK 60
#define NL 80

#define NI_NEW 40
#define NJ_NEW 20
#define NK_NEW 80
#define NL_NEW 60

#define SUB_NI 20
#define SUB_NJ 20
#define SUB_NK 20
#define SUB_NL 20

#define RUNS 100

int main(int argc, char** argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    LSB_Init("4d_transmit_manual", 0);
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

    int* originalArray = new int[NI*NJ*NK*NL];
    int* newArray = new int[NI_NEW*NJ_NEW*NK_NEW*NL_NEW];
    int* sendArray = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL];
    int* recvArray = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL];
    MPI_Request* sendreq = new MPI_Request[1];

    for (int i = 0; i < NI*NJ*NK*NL; i++)
        originalArray[i] = 1;
    for (int i = 0; i < NI_NEW*NJ_NEW*NK_NEW*NL_NEW; i++)
        newArray[i] = 0;

    srand(time(NULL));

    for (int r = 0; r < RUNS; ++r) {
        int count = 0;

        LSB_Res();

        if (rank == 0)
        {
            for (int i = 0; i < SUB_NI; i++)
            {
                for (int j = 0; j < SUB_NJ; j++)
                {
                    for (int k = 0; k < SUB_NK; k++)
                    {
                        memcpy(sendArray+ i*SUB_NJ*SUB_NK*SUB_NL + j*SUB_NK*SUB_NL + k*SUB_NL, originalArray + i*NJ*NK*NL + j*NK*NL + k*NL, sizeof(int)*SUB_NL);
                    }
                }
            }
            MPI_Isend(&(sendArray[0]), SUB_NI*SUB_NJ*SUB_NK*SUB_NL, MPI_INT, 1, 0, MPI_COMM_WORLD, &sendreq[0]);
            MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
        }

        if (rank == 1)
        {
            MPI_Recv(&(recvArray[0]), SUB_NI*SUB_NJ*SUB_NK*SUB_NL, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < SUB_NI; i++)
            {
                for (int j = 0; j < SUB_NJ; j++)
                {
                    for (int k = 0; k < SUB_NK; k++)
                    {
                        memcpy(newArray+i*NJ_NEW*NK_NEW*NL_NEW+j*NK_NEW*NL_NEW+k*NL_NEW, recvArray+i*SUB_NJ*SUB_NK*SUB_NL+j*SUB_NK*SUB_NL, sizeof(int)*SUB_NL);
                    }
                }
            }
        }

        LSB_Rec(r);
    }

    LSB_Finalize();
    MPI_Finalize();

    // if (rank == 1)
    // {
    //     for (int i = 0; i < NI_NEW; i++)
    //     {
    //         for (int j = 0; j < NJ_NEW; j++)
    //         {
    //             for (int k = 0; k < NK_NEW; k++)
    //             {
    //                 for (int l = 0; l < NL_NEW; l++) {

    //                     std::cout << newArray[i*NJ_NEW*NK_NEW*NL_NEW+j*NK_NEW*NL_NEW+k*NL_NEW+l] << " ";
    //                 }
    //                 std::cout << std::endl;
    //             }
    //             std::cout << std::endl << std::endl;
    //         }
    //         std::cout << std::endl << std::endl << std::endl;
    //     }
    // }

    delete[] originalArray; originalArray = nullptr;
    delete[] newArray; newArray = nullptr;
}
