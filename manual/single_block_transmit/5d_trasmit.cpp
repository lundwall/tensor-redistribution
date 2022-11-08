#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <time.h>

#define NI 2
#define NJ 4
#define NK 6
#define NL 8
#define NM 6

#define NI_NEW 4
#define NJ_NEW 2
#define NK_NEW 8
#define NL_NEW 6
#define NM_NEW 6

#define SUB_NI 2
#define SUB_NJ 2
#define SUB_NK 2
#define SUB_NL 2
#define SUB_NM 2

#define RUNS 100

int main(int argc, char** argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    LSB_Init("5d_transmit_manual", 0);
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

    int* originalArray = new int[NI*NJ*NK*NL*NM];
    int* newArray = new int[NI_NEW*NJ_NEW*NK_NEW*NL_NEW*NM_NEW];
    int* sendArray = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];
    int* recvArray = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];
    MPI_Request* sendreq = new MPI_Request[1];

    for (int i = 0; i < NI*NJ*NK*NL*NM; i++)
        originalArray[i] = 1;
    for (int i = 0; i < NI_NEW*NJ_NEW*NK_NEW*NL_NEW*NM_NEW; i++)
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
                        for (int l = 0; l < SUB_NL; l++)
                        {
                            memcpy(sendArray+ i*SUB_NJ*SUB_NK*SUB_NL*SUB_NM + j*SUB_NK*SUB_NL*SUB_NM + k*SUB_NL*SUB_NM + l*SUB_NM, originalArray + i*NJ*NK*NL*NM + j*NK*NL*NM + k*NL*NM + l*NM, sizeof(int)*SUB_NM);
                        }
                    }
                }
            }
            MPI_Isend(&(sendArray[0]), SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM, MPI_INT, 1, 0, MPI_COMM_WORLD, &sendreq[0]);
            MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
        }

        if (rank == 1)
        {
            MPI_Recv(&(recvArray[0]), SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < SUB_NI; i++)
            {
                for (int j = 0; j < SUB_NJ; j++)
                {
                    for (int k = 0; k < SUB_NK; k++)
                    {
                        for (int l = 0; l < SUB_NL; l++)
                        {
                            memcpy(newArray+i*NJ_NEW*NK_NEW*NL_NEW*NM_NEW+j*NK_NEW*NL_NEW*NM_NEW+k*NL_NEW*NM_NEW+l*NM_NEW, recvArray+i*SUB_NJ*SUB_NK*SUB_NL*SUB_NM+j*SUB_NK*SUB_NL*SUB_NM+k*SUB_NL*SUB_NM+l*SUB_NM, sizeof(int)*SUB_NM);
                        }
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
    //                 for (int l = 0; l < NL_NEW; l++) 
    //                 {
    //                     for (int m = 0; m < NM_NEW; m++)
    //                     {
    //                         std::cout << newArray[i*NJ_NEW*NK_NEW*NL_NEW*NM_NEW+j*NK_NEW*NL_NEW*NM_NEW+k*NL_NEW*NM_NEW+l*NM_NEW+m] << " ";
    //                     }
    //                     std::cout << std::endl;
    //                 }
    //                 std::cout << std::endl << std::endl;
    //             }
    //             std::cout << std::endl << std::endl << std::endl;
    //         }
    //         std::cout << std::endl << std::endl << std::endl << std::endl;
    //     }
    // }

    delete[] originalArray; originalArray = nullptr;
    delete[] newArray; newArray = nullptr;
}
