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
    LSB_Init("5d_transmit_datatype", 0);
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

    for (int i = 0; i < NI*NJ*NK*NL*NM; i++)
        originalArray[i] = 1;
    for (int i = 0; i < NI_NEW*NJ_NEW*NK_NEW*NL_NEW*NM_NEW; i++)
        newArray[i] = 0;
    MPI_Datatype send, recv;
    MPI_Request* sendreq = new MPI_Request[1];
    int subarray_size[5] = {SUB_NI,SUB_NJ,SUB_NK,SUB_NL,SUB_NM};

    int send_array_size[5] = {NI,NJ,NK,NL,NM};
    int send_start[5] = {0,0,0,0,0};
    MPI_Type_create_subarray(5, send_array_size, subarray_size, send_start, MPI_ORDER_C, MPI_INT, &send);
    MPI_Type_commit(&send);

    int recv_array_size[5] = {NI_NEW,NJ_NEW,NK_NEW,NL_NEW,NM_NEW};
    int recv_start[5] = {0,0,0,0,0};
    MPI_Type_create_subarray(5, recv_array_size, subarray_size, recv_start, MPI_ORDER_C, MPI_INT, &recv);
    MPI_Type_commit(&recv);

    srand(time(NULL));

    for (int k = 0; k < RUNS; ++k) {
        int count = 0;

        LSB_Res();

        if (rank == 0)
        {
            MPI_Isend(&(originalArray[0]), 1, send, 1, 0, MPI_COMM_WORLD, &sendreq[0]);
            MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
        }

        if (rank == 1)
        {

            MPI_Recv(&(newArray[0]), 1, recv, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        LSB_Rec(k);
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
