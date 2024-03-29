#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <time.h>
#include "validation.h"

#define NI 192
#define NJ 240

#define NI_NEW 240
#define NJ_NEW 192

#define SUB_NI 54
#define SUB_NJ 45
#define RUNS 100

int main(int argc, char** argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    LSB_Init("2d_transmit_datatype", 0);
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

    int* originalArray = new int[NI*NJ];
    int* newArray = new int[NI_NEW*NJ_NEW];

    for (int i = 0; i < NI*NJ; i++)
        // originalArray[i] = 1;
        originalArray[i] = i;
    for (int i = 0; i < NI_NEW*NJ_NEW; i++)
        newArray[i] = 0;

    MPI_Datatype send, recv;
    MPI_Request* sendreq = new MPI_Request[1];
    int subarray_size[2] = {SUB_NI,SUB_NJ};

    int send_array_size[2] = {NI,NJ};
    int send_start[2] = {0,0};
    MPI_Type_create_subarray(2, send_array_size, subarray_size, send_start, MPI_ORDER_C, MPI_INT, &send);
    MPI_Type_commit(&send);

    int recv_array_size[2] = {NI_NEW,NJ_NEW};
    int recv_start[2] = {0,0};
    MPI_Type_create_subarray(2, recv_array_size, subarray_size, recv_start, MPI_ORDER_C, MPI_INT, &recv);
    MPI_Type_commit(&recv);

    srand(time(NULL));
    bool correctness = 1;

    for (int k = 0; k < RUNS; ++k) {
        int count = 0;
        // bool correctness = true;

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
        if (rank == 1){
            for (int i = 0; i < NI_NEW*NJ_NEW; i++)
                correctness &= newArray[i]==getElement(i, rank, NI, NJ, NI_NEW, NJ_NEW, SUB_NI, SUB_NJ);
        }
    }
    if (rank == 1){
        if(correctness) std::cout<<"Validation passed"<<std::endl;
        else std::cout<<"Validation failed"<<std::endl;
    }

    LSB_Finalize();
    MPI_Finalize();

    // for (int i = 0; i < NI_NEW; i++)
    // {
    //     for (int j = 0; j < NJ_NEW; j++)
    //         std::cout << newArray[i*NJ_NEW+j] << " ";
    //     std::cout << std::endl;
    // }
    delete[] originalArray; originalArray = nullptr;
    delete[] newArray; newArray = nullptr;
}
