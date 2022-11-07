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
    LSB_Init("4d_transmit_datatype", 0);
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

    for (int i = 0; i < NI*NJ*NK*NL; i++)
        originalArray[i] = 1;
    for (int i = 0; i < NI_NEW*NJ_NEW*NK_NEW*NL_NEW; i++)
        newArray[i] = 0;
    MPI_Datatype send, recv;
    MPI_Request* sendreq = new MPI_Request[1];
    int subarray_size[4] = {SUB_NI,SUB_NJ,SUB_NK,SUB_NL};

    int send_array_size[4] = {NI,NJ,NK,NL};
    int send_start[4] = {0,0,0,0};
    MPI_Type_create_subarray(4, send_array_size, subarray_size, send_start, MPI_ORDER_C, MPI_INT, &send);
    MPI_Type_commit(&send);

    int recv_array_size[4] = {NI_NEW,NJ_NEW,NK_NEW,NL_NEW};
    int recv_start[4] = {0,0,0,0};
    MPI_Type_create_subarray(4, recv_array_size, subarray_size, recv_start, MPI_ORDER_C, MPI_INT, &recv);
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

    // for (int i = 0; i < NI_NEW*NJ_NEW*NK_NEW; i++)
    //     std::cout << newArray[i] << " ";
    delete[] originalArray; originalArray = nullptr;
    delete[] newArray; newArray = nullptr;
}
