#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>

#define NI 4
#define NJ 6 
#define NK 8

#define NI_NEW 8
#define NJ_NEW 12
#define NK_NEW 2

#define SUB_NI 2
#define SUB_NJ 2
#define SUB_NK 2


#define CALIBRATE
#define NUM_RUNS 1
#define TIME_REQUIRED 0.5

int main(int argc, char** argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

    for (int i = 0; i < NI*NJ*NK; i++)
        originalArray[i] = 1;
    for (int i = 0; i < NI_NEW*NJ_NEW*NK_NEW; i++)
        newArray[i] = 0;
    MPI_Datatype send, recv;
    MPI_Request* sendreq = new MPI_Request[1];
    int subarray_size[3] = {SUB_NI,SUB_NJ,SUB_NK};

    int send_array_size[3] = {NI,NJ,NK};
    int send_start[3] = {0,0,0};
    MPI_Type_create_subarray(3, send_array_size, subarray_size, send_start, MPI_ORDER_C, MPI_INT, &send);
    MPI_Type_commit(&send);

    int recv_array_size[3] = {NI_NEW,NJ_NEW,NK_NEW};
    int recv_start[3] = {0,0,0};
    MPI_Type_create_subarray(3, recv_array_size, subarray_size, recv_start, MPI_ORDER_C, MPI_INT, &recv);
    MPI_Type_commit(&recv);

    double start, end;
    int num_runs;
    num_runs = NUM_RUNS;

#ifdef CALIBRATE
    while (num_runs < (1 << 14)) {
        start = MPI_Wtime();
        for (int k = 0; k < num_runs; ++k) {
            int count = 0;
            if (rank == 0)
            {
                MPI_Isend(&(originalArray[0]), 1, send, 1, 0, MPI_COMM_WORLD, &sendreq[0]);
            }

            if (rank == 1)
            {

                MPI_Recv(&(newArray[0]), 1, recv, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        end = MPI_Wtime();

        double ttime = (double)(end - start);

        if (ttime > TIME_REQUIRED) break;
        num_runs *= 2;
    }
#endif

    start = MPI_Wtime();
    for (int k = 0; k < num_runs; ++k) {
        int count = 0;
        if (rank == 0)
        {
            MPI_Isend(&(originalArray[0]), 1, send, 1, 0, MPI_COMM_WORLD, &sendreq[0]);
            MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
        }

        if (rank == 1)
        {

            MPI_Recv(&(newArray[0]), 1, recv, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    end = MPI_Wtime();

    std::cout << rank << " time : " << (end - start) / (double)num_runs * 1000000.0 << "us"<< std::endl;

    MPI_Finalize();

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
    delete[] originalArray; originalArray = nullptr;
    delete[] newArray; newArray = nullptr;
}
