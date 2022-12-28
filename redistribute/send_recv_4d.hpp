#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

void send_4d(int* source, int other_rank, int* current_size, int* from, int* to, MPI_Request *request, int* buffer) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int SUB_NK = to[2] - from[2];
    int SUB_NL = to[3] - from[3];
    int NI = current_size[0];
    int NJ = current_size[1];
    int NK = current_size[2];
    int NL = current_size[3];
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < SUB_NI; i++)
    {
        for (int j = 0; j < SUB_NJ; j++)
        {
            for (int k = 0; k < SUB_NK; k++)
            {
                memcpy(buffer + i*SUB_NJ*SUB_NK*SUB_NL + j*SUB_NK*SUB_NL + k*SUB_NL , source + (i+from[0])*NJ*NK*NL + (j+from[1])*NK*NL+ (k+from[2])*NL + from[3], sizeof(int)*SUB_NL);
            }
        }
    }
    MPI_Isend(buffer, SUB_NI*SUB_NJ*SUB_NK*SUB_NL, MPI_INT, other_rank, 0, MPI_COMM_WORLD, &request[0]);
}

void recv_4d(int* target, int other_rank, int* new_size, int* from, int* to) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int SUB_NK = to[2] - from[2];
    int SUB_NL = to[3] - from[3];
    int NI_NEW = new_size[0];
    int NJ_NEW = new_size[1];
    int NK_NEW = new_size[2];
    int NL_NEW = new_size[3];
    int* buffer = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL];
    MPI_Recv(&(buffer[0]), SUB_NI*SUB_NJ*SUB_NK*SUB_NL, MPI_INT, other_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < SUB_NI; i++)
    {
        for (int j = 0; j < SUB_NJ; j++)
        {
            for (int k = 0; k < SUB_NK; k++)
            {
                memcpy(target + (i+from[0])*NJ_NEW*NK_NEW*NL_NEW + (j+from[1])*NK_NEW*NL_NEW + (k+from[2])*NL_NEW + from[3], buffer + i*SUB_NJ*SUB_NK*SUB_NL + j*SUB_NK*SUB_NL + k*SUB_NL, sizeof(int)*SUB_NL);
            }
        }
    }
}