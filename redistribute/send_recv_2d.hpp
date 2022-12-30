#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

void send_2d(int* source, int other_rank, int* current_size, int* from, int* to, MPI_Request *request, int* buffer) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int NI = current_size[0];
    int NJ = current_size[1];
    #pragma omp parallel for collapse(4)
    for (int i = 0; i < SUB_NI; i++)
    {
        memcpy(buffer + i*SUB_NJ, source + (i+from[0])*NJ + from[1], sizeof(int)*SUB_NJ);
    }
    MPI_Isend(buffer, SUB_NI*SUB_NJ, MPI_INT, other_rank, 0, MPI_COMM_WORLD, &request[0]);
}

void recv_2d(int* target, int other_rank, int* new_size, int* from, int* to) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int NI_NEW = new_size[0];
    int NJ_NEW = new_size[1];
    int* buffer = new int[SUB_NI*SUB_NJ];
    MPI_Recv(&(buffer[0]), SUB_NI*SUB_NJ, MPI_INT, other_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp parallel for collapse(4)
    for (int i = 0; i < SUB_NI; i++)
    {
        memcpy(target + (i+from[0])*NJ_NEW + from[1], buffer + i*SUB_NJ, sizeof(int)*SUB_NJ);
    }
}