#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

void send_3d(int* source, int other_rank, int* current_size, int* from, int* to, MPI_Request *request, int* buffer) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int SUB_NK = to[2] - from[2];
    int NI = current_size[0];
    int NJ = current_size[1];
    int NK = current_size[2];
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < SUB_NI; i++)
    {
        for (int j = 0; j < SUB_NJ; j++)
        {
            memcpy(buffer + i*SUB_NJ*SUB_NK + j*SUB_NK, source + (i+from[0])*NJ*NK + (j+from[1])*NK + (from[2]), sizeof(int)*SUB_NK);
        }
    }
    MPI_Isend(buffer, SUB_NI*SUB_NJ*SUB_NK, MPI_INT, other_rank, 0, MPI_COMM_WORLD, &request[0]);
}

void recv_3d(int* target, int other_rank, int* new_size, int* from, int* to) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int SUB_NK = to[2] - from[2];
    int NI_NEW = new_size[0];
    int NJ_NEW = new_size[1];
    int NK_NEW = new_size[2];
    int* buffer = new int[SUB_NI*SUB_NJ*SUB_NK];
    MPI_Recv(&(buffer[0]), SUB_NI*SUB_NJ*SUB_NK, MPI_INT, other_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp parallel for collapse(4)
    for (int i = 0; i < SUB_NI; i++)
    {
        for (int j = 0; j < SUB_NJ; j++)
        {
            memcpy(target + (i+from[0])*NJ_NEW*NK_NEW + (j+from[1])*NK_NEW + from[2], buffer + i*SUB_NJ*SUB_NK + j*SUB_NK, sizeof(int)*SUB_NK);
        }
    }
}