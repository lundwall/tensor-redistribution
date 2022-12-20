#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

void send_5d(int* source, int other_rank, int* current_size, int* from, int* to, MPI_Request *request, int* buffer) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int SUB_NK = to[2] - from[2];
    int SUB_NL = to[3] - from[3];
    int SUB_NM = to[4] - from[4];
    int NI = current_size[0];
    int NJ = current_size[1];
    int NK = current_size[2];
    int NL = current_size[3];
    int NM = current_size[4];
    #pragma omp parallel for collapse(4)
    for (int i = 0; i < SUB_NI; i++)
    {
        for (int j = 0; j < SUB_NJ; j++)
        {
            for (int k = 0; k < SUB_NK; k++)
            {
                for (int l = 0; l < SUB_NL; l++)
                {
                    memcpy(buffer + i*SUB_NJ*SUB_NK*SUB_NL*SUB_NM + j*SUB_NK*SUB_NL*SUB_NM + k*SUB_NL*SUB_NM + l*SUB_NM, source + (i+from[0])*NJ*NK*NL*NM + (j+from[1])*NK*NL*NM + (k+from[2])*NL*NM + (l+from[3])*NM + (from[4]), sizeof(int)*SUB_NM);
                }
            }
        }
    }
    MPI_Isend(buffer, SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM, MPI_INT, other_rank, 0, MPI_COMM_WORLD, &request[0]);
}

void recv_5d(int* target, int other_rank, int* new_size, int* from, int* to) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int SUB_NK = to[2] - from[2];
    int SUB_NL = to[3] - from[3];
    int SUB_NM = to[4] - from[4];
    int NI_NEW = new_size[0];
    int NJ_NEW = new_size[1];
    int NK_NEW = new_size[2];
    int NL_NEW = new_size[3];
    int NM_NEW = new_size[4];
    int* buffer = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];
    MPI_Recv(&(buffer[0]), SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM, MPI_INT, other_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #pragma omp parallel for collapse(4)
    for (int i = 0; i < SUB_NI; i++)
    {
        for (int j = 0; j < SUB_NJ; j++)
        {
            for (int k = 0; k < SUB_NK; k++)
            {
                for (int l = 0; l < SUB_NL; l++)
                {
                    memcpy(target + (i+from[0])*NJ_NEW*NK_NEW*NL_NEW*NM_NEW + (j+from[1])*NK_NEW*NL_NEW*NM_NEW + (k+from[2])*NL_NEW*NM_NEW + (l+from[3])*NM_NEW + (from[4]), buffer + i*SUB_NJ*SUB_NK*SUB_NL*SUB_NM + j*SUB_NK*SUB_NL*SUB_NM + k*SUB_NL*SUB_NM + l*SUB_NM, sizeof(int)*SUB_NM);
                }
            }
        }
    }
}