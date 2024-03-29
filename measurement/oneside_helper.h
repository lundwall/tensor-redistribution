#ifndef ONESIDE_HELPER_H
#define ONESIDE_HELPER_H

#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

void prepare_send_buffer(int* source, int* current_size, int* from, int* to, int* buffer) {
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
}

void unpack_recv_buffer(int* target, int* new_size, int* from, int* to, int* buffer) {
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

void one_sided_send_with_chunks(int* source, int* current_size, int* from, int* to, int* buffer, int num_chunks, MPI_Win& window, int other_rank) {
    if (num_chunks > 10)
    {
        std::cout << "Cannot have more than 10 chunks" << std::endl;
        std::exit(EXIT_FAILURE);
    }
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

    int chunk_size = (int) SUB_NI / num_chunks;
    int rem = SUB_NI % num_chunks;
    for (int c = 0; c < num_chunks; ++c)
    {
        int start_SUB_NI = c*chunk_size + std::min(c, rem);
        int end_SUB_NI = (c + 1)*chunk_size + std::min(c + 1, rem);
        #pragma omp parallel for collapse(4)
        for (int i = start_SUB_NI; i < end_SUB_NI; i++)
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
        int size = (end_SUB_NI - start_SUB_NI)*SUB_NJ*SUB_NK*SUB_NL*SUB_NM;
        MPI_Put(buffer + start_SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM, size, MPI_INT, other_rank, start_SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM, size, MPI_INT, window);
    }
}

#endif
