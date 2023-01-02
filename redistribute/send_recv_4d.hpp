#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

void send_4d(int* source, int other_rank, int* current_size, int* from, int* to, int num_chunks, MPI_Request *request, int* buffer) {
	if (num_chunks > 10)
	{
		std::cout << "Cannot have more than 10 chunks" << std::endl;
		std::exit(EXIT_FAILURE);
	}
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int SUB_NK = to[2] - from[2];
    int SUB_NL = to[3] - from[3];
    int NI = current_size[0];
    int NJ = current_size[1];
    int NK = current_size[2];
    int NL = current_size[3];

    for (int c = 0; c < num_chunks; ++c) { request[c] = MPI_REQUEST_NULL; }

    int chunk_size = (int) SUB_NI / num_chunks;
    int rem = SUB_NI % num_chunks;
    for (int c = 0; c < num_chunks; ++c)
    {
        int start_SUB_NI = c*chunk_size + std::min(c, rem);
        int end_SUB_NI = (c + 1)*chunk_size + std::min(c + 1, rem);
        #pragma omp parallel for collapse(3)
        for (int i = start_SUB_NI; i < end_SUB_NI; i++)
        {
            for (int j = 0; j < SUB_NJ; j++)
            {
                for (int k = 0; k < SUB_NK; k++)
                {
                    memcpy(buffer + i*SUB_NJ*SUB_NK*SUB_NL + j*SUB_NK*SUB_NL + k*SUB_NL , source + (i+from[0])*NJ*NK*NL + (j+from[1])*NK*NL+ (k+from[2])*NL + from[3], sizeof(int)*SUB_NL);
                }
            }
        }
        MPI_Isend(buffer + start_SUB_NI*SUB_NJ*SUB_NK*SUB_NL, (end_SUB_NI - start_SUB_NI)*SUB_NJ*SUB_NK*SUB_NL, MPI_INT, other_rank, 0, MPI_COMM_WORLD, &request[c]);
    }
}

void recv_4d(int* target, int other_rank, int* new_size, int* from, int* to, int num_chunks) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int SUB_NK = to[2] - from[2];
    int SUB_NL = to[3] - from[3];
    int NI_NEW = new_size[0];
    int NJ_NEW = new_size[1];
    int NK_NEW = new_size[2];
    int NL_NEW = new_size[3];
    int* buffer = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL];

    int chunk_size = (int) SUB_NI / num_chunks;
    int rem = SUB_NI % num_chunks;
    MPI_Recv(&(buffer[0]), (chunk_size + std::min(1, rem))*SUB_NJ*SUB_NK*SUB_NL, MPI_INT, other_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int c = 1; c < num_chunks; ++c)
    {
        int start_SUB_NI = c*chunk_size + std::min(c, rem);
        int end_SUB_NI = (c + 1)*chunk_size + std::min(c + 1, rem);
        MPI_Request request;
        MPI_Irecv(buffer + start_SUB_NI*SUB_NJ*SUB_NK*SUB_NL, (end_SUB_NI - start_SUB_NI)*SUB_NJ*SUB_NK*SUB_NL, MPI_INT, other_rank, 0, MPI_COMM_WORLD, &request);
        #pragma omp parallel for collapse(3)
        for (int i = (c-1)*chunk_size + std::min(c-1, rem); i < start_SUB_NI; i++)
        {
            for (int j = 0; j < SUB_NJ; j++)
            {
                for (int k = 0; k < SUB_NK; k++)
                {
                    memcpy(target + (i+from[0])*NJ_NEW*NK_NEW*NL_NEW + (j+from[1])*NK_NEW*NL_NEW + (k+from[2])*NL_NEW + from[3], buffer + i*SUB_NJ*SUB_NK*SUB_NL + j*SUB_NK*SUB_NL + k*SUB_NL, sizeof(int)*SUB_NL);
                }
            }
        }
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    #pragma omp parallel for collapse(3)
    for (int i = (num_chunks - 1)*chunk_size + std::min(num_chunks - 1, rem); i < SUB_NI; i++)
    {
        for (int j = 0; j < SUB_NJ; j++)
        {
            for (int k = 0; k < SUB_NK; k++)
            {
                    memcpy(target + (i+from[0])*NJ_NEW*NK_NEW*NL_NEW + (j+from[1])*NK_NEW*NL_NEW + (k+from[2])*NL_NEW + from[3], buffer + i*SUB_NJ*SUB_NK*SUB_NL + j*SUB_NK*SUB_NL + k*SUB_NL, sizeof(int)*SUB_NL);
            }
        }
    }

    delete[] buffer;
}