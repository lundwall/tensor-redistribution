#include <omp.h>
#include <mpi.h>
#include <stdlib.h>

void send_2d(int* source, int other_rank, int* current_size, int* from, int* to, int num_chunks, MPI_Request *request, int* buffer) {
	if (num_chunks > 10)
	{
		std::cout << "Cannot have more than 10 chunks" << std::endl;
		std::exit(EXIT_FAILURE);
	}
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int NI = current_size[0];
    int NJ = current_size[1];

    for (int c = 0; c < num_chunks; ++c) { request[c] = MPI_REQUEST_NULL; }

    int chunk_size = (int) SUB_NI / num_chunks;
    int rem = SUB_NI % num_chunks;
    for (int c = 0; c < num_chunks; ++c)
    {
        int start_SUB_NI = c*chunk_size + std::min(c, rem);
        int end_SUB_NI = (c + 1)*chunk_size + std::min(c + 1, rem);
        #pragma omp parallel for
        for (int i = start_SUB_NI; i < end_SUB_NI; i++)
        {
            memcpy(buffer + i*SUB_NJ, source + (i+from[0])*NJ + from[1], sizeof(int)*SUB_NJ);
        }
        MPI_Isend(buffer + start_SUB_NI*SUB_NJ, (end_SUB_NI - start_SUB_NI)*SUB_NJ, MPI_INT, other_rank, 0, MPI_COMM_WORLD, &request[c]);
    }
}

void recv_2d(int* target, int other_rank, int* new_size, int* from, int* to, int num_chunks) {
    int SUB_NI = to[0] - from[0];
    int SUB_NJ = to[1] - from[1];
    int NI_NEW = new_size[0];
    int NJ_NEW = new_size[1];
    int* buffer = new int[SUB_NI*SUB_NJ];

    int chunk_size = (int) SUB_NI / num_chunks;
    int rem = SUB_NI % num_chunks;
    MPI_Recv(&(buffer[0]), (chunk_size + std::min(1, rem))*SUB_NJ, MPI_INT, other_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int c = 1; c < num_chunks; ++c)
    {
        int start_SUB_NI = c*chunk_size + std::min(c, rem);
        int end_SUB_NI = (c + 1)*chunk_size + std::min(c + 1, rem);
        MPI_Request request;
        MPI_Irecv(buffer + start_SUB_NI*SUB_NJ, (end_SUB_NI - start_SUB_NI)*SUB_NJ, MPI_INT, other_rank, 0, MPI_COMM_WORLD, &request);
        #pragma omp parallel for
        for (int i = (c-1)*chunk_size + std::min(c-1, rem); i < start_SUB_NI; i++)
        {
            memcpy(target + (i+from[0])*NJ_NEW + from[1], buffer + i*SUB_NJ, sizeof(int)*SUB_NJ);
        }
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    #pragma omp parallel for
    for (int i = (num_chunks - 1)*chunk_size + std::min(num_chunks - 1, rem); i < SUB_NI; i++)
    {
        memcpy(target + (i+from[0])*NJ_NEW + from[1], buffer + i*SUB_NJ, sizeof(int)*SUB_NJ);
    }

    delete[] buffer;
}