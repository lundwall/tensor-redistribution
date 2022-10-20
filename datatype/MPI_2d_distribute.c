#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int** allocate2dArray(int row_num, int col_num)
{
	int *array_begin = (int*)malloc(sizeof(int) * row_num * col_num);
	int **array = (int**)malloc(sizeof(int*) * row_num);

	for (int i = 0; i < row_num; i++)
	{
		array[i] = array_begin + col_num*i;
	}
	return array;
}

int main(int argc, char** argv)
{
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int P = size;
	int N = atoi(argv[1]);

	if (N % P != 0 || N < P)
	{
		printf("Number of processes should divide array length\n");
		return 0;
	}

	int **originalArray = allocate2dArray(N/P, N);
	int **newArray = allocate2dArray(N, N/P);

	int send_array_size[2] = {N/P, N};
	int recv_array_size[2] = {N, N/P};
	int subarray_size[2] = {N/P, N/P};

	for (int i = 0; i < N/P; i++)
		for (int j = 0; j < N; j++)
			originalArray[i][j] = (N*N/P)*rank+i*N+j;

	MPI_Datatype* sendtype = (MPI_Datatype*)malloc(sizeof(MPI_Datatype)*P);
	MPI_Datatype* recvtype = (MPI_Datatype*)malloc(sizeof(MPI_Datatype)*P);

	for (int i = 0; i < P; i++)
	{
		int send_start[2] = {0, (N/P)*i};
		MPI_Type_create_subarray(2, send_array_size, subarray_size, send_start, MPI_ORDER_C, MPI_INT, &sendtype[i]);
		MPI_Type_commit(&sendtype[i]);
		if (rank != i)
		{
			MPI_Send(&(originalArray[0][0]), 1, sendtype[i], i, 0, MPI_COMM_WORLD);
		}
	}


	for (int i = 0; i < P; i++)
	{
		int recv_start[2] = {(N/P)*i, 0};
		MPI_Type_create_subarray(2, recv_array_size, subarray_size, recv_start, MPI_ORDER_C, MPI_INT, &recvtype[i]);
		MPI_Type_commit(&recvtype[i]);
		if (rank != i)
		{
			MPI_Recv(&(newArray[0][0]), 1, recvtype[i], i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	// self-copy
	for (int i = 0; i < N/P; i++)
	{
		memcpy(&newArray[(N/P)*rank+i][0], &originalArray[i][(N/P)*rank], sizeof(int)*(N/P));
	}

	if (rank == 1)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N/P; j++)
			{
				printf("%d_%d ", rank, newArray[i][j]);
			}
			printf("\n");
		}
	}

	MPI_Finalize();
}