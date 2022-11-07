#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
// #include <liblsb.h>

#define CALIBRATE
#define NUM_RUNS 1
#define TIME_REQUIRED 0.5

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

    char* processor_name = new char[256]; int len_processor_name = 0;
    MPI_Get_processor_name(processor_name, &len_processor_name);
    std::cout << processor_name << std::endl;

    int **originalArray = new int*[N/P];
    int *originalBuffer = new int[N*N/P];
    for(int i = 0; i < N/P; ++i)
    {
        originalArray[i] = originalBuffer + N * i;
    }
    int **newArray = new int*[N];
    int *newBuffer = new int[N*N/P];
    for(int i = 0; i < N; ++i)
    {
        newArray[i] = newBuffer + N/P * i;
    }

    int send_array_size[2] = {N/P, N};
    int recv_array_size[2] = {N, N/P};
    int subarray_size[2] = {N/P, N/P};

    for (int i = 0; i < N/P; i++)
        for (int j = 0; j < N; j++)
            originalArray[i][j] = (N*N/P)*rank+i*N+j;

    double start, end;
    int num_runs;
    num_runs = NUM_RUNS;
    MPI_Request* sendreq = new MPI_Request[(P-1)*(N/P)];

#ifdef CALIBRATE
    while (num_runs < (1 << 14)) {
        start = MPI_Wtime();
        for (int k = 0; k < num_runs; ++k) {
            int count = 0;
            for (int i = 0; i < P; i++) {
                if (rank != i) {
                    for (int j = 0; j < N/P; j++) {
                        MPI_Isend(originalBuffer + j*N + i*(N/P), N/P, MPI_INT, i, 0, MPI_COMM_WORLD, &sendreq[count++]);
                    }
                }
            }

            for (int i = 0; i < P; i++) {
                if (rank != i) {
                    for (int j = 0; j < N/P; j++) {
                        MPI_Recv(newBuffer + j*(N/P) + i*(N/P)*(N/P), N/P, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }

            MPI_Waitall((P-1)*(N/P), sendreq, MPI_STATUSES_IGNORE);
        }
        end = MPI_Wtime();

        double ttime = (double)(end - start);

        if (ttime > TIME_REQUIRED) break;
        num_runs *= 2;
    }
#endif

    start = MPI_Wtime();
    for(int k = 0; k < num_runs; ++k) {
        int count = 0;
        for (int i = 0; i < P; i++) {
            if (rank != i) {
                for (int j = 0; j < N/P; j++) {
                    MPI_Isend(originalBuffer + j*N + i*(N/P), N/P, MPI_INT, i, 0, MPI_COMM_WORLD, &sendreq[count++]);
                }
            }
        }

        for (int i = 0; i < P; i++) {
            if (rank != i) {
                for (int j = 0; j < N/P; j++) {
                    MPI_Recv(newBuffer + j*(N/P) + i*(N/P)*(N/P), N/P, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                }
            }
        }

        MPI_Waitall((P-1)*(N/P), sendreq, MPI_STATUSES_IGNORE);
    }
    end = MPI_Wtime();
    std::cout << rank << " time : " << (end - start) / (double)num_runs * 1000000.0 << "us"<< std::endl;

    // self-copy
    for (int i = 0; i < N/P; i++) {
        memcpy(&newArray[(N / P) * rank + i][0], &originalArray[i][(N / P) * rank], sizeof(int) * (N / P));
    }

    MPI_Finalize();
    
    delete[] originalBuffer; originalBuffer = nullptr;
    delete[] originalArray; originalArray = nullptr;
    delete[] newBuffer; newBuffer = nullptr;
    delete[] newArray; newArray = nullptr;
    delete[] sendreq; sendreq = nullptr;
}
