#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <time.h>

#define RUNS 10

int main(int argc, char** argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    LSB_Init("test_manual_2d", 0);

    int P = size;
    int N = atoi(argv[1]);

    if (N % P != 0 || N < P)
    {
        printf("Number of processes should divide array length\n");
        return 0;
    }

    LSB_Set_Rparam_int("rank", rank);
    LSB_Set_Rparam_int("N", N);
    LSB_Set_Rparam_int("runs", RUNS);

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

    for (int i = 0; i < N/P; i++)
        for (int j = 0; j < N; j++)
            originalArray[i][j] = (N*N/P)*rank+i*N+j;

    MPI_Request* sendreq = new MPI_Request[P-1];

    srand(time(NULL));

    for(int k = 0; k < RUNS; ++k) {
        int count = 0;

        // Manual packing
        int **original_packs = new int*[P];
        int *original_pack = new int[N*(N/P)];
        for (int i = 0; i < P; ++i) {
            if (rank != i) {
                original_packs[i] = original_pack + (N/P) * (N/P) * i;
                for (int j = 0; j < N/P; ++j) {
                    memcpy(&original_packs[i][j * (N/P)], &originalArray[j][i * (N/P)], sizeof(int)*(N/P));
                }
            }
        }
        int **new_packs = new int*[P];
        int *new_pack = new int[N*(N/P)];
        for (int i = 0; i < P; ++i) {
            if (rank != i) {
                new_packs[i] = new_pack + (N/P) * (N/P) * i;
            }
        }

        LSB_Res();

        // Send
        for (int i = 0; i < P; i++) {
            if (rank != i) {
                MPI_Isend(&original_packs[i][0], (N/P)*(N/P), MPI_INT, i, 0, MPI_COMM_WORLD, &sendreq[count++]);
            }
        }

        // Receive
        for (int i = 0; i < P; i++) {
            if (rank != i) {
                MPI_Recv(&new_packs[i][0], (N/P)*(N/P), MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        LSB_Rec(k);

        for (int i = 0; i < P; i++) {
            if (rank != i) {
                for (int j = 0; j < N/P; j++) {
                    memcpy(&newArray[i * (N/P) + j][0], &new_packs[i][j * (N/P)], sizeof(int) * (N/P));
                }
            }
        }

        // Self-copy
        for (int i = 0; i < N/P; i++) {
            memcpy(&newArray[(N / P) * rank + i][0], &originalArray[i][(N / P) * rank], sizeof(int) * (N / P));
        }

        MPI_Waitall(P-1, sendreq, MPI_STATUSES_IGNORE);

        delete[] original_packs; original_packs = nullptr;
        delete[] original_pack; original_pack = nullptr;
        delete[] new_packs; new_packs = nullptr;
        delete[] new_pack; new_pack = nullptr;
    }

    LSB_Finalize();
    MPI_Finalize();
    
    delete[] originalBuffer; originalBuffer = nullptr;
    delete[] originalArray; originalArray = nullptr;
    delete[] newBuffer; newBuffer = nullptr;
    delete[] newArray; newArray = nullptr;
    delete[] sendreq; sendreq = nullptr;
}
