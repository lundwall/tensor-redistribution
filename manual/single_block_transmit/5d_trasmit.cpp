#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <time.h>
#include <omp.h>
#include <string>
#define NI 40
#define NJ 40
#define NK 60
#define NL 80
#define NM 60

#define NI_NEW 40
#define NJ_NEW 40
#define NK_NEW 80
#define NL_NEW 60
#define NM_NEW 60

#define SUB_NI 30
#define SUB_NJ 30
#define SUB_NK 30
#define SUB_NL 30
#define SUB_NM 30

#define RUNS 100
#define CHUNK_SIZE 4050000 // 5*30*30*30, 1/6 of all the data
#define NUM_CHUNKS SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM/CHUNK_SIZE

int main(int argc, char** argv)
{
    int size, rank;
    std::string name_string = "5d_transmit_manual"+std::string(std::getenv("OMP_NUM_THREADS"));
    const char* liblsb_fname = name_string.c_str();
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    LSB_Init(liblsb_fname, 0);
    LSB_Set_Rparam_int("rank", rank);
    int max_threads = omp_get_max_threads();
    LSB_Set_Rparam_int("threads", max_threads);

    // size should be 2!
    if (size != 2)
    {
        printf("When testing only 1 block transmission, only 2 processors are needed\n");
        return 0;
    }

    char* processor_name = new char[256]; int len_processor_name = 0;
    MPI_Get_processor_name(processor_name, &len_processor_name);
    std::cout << processor_name << std::endl;

    int* originalArray = new int[NI*NJ*NK*NL*NM];
    int* newArray = new int[NI_NEW*NJ_NEW*NK_NEW*NL_NEW*NM_NEW];
    int* sendArray = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];
    int* recvArray = new int[SUB_NI*SUB_NJ*SUB_NK*SUB_NL*SUB_NM];
    MPI_Request* sendreq = new MPI_Request[1];
    MPI_Request* recvreq = new MPI_Request[1];

    for (int i = 0; i < NI*NJ*NK*NL*NM; i++)
        originalArray[i] = 1;
    for (int i = 0; i < NI_NEW*NJ_NEW*NK_NEW*NL_NEW*NM_NEW; i++)
        newArray[i] = 0;

    srand(time(NULL));

    for (int r = 0; r < RUNS; ++r) {
        int count = 0;

        if (rank == 0)
        {
            LSB_Res();
            for (int chunk = 0; chunk < NUM_CHUNKS; chunk++)
            {
                if (chunk > 0) {
                    MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
                }
                for (int i = 0; i < CHUNK_SIZE / (SUB_NJ*SUB_NK*SUB_NL*SUB_NM); i++)
                {
                    for (int j = 0; j < SUB_NJ; j++)
                    {
                        for (int k = 0; k < SUB_NK; k++)
                        {
                            #pragma omp parallel
                            #pragma omp single
                            for (int l = 0; l < SUB_NL; l++)
                            {
                                #pragma omp task
                                memcpy(sendArray + chunk*CHUNK_SIZE + i*SUB_NJ*SUB_NK*SUB_NL*SUB_NM + j*SUB_NK*SUB_NL*SUB_NM + k*SUB_NL*SUB_NM + l*SUB_NM, originalArray + chunk*(CHUNK_SIZE/(SUB_NJ*SUB_NK*SUB_NL*SUB_NM))*(NJ*NK*NL*NM) + i*NJ*NK*NL*NM + j*NK*NL*NM + k*NL*NM + l*NM, sizeof(int)*SUB_NM);
                            }
                        }
                    }
                }
                MPI_Isend(&(sendArray[chunk*CHUNK_SIZE]), CHUNK_SIZE, MPI_INT, 1, chunk, MPI_COMM_WORLD, sendreq);
            }
            MPI_Waitall(1, sendreq, MPI_STATUSES_IGNORE);
            LSB_Rec(r);
        }

        if (rank == 1)
        {
            LSB_Res();
            MPI_Irecv(&(recvArray[0]), CHUNK_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD, recvreq);

            for (int chunk = 0; chunk < NUM_CHUNKS; chunk++)
            {
                MPI_Waitall(1, recvreq, MPI_STATUSES_IGNORE);
                if (chunk != NUM_CHUNKS - 1)
                {
                    MPI_Irecv(&(recvArray[(chunk+1)*CHUNK_SIZE]), CHUNK_SIZE, MPI_INT, 0, chunk+1, MPI_COMM_WORLD, recvreq);
                }
                for (int i = 0; i < CHUNK_SIZE / (SUB_NJ*SUB_NK*SUB_NL); i++)
                {
                    for (int j = 0; j < SUB_NJ; j++)
                    {
                        for (int k = 0; k < SUB_NK; k++)
                        {
                            #pragma omp parallel
                            #pragma omp single
                            for (int l = 0; l < SUB_NL; l++)
                            {
                                #pragma omp task
                                memcpy(newArray + chunk*(CHUNK_SIZE / (SUB_NJ*SUB_NK*SUB_NL*SUB_NM)) + i*NJ_NEW*NK_NEW*NL_NEW*NM_NEW + j*NK_NEW*NL_NEW*NM_NEW + k*NL_NEW*NM_NEW + l*NM_NEW, recvArray + chunk*CHUNK_SIZE + i*SUB_NJ*SUB_NK*SUB_NL*SUB_NM + j*SUB_NK*SUB_NL*SUB_NM + k*SUB_NL*SUB_NM + l*SUB_NM, sizeof(int)*SUB_NM);
                            }
                        }
                    }
                }
            }
            LSB_Rec(r);
        }
    }

    LSB_Finalize();
    MPI_Finalize();

    // if (rank == 1)
    // {
    //     for (int i = 0; i < NI_NEW; i++)
    //     {
    //         for (int j = 0; j < NJ_NEW; j++)
    //         {
    //             for (int k = 0; k < NK_NEW; k++)
    //             {
    //                 for (int l = 0; l < NL_NEW; l++) 
    //                 {
    //                     for (int m = 0; m < NM_NEW; m++)
    //                     {
    //                         std::cout << newArray[i*NJ_NEW*NK_NEW*NL_NEW*NM_NEW+j*NK_NEW*NL_NEW*NM_NEW+k*NL_NEW*NM_NEW+l*NM_NEW+m] << " ";
    //                     }
    //                     std::cout << std::endl;
    //                 }
    //                 std::cout << std::endl << std::endl;
    //             }
    //             std::cout << std::endl << std::endl << std::endl;
    //         }
    //         std::cout << std::endl << std::endl << std::endl << std::endl;
    //     }
    // }

    delete[] originalArray; originalArray = nullptr;
    delete[] newArray; newArray = nullptr;
}
