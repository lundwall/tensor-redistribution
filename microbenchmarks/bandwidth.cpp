#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

void zero(char *buf, size_t size)
{
    size_t my_start, my_size;

    if (omp_in_parallel())
    {
        int id = omp_get_thread_num();
        int num = omp_get_num_threads();

        my_start = (id*size)/num;
        my_size = ((id+1)*size)/num - my_start;
    }
    else
    {
        my_start = 0;
        my_size = size;
    }

    memset(buf + my_start, 0, my_size);
}

int main (void)
{
    char *buf;
    size_t size = 1L << 31; // 2 GiB
    double tmr;

    buf = (char*)malloc(size);

    // Touch
    tmr = -omp_get_wtime();
    #pragma omp parallel
    {
        zero(buf, size);
    }
    tmr += omp_get_wtime();
    printf("Touch:   %.3f MB/s\n", size/(1.e+6*tmr));

    // Rewrite
    tmr = -omp_get_wtime();
    #pragma omp parallel
    {
        zero(buf, size);
    }
    tmr += omp_get_wtime();
    printf("Rewrite: %.3f MB/s\n", size/(1.e+6*tmr));

    free(buf);

    return 0;
}