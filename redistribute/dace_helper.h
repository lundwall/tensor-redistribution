#ifndef _DACE_HELPER_H
#define _DACE_HELPER_H
#include "mpi.h"
#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>

template <typename T, int VECLEN, int ALIGNED, int N>
struct CopyNDDynamic {
    struct Dynamic{

        static void Copy(const T *src, T *dst, const int &copydim, const int &src_stride, const int &dst_stride)
        {
            if (N == 1 && src_stride == 1 && dst_stride == 1) {
                memcpy(dst, src, copydim * sizeof(T) * VECLEN);
                return;
            }
        }

        template<typename ...Args>
        static void Copy(const T *src, T *dst, const int &copydim, const int &src_stride, const int &dst_stride,
                         const Args &... otherdims) {
            static_assert(sizeof...(otherdims) == (N - 1) * 3, "Dimensionality mismatch in dynamic copy");
            // Memcpy specialization
            if (N == 1 && src_stride == 1 && dst_stride == 1) {
                memcpy(dst, src, copydim * sizeof(T) * VECLEN);
                return;
            }

            for (int i = 0; i < copydim; ++i) {
                CopyNDDynamic<T, VECLEN, ALIGNED, N - 1>::Dynamic::Copy(src + i * src_stride, dst + i * dst_stride, otherdims...);
            }
        }
    };
};

int int_ceil(int x, int y)
{
    return (x + y - 1) / y;
}

int get_cart_rank(int grid_length, const int* grid, const int* coords) {
    int rank = coords[0];
    for (auto i = 1; i < grid_length; ++i) {
        rank *= grid[i];
        rank += coords[i];
    }
    return rank;
}
#endif
