#include <mpi.h>
#include <iostream>
#include <cstring>

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

struct matrix_1d_1d_t {
    MPI_Comm __pgrid_0_comm;
    MPI_Group __pgrid_0_group;
    int __pgrid_0_coords[2];
    int __pgrid_0_dims[2];
    int __pgrid_0_rank;
    int __pgrid_0_size;
    bool __pgrid_0_valid;
    MPI_Comm __pgrid_1_comm;
    MPI_Group __pgrid_1_group;
    int __pgrid_1_coords[2];
    int __pgrid_1_dims[2];
    int __pgrid_1_rank;
    int __pgrid_1_size;
    bool __pgrid_1_valid;
    MPI_Datatype __subarray_0;
    int* __subarray_0_counts;
    int* __subarray_0_displs;
    MPI_Datatype __subarray_1;
    int* __subarray_1_counts;
    int* __subarray_1_displs;
    MPI_Datatype __rdistrarray_0;
    int __rdistrarray_0_sends;
    MPI_Datatype* __rdistrarray_0_send_types;
    int* __rdistrarray_0_dst_ranks;
    int __rdistrarray_0_recvs;
    MPI_Datatype* __rdistrarray_0_recv_types;
    int* __rdistrarray_0_src_ranks;
    int __rdistrarray_0_self_copies;
    int* __rdistrarray_0_self_src;
    int* __rdistrarray_0_self_dst;
    int* __rdistrarray_0_self_size;
};

void __program_matrix_1d_1d_internal(matrix_1d_1d_t *__state, int * __restrict__ A, int * __restrict__ __return, int P, int m)
{

    {
        int a_grid;
        int b_grid;
        int a_arr;
        int b_arr;
        int rdistr;

        {
            int* _inp_buffer = &A[0];
            int* _out_buffer = __return;

            ///////////////////

            int myrank;
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            MPI_Request* req = new MPI_Request[__state->__rdistrarray_0_sends];
            MPI_Status* status = new MPI_Status[__state->__rdistrarray_0_sends];
            MPI_Status recv_status;
            if (__state->__pgrid_0_valid) {
                for (auto __idx = 0; __idx < __state->__rdistrarray_0_sends; ++__idx) {
                    // printf("(__subarray_0 -> __subarray_1) I am rank %d and I send to %d\n", myrank, __state->__rdistrarray_0_dst_ranks[__idx]);
                    // fflush(stdout);
                    MPI_Isend(_inp_buffer, 1, __state->__rdistrarray_0_send_types[__idx], __state->__rdistrarray_0_dst_ranks[__idx], 0, MPI_COMM_WORLD, &req[__idx]);
                }
            }
            if (__state->__pgrid_1_valid) {
                for (auto __idx = 0; __idx < __state->__rdistrarray_0_self_copies; ++__idx) {
                    // printf("(__subarray_0 -> __subarray_1) I am rank %d and I self-copy\n", myrank);
                    // fflush(stdout);
                    int __inp_s0 = __state->__rdistrarray_0_self_src[__idx * 2 + 0];
                    int __inp_s1 = __state->__rdistrarray_0_self_src[__idx * 2 + 1];

                    int __out_s0 = __state->__rdistrarray_0_self_dst[__idx * 2 + 0];
                    int __out_s1 = __state->__rdistrarray_0_self_dst[__idx * 2 + 1];

                    CopyNDDynamic<int, 1, false, 2>::Dynamic::Copy(
                            _inp_buffer + (((P * __inp_s0) * m) + __inp_s1), _out_buffer + ((__out_s0 * m) + __out_s1), __state->__rdistrarray_0_self_size[__idx * 2 + 0], P*m, m, __state->__rdistrarray_0_self_size[__idx * 2 + 1], 1, 1
                    );
                }
                for (auto __idx = 0; __idx < __state->__rdistrarray_0_recvs; ++__idx) {
                    // printf("(__subarray_0 -> __subarray_1) I am rank %d and I receive from %d\n", myrank, __state->__rdistrarray_0_src_ranks[__idx]);
                    // fflush(stdout);
                    MPI_Recv(_out_buffer, 1, __state->__rdistrarray_0_recv_types[__idx], __state->__rdistrarray_0_src_ranks[__idx], 0, MPI_COMM_WORLD, &recv_status);
                }
            }
            if (__state->__pgrid_0_valid) {
                MPI_Waitall(__state->__rdistrarray_0_sends, req, status);
                delete[] req;
                delete[] status;
            }
            // printf("I am rank %d and I finished the redistribution __subarray_0 -> __subarray_1\n", myrank);
            // fflush(stdout);


            ///////////////////

        }
        {
            int __out;

            ///////////////////
            // Tasklet code (__pgrid_0)
            ///////////////////

            a_grid = __out;
        }
        {
            int __out;

            ///////////////////
            // Tasklet code (__pgrid_1)
            ///////////////////

            b_grid = __out;
        }
        {
            int __out;

            ///////////////////
            // Tasklet code (__subarray_0)
            ///////////////////

            a_arr = __out;
        }
        {
            int __out;

            ///////////////////
            // Tasklet code (__subarray_1)
            ///////////////////

            b_arr = __out;
        }
        {
            int __out;

            ///////////////////
            // Tasklet code (__rdistrarray_0)
            ///////////////////

            rdistr = __out;
        }

    }
}

void __program_matrix_1d_1d(matrix_1d_1d_t *__state, int * __restrict__ A, int * __restrict__ __return, int P, int m)
{
    __program_matrix_1d_1d_internal(__state, A, __return, P, m);
}

// the total data m * P, m * P
// grid_0: P * 1, block_size_0 is (m, m * P)
// grid_1: 1 * P, block_size_1 is (m * p, m)
matrix_1d_1d_t *__dace_init_matrix_1d_1d(int P, int m)
{
    int __result = 0;
    matrix_1d_1d_t *__state = new matrix_1d_1d_t;


    {  // Environment: MPI
        int t; MPI_Initialized(&t);  if (!t) MPI_Init(NULL, NULL);
    }
    __state->__pgrid_0_dims[0] = P;
    __state->__pgrid_0_dims[1] = 1;

    // create cartesian grid for pgrid0 and pgrid1
    int __pgrid_0_periods[2] = {0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, __state->__pgrid_0_dims, __pgrid_0_periods, 0, &__state->__pgrid_0_comm);
    if (__state->__pgrid_0_comm != MPI_COMM_NULL) {
        MPI_Comm_group(__state->__pgrid_0_comm, &__state->__pgrid_0_group);
        MPI_Comm_rank(__state->__pgrid_0_comm, &__state->__pgrid_0_rank);
        MPI_Comm_size(__state->__pgrid_0_comm, &__state->__pgrid_0_size);
        MPI_Cart_coords(__state->__pgrid_0_comm, __state->__pgrid_0_rank, 2, __state->__pgrid_0_coords);
        __state->__pgrid_0_valid = true;
    } else {
        __state->__pgrid_0_group = MPI_GROUP_NULL;
        __state->__pgrid_0_rank = MPI_PROC_NULL;
        __state->__pgrid_0_size = 0;
        __state->__pgrid_0_valid = false;
    }
    __state->__pgrid_1_dims[0] = 1;
    __state->__pgrid_1_dims[1] = P;

    int __pgrid_1_periods[2] = {0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, __state->__pgrid_1_dims, __pgrid_1_periods, 0, &__state->__pgrid_1_comm);
    if (__state->__pgrid_1_comm != MPI_COMM_NULL) {
        MPI_Comm_group(__state->__pgrid_1_comm, &__state->__pgrid_1_group);
        MPI_Comm_rank(__state->__pgrid_1_comm, &__state->__pgrid_1_rank);
        MPI_Comm_size(__state->__pgrid_1_comm, &__state->__pgrid_1_size);
        MPI_Cart_coords(__state->__pgrid_1_comm, __state->__pgrid_1_rank, 2, __state->__pgrid_1_coords);
        __state->__pgrid_1_valid = true;
    } else {
        __state->__pgrid_1_group = MPI_GROUP_NULL;
        __state->__pgrid_1_rank = MPI_PROC_NULL;
        __state->__pgrid_1_size = 0;
        __state->__pgrid_1_valid = false;
    }


    {
        __state->__rdistrarray_0_sends = 0;
        __state->__rdistrarray_0_recvs = 0;
        __state->__rdistrarray_0_self_copies = 0;
        int max_sends = 1;
        int max_recvs = 1;

        int kappa[2];
        int lambda[2];
        int xi[2];
        int pcoords[2];

        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        max_sends *= int_ceil(int(m) + int(P*m) - 1, int(P*m));
        max_recvs *= int_ceil(int(P*m) - 1, int(m)) + 1;

        max_sends *= int_ceil(int(P*m) + int(m) - 1, int(m));
        max_recvs *= int_ceil(int(m) - 1, int(P*m)) + 1;

        __state->__rdistrarray_0_send_types = new MPI_Datatype[max_sends];
        __state->__rdistrarray_0_dst_ranks = new int[max_sends];
        __state->__rdistrarray_0_recv_types = new MPI_Datatype[max_recvs];
        __state->__rdistrarray_0_src_ranks = new int[max_recvs];
        __state->__rdistrarray_0_self_src = new int[max_sends * 2];
        __state->__rdistrarray_0_self_dst = new int[max_sends * 2];
        __state->__rdistrarray_0_self_size = new int[max_sends * 2];

        if (__state->__pgrid_1_valid) {

            int sizes[2] = {P*m, m};
            int subsizes[2];
            int origin[2];

            xi[0] = (__state->__pgrid_1_coords[0] * int(P*m)) / int(m);
            lambda[0] = __state->__pgrid_1_coords[0] * int(P*m) % int(m);
            kappa[0] = int_ceil(int(P*m) + lambda[0], int(m));

            xi[1] = (__state->__pgrid_1_coords[1] * int(m)) / int(P*m);
            lambda[1] = __state->__pgrid_1_coords[1] * int(m) % int(P*m);
            kappa[1] = int_ceil(int(m) + lambda[1], int(P*m));

            int rem0 = P*m;
            for (auto idx0 = 0; idx0 < kappa[0]; ++idx0) {
                int actual_idx0 = 0;
                pcoords[actual_idx0] = xi[0] + idx0;
                int lo0 = (idx0 == 0 ? lambda[0] : 0);
                int uo0 = std::min(int(m), lo0 + rem0);
                subsizes[0] = uo0 - lo0;
                origin[0] = P*m - rem0;
                rem0 -= uo0 - lo0;

                int rem1 = m;
                for (auto idx1 = 0; idx1 < kappa[1]; ++idx1) {
                    int actual_idx1 = 1;
                    pcoords[actual_idx1] = xi[1] + idx1;
                    int lo1 = (idx1 == 0 ? lambda[1] : 0);
                    int uo1 = std::min(int(P*m), lo1 + rem1);
                    subsizes[1] = uo1 - lo1;
                    origin[1] = m - rem1;
                    rem1 -= uo1 - lo1;
                    int cart_rank = get_cart_rank(2, __state->__pgrid_0_dims, pcoords);
                    if (myrank == cart_rank) { // self-copy
                        __state->__rdistrarray_0_self_src[__state->__rdistrarray_0_self_copies * 2 + 0] = lo0;
                        __state->__rdistrarray_0_self_dst[__state->__rdistrarray_0_self_copies * 2 + 0] = origin[0];
                        __state->__rdistrarray_0_self_size[__state->__rdistrarray_0_self_copies * 2 + 0] = subsizes[0];

                        __state->__rdistrarray_0_self_src[__state->__rdistrarray_0_self_copies * 2 + 1] = lo1;
                        __state->__rdistrarray_0_self_dst[__state->__rdistrarray_0_self_copies * 2 + 1] = origin[1];
                        __state->__rdistrarray_0_self_size[__state->__rdistrarray_0_self_copies * 2 + 1] = subsizes[1];

                        __state->__rdistrarray_0_self_copies++;
                        // printf("(__subarray_0 -> __subarray_1) I am rank %d and I self-copy {I receive from %d%d (%d - %d) in (%d, %d) size (%d, %d)} \n", myrank, pcoords[0], pcoords[1], cart_rank, cart_rank, origin[0], origin[1], subsizes[0], subsizes[1]);
                    } else {
                        MPI_Type_create_subarray(2,  sizes, subsizes, origin, MPI_ORDER_C, MPI_INT, &__state->__rdistrarray_0_recv_types[__state->__rdistrarray_0_recvs]);
                        MPI_Type_commit(&__state->__rdistrarray_0_recv_types[__state->__rdistrarray_0_recvs]);
                        __state->__rdistrarray_0_src_ranks[__state->__rdistrarray_0_recvs] = cart_rank;
                        // printf("(__subarray_0 -> __subarray_1) I am rank %d and I receive from %d%d (%d - %d) in (%d, %d) size (%d, %d) \n", myrank, pcoords[0], pcoords[1], cart_rank, __state->__rdistrarray_0_src_ranks[__state->__rdistrarray_0_recvs], origin[0], origin[1], subsizes[0], subsizes[1]);
                        __state->__rdistrarray_0_recvs++;
                    }
                }}}
        if (__state->__pgrid_0_valid) {

            int sizes[2] = {m, P*m};
            int subsizes[2];
            int origin[2];

            // int_ceil(x, y) := (x + y - 1) / y
            // int_ceil(pcoord * sa - sb + 1, sb) = (pcoord * sa) / sb
            int lp0 = std::max(0, (__state->__pgrid_0_coords[0] * int(m)) / int(P*m)); // int_ceil(x, y) := (x + y - 1) / y
            int up0 = std::min(__state->__pgrid_1_dims[0], int_ceil((__state->__pgrid_0_coords[0] + 1) * int(m), int(P*m)));
            // printf("I am rank %d and I have 0-th bounds [%d, %d)\n", myrank, lp0, up0);
            for (auto idx0 = lp0; idx0 < up0; ++idx0) {
                int actual_idx0 = 0;

                xi[0] = (idx0 * int(P*m)) / int(m);
                lambda[0] = idx0 * int(P*m) % int(m);
                kappa[0] = int_ceil(int(P*m) + lambda[0], int(m));
                int idx0_dst = __state->__pgrid_0_coords[0] - xi[0];

                if (idx0_dst < 0 || idx0_dst >= kappa[0]) continue;
                int lo0 = (idx0_dst == 0 ? lambda[0] : 0);
                int uo0 = (idx0_dst == kappa[0] - 1 ? int(P*m) + lambda[0] - idx0_dst * int(m) : int(m));
                subsizes[0] = uo0 - lo0;
                origin[0] = lo0;
                pcoords[actual_idx0] = idx0;


                // int_ceil(x, y) := (x + y - 1) / y
                // int_ceil(pcoord * sa - sb + 1, sb) = (pcoord * sa) / sb
                int lp1 = std::max(0, (__state->__pgrid_0_coords[1] * int(P*m)) / int(m)); // int_ceil(x, y) := (x + y - 1) / y
                int up1 = std::min(__state->__pgrid_1_dims[1], int_ceil((__state->__pgrid_0_coords[1] + 1) * int(P*m), int(m)));
                // printf("I am rank %d and I have 1-th bounds [%d, %d)\n", myrank, lp1, up1);
                for (auto idx1 = lp1; idx1 < up1; ++idx1) {
                    int actual_idx1 = 1;

                    xi[1] = (idx1 * int(m)) / int(P*m);
                    lambda[1] = idx1 * int(m) % int(P*m);
                    kappa[1] = int_ceil(int(m) + lambda[1], int(P*m));
                    int idx1_dst = __state->__pgrid_0_coords[1] - xi[1];

                    if (idx1_dst < 0 || idx1_dst >= kappa[1]) continue;
                    int lo1 = (idx1_dst == 0 ? lambda[1] : 0);
                    int uo1 = (idx1_dst == kappa[1] - 1 ? int(m) + lambda[1] - idx1_dst * int(P*m) : int(P*m));
                    subsizes[1] = uo1 - lo1;
                    origin[1] = lo1;
                    pcoords[actual_idx1] = idx1;

                    int cart_rank = get_cart_rank(2, __state->__pgrid_1_dims, pcoords);

                    if (myrank != cart_rank) { // not self-copy
                        MPI_Type_create_subarray(2,  sizes, subsizes, origin, MPI_ORDER_C, MPI_INT, &__state->__rdistrarray_0_send_types[__state->__rdistrarray_0_sends]);
                        MPI_Type_commit(&__state->__rdistrarray_0_send_types[__state->__rdistrarray_0_sends]);
                        __state->__rdistrarray_0_dst_ranks[__state->__rdistrarray_0_sends] = cart_rank;
                        // printf("(__subarray_0 -> __subarray_1) I am rank %d and I send to %d%d (%d - %d) from (%d, %d) size (%d, %d)\n", myrank, pcoords[0], pcoords[1], cart_rank, __state->__rdistrarray_0_dst_ranks[__state->__rdistrarray_0_sends], origin[0], origin[1], subsizes[0], subsizes[1]);
                        __state->__rdistrarray_0_sends++;
                    }
                }}}}

    if (__result) {
        delete __state;
        return nullptr;
    }
    return __state;
}

void __dace_exit_matrix_1d_1d(matrix_1d_1d_t *__state)
{

    if (__state->__pgrid_0_valid) {
        MPI_Group_free(&__state->__pgrid_0_group);
        MPI_Comm_free(&__state->__pgrid_0_comm);
    }

    if (__state->__pgrid_1_valid) {
        MPI_Group_free(&__state->__pgrid_1_group);
        MPI_Comm_free(&__state->__pgrid_1_comm);
    }

    if (__state->__pgrid_0_valid) {
        for (auto __idx = 0; __idx < __state->__rdistrarray_0_sends; ++__idx) {
            MPI_Type_free(&__state->__rdistrarray_0_send_types[__idx]);
        }
    }
    delete[] __state->__rdistrarray_0_send_types;
    delete[] __state->__rdistrarray_0_dst_ranks;
    delete[] __state->__rdistrarray_0_recv_types;
    delete[] __state->__rdistrarray_0_src_ranks;
    delete[] __state->__rdistrarray_0_self_src;
    delete[] __state->__rdistrarray_0_self_dst;
    delete[] __state->__rdistrarray_0_self_size;

    {  // Environment: MPI
        // MPI_Finalize();
    }
    delete __state;
}


int main(int argc, char** argv)
{
    int P = 4;
    int m = 3;
    matrix_1d_1d_t* state = __dace_init_matrix_1d_1d(P, m);

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // the total data m * P, m * P
    // grid_0: P * 1, block_size_0 is (m, m * P)
    // grid_1: 1 * P, block_size_1 is (m * p, m)

    int *originalArray = new int [P * m * m];
    int *newArray = new int [P * m * m];
    for(int i = 0; i < m; ++i)
    {
        for(int j = 0; j < m * P; ++j)
        {
            originalArray[i * (m * P) + j] = (m * m * P) * myrank + i * m * P + j;
        }
    }
    __program_matrix_1d_1d(state, originalArray, newArray, P, m);

    if(myrank == 1) {
        for (int i = 0; i < P * m * m; ++i) {
            std::cout << newArray[i] << " ";
        }
    }

    __dace_exit_matrix_1d_1d(state);
}