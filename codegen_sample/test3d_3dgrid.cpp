/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"
#include "mpi.h"

struct matrix_1d_1d_t {
    MPI_Comm __pgrid_0_comm;
    MPI_Group __pgrid_0_group;
    int __pgrid_0_coords[3];
    int __pgrid_0_dims[3];
    int __pgrid_0_rank;
    int __pgrid_0_size;
    bool __pgrid_0_valid;
    MPI_Comm __pgrid_1_comm;
    MPI_Group __pgrid_1_group;
    int __pgrid_1_coords[3];
    int __pgrid_1_dims[3];
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

        #pragma omp parallel sections
        {
            #pragma omp section
            {
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
                            int __inp_s0 = __state->__rdistrarray_0_self_src[__idx * 3 + 0];
                            int __inp_s1 = __state->__rdistrarray_0_self_src[__idx * 3 + 1];
                            int __inp_s2 = __state->__rdistrarray_0_self_src[__idx * 3 + 2];

                            int __out_s0 = __state->__rdistrarray_0_self_dst[__idx * 3 + 0];
                            int __out_s1 = __state->__rdistrarray_0_self_dst[__idx * 3 + 1];
                            int __out_s2 = __state->__rdistrarray_0_self_dst[__idx * 3 + 2];

                            dace::CopyNDDynamic<int, 1, false, 3>::Dynamic::Copy(
                            _inp_buffer + ((((((P * P) * __inp_s0) * m) * m) + ((P * __inp_s1) * m)) + __inp_s2), _out_buffer + (((((P * __out_s0) * m) * m) + ((P * __out_s1) * m)) + __out_s2), __state->__rdistrarray_0_self_size[__idx * 3 + 0], P**2*m**2, P*m**2, __state->__rdistrarray_0_self_size[__idx * 3 + 1], P*m, P*m, __state->__rdistrarray_0_self_size[__idx * 3 + 2], 1, 1
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
            } // End omp section
            #pragma omp section
            {
                {
                    int __out;

                    ///////////////////
                    // Tasklet code (__pgrid_0)
                    ///////////////////

                    a_grid = __out;
                }
            } // End omp section
            #pragma omp section
            {
                {
                    int __out;

                    ///////////////////
                    // Tasklet code (__pgrid_1)
                    ///////////////////

                    b_grid = __out;
                }
            } // End omp section
            #pragma omp section
            {
                {
                    int __out;

                    ///////////////////
                    // Tasklet code (__subarray_0)
                    ///////////////////

                    a_arr = __out;
                }
            } // End omp section
            #pragma omp section
            {
                {
                    int __out;

                    ///////////////////
                    // Tasklet code (__subarray_1)
                    ///////////////////

                    b_arr = __out;
                }
            } // End omp section
            #pragma omp section
            {
                {
                    int __out;

                    ///////////////////
                    // Tasklet code (__rdistrarray_0)
                    ///////////////////

                    rdistr = __out;
                }
            } // End omp section
        } // End omp sections

    }
}

DACE_EXPORTED void __program_matrix_1d_1d(matrix_1d_1d_t *__state, int * __restrict__ A, int * __restrict__ __return, int P, int m)
{
    __program_matrix_1d_1d_internal(__state, A, __return, P, m);
}

DACE_EXPORTED matrix_1d_1d_t *__dace_init_matrix_1d_1d(int P, int m)
{
    int __result = 0;
    matrix_1d_1d_t *__state = new matrix_1d_1d_t;


    {  // Environment: MPI
        int t; MPI_Initialized(&t);  if (!t) MPI_Init(NULL, NULL);
    }
    __state->__pgrid_0_dims[0] = P;
    __state->__pgrid_0_dims[1] = 1;
    __state->__pgrid_0_dims[2] = 1;

    int __pgrid_0_periods[3] = {0};
    MPI_Cart_create(MPI_COMM_WORLD, 3, __state->__pgrid_0_dims, __pgrid_0_periods, 0, &__state->__pgrid_0_comm);
    if (__state->__pgrid_0_comm != MPI_COMM_NULL) {
        MPI_Comm_group(__state->__pgrid_0_comm, &__state->__pgrid_0_group);
        MPI_Comm_rank(__state->__pgrid_0_comm, &__state->__pgrid_0_rank);
        MPI_Comm_size(__state->__pgrid_0_comm, &__state->__pgrid_0_size);
        MPI_Cart_coords(__state->__pgrid_0_comm, __state->__pgrid_0_rank, 3, __state->__pgrid_0_coords);
        __state->__pgrid_0_valid = true;
    } else {
        __state->__pgrid_0_group = MPI_GROUP_NULL;
        __state->__pgrid_0_rank = MPI_PROC_NULL;
        __state->__pgrid_0_size = 0;
        __state->__pgrid_0_valid = false;
    }
    __state->__pgrid_1_dims[0] = 1;
    __state->__pgrid_1_dims[1] = P;
    __state->__pgrid_1_dims[2] = 1;

    int __pgrid_1_periods[3] = {0};
    MPI_Cart_create(MPI_COMM_WORLD, 3, __state->__pgrid_1_dims, __pgrid_1_periods, 0, &__state->__pgrid_1_comm);
    if (__state->__pgrid_1_comm != MPI_COMM_NULL) {
        MPI_Comm_group(__state->__pgrid_1_comm, &__state->__pgrid_1_group);
        MPI_Comm_rank(__state->__pgrid_1_comm, &__state->__pgrid_1_rank);
        MPI_Comm_size(__state->__pgrid_1_comm, &__state->__pgrid_1_size);
        MPI_Cart_coords(__state->__pgrid_1_comm, __state->__pgrid_1_rank, 3, __state->__pgrid_1_coords);
        __state->__pgrid_1_valid = true;
    } else {
        __state->__pgrid_1_group = MPI_GROUP_NULL;
        __state->__pgrid_1_rank = MPI_PROC_NULL;
        __state->__pgrid_1_size = 0;
        __state->__pgrid_1_valid = false;
    }

    if (__state->__pgrid_0_valid) {
        int sizes[3] = {P*m, P*m, P*m};
        int subsizes[3] = {m, P*m, P*m};
        int corr[3] = {0, 1, 2};

        int basic_stride = subsizes[3 - 1];

        int process_strides[3];
        int block_strides[3];
        int data_strides[3];

        process_strides[3 - 1] = 1;
        block_strides[3 - 1] = subsizes[3 - 1];
        data_strides[3 - 1] = 1;

        for (auto i = 3 - 2; i >= 0; --i) {
            block_strides[i] = block_strides[i+1] * subsizes[i];
            process_strides[i] = process_strides[i+1] * __state->__pgrid_0_dims[corr[i+1]];
            data_strides[i] = block_strides[i] * process_strides[i] / basic_stride;
        }

        MPI_Datatype type;
        int origin[3] = {0,0,0};
        MPI_Type_create_subarray(3, sizes, subsizes, origin, MPI_ORDER_C, MPI_INT, &type);
        MPI_Type_create_resized(type, 0, basic_stride*sizeof(int), &__state->__subarray_0);
        MPI_Type_commit(&__state->__subarray_0);

        __state->__subarray_0_counts = new int[__state->__pgrid_0_size];
        __state->__subarray_0_displs = new int[__state->__pgrid_0_size];
        int block_id[3] = {0};
        int displ = 0;
        for (auto i = 0; i < __state->__pgrid_0_size; ++i) {
            __state->__subarray_0_counts[i] = 1;
            __state->__subarray_0_displs[i] = displ;
            int idx = 3 - 1;
            while (idx >= 0 && block_id[idx] + 1 >= __state->__pgrid_0_dims[corr[idx]]) {
                block_id[idx] = 0;
                displ -= data_strides[idx] * (__state->__pgrid_0_dims[corr[idx]] - 1);
                idx--;
            }
            if (idx >= 0) {
                block_id[idx] += 1;
                displ += data_strides[idx];
            } else {
                assert(i == __state->__pgrid_0_size - 1);
            }
        }
    }

    if (__state->__pgrid_1_valid) {
        int sizes[3] = {P*m, P*m, P*m};
        int subsizes[3] = {P*m, m, P*m};
        int corr[3] = {0, 1, 2};

        int basic_stride = subsizes[3 - 1];

        int process_strides[3];
        int block_strides[3];
        int data_strides[3];

        process_strides[3 - 1] = 1;
        block_strides[3 - 1] = subsizes[3 - 1];
        data_strides[3 - 1] = 1;

        for (auto i = 3 - 2; i >= 0; --i) {
            block_strides[i] = block_strides[i+1] * subsizes[i];
            process_strides[i] = process_strides[i+1] * __state->__pgrid_1_dims[corr[i+1]];
            data_strides[i] = block_strides[i] * process_strides[i] / basic_stride;
        }

        MPI_Datatype type;
        int origin[3] = {0,0,0};
        MPI_Type_create_subarray(3, sizes, subsizes, origin, MPI_ORDER_C, MPI_INT, &type);
        MPI_Type_create_resized(type, 0, basic_stride*sizeof(int), &__state->__subarray_1);
        MPI_Type_commit(&__state->__subarray_1);

        __state->__subarray_1_counts = new int[__state->__pgrid_1_size];
        __state->__subarray_1_displs = new int[__state->__pgrid_1_size];
        int block_id[3] = {0};
        int displ = 0;
        for (auto i = 0; i < __state->__pgrid_1_size; ++i) {
            __state->__subarray_1_counts[i] = 1;
            __state->__subarray_1_displs[i] = displ;
            int idx = 3 - 1;
            while (idx >= 0 && block_id[idx] + 1 >= __state->__pgrid_1_dims[corr[idx]]) {
                block_id[idx] = 0;
                displ -= data_strides[idx] * (__state->__pgrid_1_dims[corr[idx]] - 1);
                idx--;
            }
            if (idx >= 0) {
                block_id[idx] += 1;
                displ += data_strides[idx];
            } else {
                assert(i == __state->__pgrid_1_size - 1);
            }
        }
    }
    {
        __state->__rdistrarray_0_sends = 0;
        __state->__rdistrarray_0_recvs = 0;
        __state->__rdistrarray_0_self_copies = 0;
        int max_sends = 1;
        int max_recvs = 1;

        int kappa[3];
        int lambda[3];
        int xi[3];
        int pcoords[3];

        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        max_sends *= int_ceil(int(m) + int(P*m) - 1, int(P*m));
        max_recvs *= int_ceil(int(P*m) - 1, int(m)) + 1;

        max_sends *= int_ceil(int(P*m) + int(m) - 1, int(m));
        max_recvs *= int_ceil(int(m) - 1, int(P*m)) + 1;

        max_sends *= int_ceil(int(P*m) + int(P*m) - 1, int(P*m));
        max_recvs *= int_ceil(int(P*m) - 1, int(P*m)) + 1;

        __state->__rdistrarray_0_send_types = new MPI_Datatype[max_sends];
        __state->__rdistrarray_0_dst_ranks = new int[max_sends];
        __state->__rdistrarray_0_recv_types = new MPI_Datatype[max_recvs];
        __state->__rdistrarray_0_src_ranks = new int[max_recvs];
        __state->__rdistrarray_0_self_src = new int[max_sends * 3];
        __state->__rdistrarray_0_self_dst = new int[max_sends * 3];
        __state->__rdistrarray_0_self_size = new int[max_sends * 3];

        if (__state->__pgrid_1_valid) {

            int sizes[3] = {P*m, m, P*m};
            int subsizes[3];
            int origin[3];

            xi[0] = (__state->__pgrid_1_coords[0] * int(P*m)) / int(m);
            lambda[0] = __state->__pgrid_1_coords[0] * int(P*m) % int(m);
            kappa[0] = int_ceil(int(P*m) + lambda[0], int(m));

            xi[1] = (__state->__pgrid_1_coords[1] * int(m)) / int(P*m);
            lambda[1] = __state->__pgrid_1_coords[1] * int(m) % int(P*m);
            kappa[1] = int_ceil(int(m) + lambda[1], int(P*m));

            xi[2] = (__state->__pgrid_1_coords[2] * int(P*m)) / int(P*m);
            lambda[2] = __state->__pgrid_1_coords[2] * int(P*m) % int(P*m);
            kappa[2] = int_ceil(int(P*m) + lambda[2], int(P*m));

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

                    int rem2 = P*m;
                    for (auto idx2 = 0; idx2 < kappa[2]; ++idx2) {
                        int actual_idx2 = 2;
                        pcoords[actual_idx2] = xi[2] + idx2;
                        int lo2 = (idx2 == 0 ? lambda[2] : 0);
                        int uo2 = std::min(int(P*m), lo2 + rem2);
                        subsizes[2] = uo2 - lo2;
                        origin[2] = P*m - rem2;
                        rem2 -= uo2 - lo2;
                        int cart_rank = dace::comm::cart_rank(3, __state->__pgrid_0_dims, pcoords);
                        if (myrank == cart_rank) { // self-copy
                            __state->__rdistrarray_0_self_src[__state->__rdistrarray_0_self_copies * 3 + 0] = lo0;
                            __state->__rdistrarray_0_self_dst[__state->__rdistrarray_0_self_copies * 3 + 0] = origin[0];
                            __state->__rdistrarray_0_self_size[__state->__rdistrarray_0_self_copies * 3 + 0] = subsizes[0];

                            __state->__rdistrarray_0_self_src[__state->__rdistrarray_0_self_copies * 3 + 1] = lo1;
                            __state->__rdistrarray_0_self_dst[__state->__rdistrarray_0_self_copies * 3 + 1] = origin[1];
                            __state->__rdistrarray_0_self_size[__state->__rdistrarray_0_self_copies * 3 + 1] = subsizes[1];

                            __state->__rdistrarray_0_self_src[__state->__rdistrarray_0_self_copies * 3 + 2] = lo2;
                            __state->__rdistrarray_0_self_dst[__state->__rdistrarray_0_self_copies * 3 + 2] = origin[2];
                            __state->__rdistrarray_0_self_size[__state->__rdistrarray_0_self_copies * 3 + 2] = subsizes[2];

                            __state->__rdistrarray_0_self_copies++;
                            // printf("(__subarray_0 -> __subarray_1) I am rank %d and I self-copy {I receive from %d%d (%d - %d) in (%d, %d) size (%d, %d)} \n", myrank, pcoords[0], pcoords[1], cart_rank, cart_rank, origin[0], origin[1], subsizes[0], subsizes[1]);
                        } else {
                            MPI_Type_create_subarray(3,  sizes, subsizes, origin, MPI_ORDER_C, MPI_INT, &__state->__rdistrarray_0_recv_types[__state->__rdistrarray_0_recvs]);
                            MPI_Type_commit(&__state->__rdistrarray_0_recv_types[__state->__rdistrarray_0_recvs]);
                            __state->__rdistrarray_0_src_ranks[__state->__rdistrarray_0_recvs] = cart_rank;
                            // printf("(__subarray_0 -> __subarray_1) I am rank %d and I receive from %d%d (%d - %d) in (%d, %d) size (%d, %d) \n", myrank, pcoords[0], pcoords[1], cart_rank, __state->__rdistrarray_0_src_ranks[__state->__rdistrarray_0_recvs], origin[0], origin[1], subsizes[0], subsizes[1]);
                            __state->__rdistrarray_0_recvs++;
                        }
        }}}}
        if (__state->__pgrid_0_valid) {

            int sizes[3] = {m, P*m, P*m};
            int subsizes[3];
            int origin[3];

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


                    // int_ceil(x, y) := (x + y - 1) / y
                    // int_ceil(pcoord * sa - sb + 1, sb) = (pcoord * sa) / sb
                    int lp2 = std::max(0, (__state->__pgrid_0_coords[2] * int(P*m)) / int(P*m)); // int_ceil(x, y) := (x + y - 1) / y
                    int up2 = std::min(__state->__pgrid_1_dims[2], int_ceil((__state->__pgrid_0_coords[2] + 1) * int(P*m), int(P*m)));
                    // printf("I am rank %d and I have 2-th bounds [%d, %d)\n", myrank, lp2, up2);
                    for (auto idx2 = lp2; idx2 < up2; ++idx2) {
                        int actual_idx2 = 2;

                        xi[2] = (idx2 * int(P*m)) / int(P*m);
                        lambda[2] = idx2 * int(P*m) % int(P*m);
                        kappa[2] = int_ceil(int(P*m) + lambda[2], int(P*m));
                        int idx2_dst = __state->__pgrid_0_coords[2] - xi[2];

                        if (idx2_dst < 0 || idx2_dst >= kappa[2]) continue;
                        int lo2 = (idx2_dst == 0 ? lambda[2] : 0);
                        int uo2 = (idx2_dst == kappa[2] - 1 ? int(P*m) + lambda[2] - idx2_dst * int(P*m) : int(P*m));
                        subsizes[2] = uo2 - lo2;
                        origin[2] = lo2;
                        pcoords[actual_idx2] = idx2;

                        int cart_rank = dace::comm::cart_rank(3, __state->__pgrid_1_dims, pcoords);

                        if (myrank != cart_rank) { // not self-copy
                            MPI_Type_create_subarray(3,  sizes, subsizes, origin, MPI_ORDER_C, MPI_INT, &__state->__rdistrarray_0_send_types[__state->__rdistrarray_0_sends]);
                            MPI_Type_commit(&__state->__rdistrarray_0_send_types[__state->__rdistrarray_0_sends]);
                            __state->__rdistrarray_0_dst_ranks[__state->__rdistrarray_0_sends] = cart_rank;
                            // printf("(__subarray_0 -> __subarray_1) I am rank %d and I send to %d%d (%d - %d) from (%d, %d) size (%d, %d)\n", myrank, pcoords[0], pcoords[1], cart_rank, __state->__rdistrarray_0_dst_ranks[__state->__rdistrarray_0_sends], origin[0], origin[1], subsizes[0], subsizes[1]);
                            __state->__rdistrarray_0_sends++;
                        }
    }}}}}

    if (__result) {
        delete __state;
        return nullptr;
    }
    return __state;
}

DACE_EXPORTED void __dace_exit_matrix_1d_1d(matrix_1d_1d_t *__state)
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
        delete[] __state->__subarray_0_counts;
        delete[] __state->__subarray_0_displs;
        MPI_Type_free(&__state->__subarray_0);
    }

    if (__state->__pgrid_1_valid) {
        delete[] __state->__subarray_1_counts;
        delete[] __state->__subarray_1_displs;
        MPI_Type_free(&__state->__subarray_1);
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

