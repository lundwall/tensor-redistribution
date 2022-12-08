#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <cstring>

#include "dace_helper.h"
#include "redistribute_by_dimension.hpp"

void calculate_recv_redistribution_blocks(int current_dim, int total_dim, int* pcoords, int* subsize, int* from, int* to, int* lo_send, redistribution_info* state, block_description* recv_block_descriptions, int* recv_from_ranks, int myrank)
{
	if (total_dim == current_dim)
	{
		int cart_rank = get_cart_rank(state->send_dimension, state->send_pgrid, pcoords);
		if (cart_rank != myrank)
		{
			recv_block_descriptions[state->recv_count].from = new int[state->recv_dimension]();
			recv_block_descriptions[state->recv_count].to = new int[state->recv_dimension]();
			std::memcpy(recv_block_descriptions[state->recv_count].from, from, sizeof(int)*total_dim);
			std::memcpy(recv_block_descriptions[state->recv_count].to, to, sizeof(int)*total_dim);

			// These are only used if MODE == datatype
			MPI_Type_create_subarray(total_dim, state->recv_array_dimension, subsize, from, MPI_ORDER_C, MPI_INT, &state->recv_types[state->recv_count]);
	        MPI_Type_commit(&state->recv_types[state->recv_count]);
	        recv_from_ranks[state->recv_count] = cart_rank;
	        state->recv_count++;
		}
		else
		{
			std::memcpy(state->self_src+state->copy_count*state->send_dimension, lo_send, state->send_dimension*sizeof(int));
			std::memcpy(state->self_dst+state->copy_count*state->recv_dimension, from, state->recv_dimension*sizeof(int));
			std::memcpy(state->self_size+state->copy_count*state->send_dimension, subsize, state->send_dimension*sizeof(int));
			state->copy_count++;
		}
		return;
	}

	int xi = state->recv_coords[current_dim] * state->recv_array_dimension[current_dim] / state->send_array_dimension[current_dim];  // (__state->__pgrid_1_coords[0] * int(P*m)) / int(m);
    int lambda = state->recv_coords[current_dim] * state->recv_array_dimension[current_dim] % state->send_array_dimension[current_dim]; // __state->__pgrid_1_coords[0] * int(P*m) % int(m);
    int kappa = int_ceil(int(state->recv_array_dimension[current_dim]) + lambda, int(state->send_array_dimension[current_dim])); // int_ceil(int(P*m) + lambda[0], int(m));

    int rem = state->recv_array_dimension[current_dim]; // By
	for (auto idx = 0; idx < kappa; ++idx)
	{
		pcoords[current_dim] = xi + idx;
        int lo = (idx == 0 ? lambda : 0);
        int uo = std::min(int(state->send_array_dimension[current_dim]), lo + rem);
        subsize[current_dim] = uo - lo;
        from[current_dim] = state->recv_array_dimension[current_dim] - rem;
        to[current_dim] = (state->recv_array_dimension[current_dim] - rem) + (uo - lo);
        lo_send[current_dim] = lo;
       	rem -= uo - lo;
       	calculate_recv_redistribution_blocks(current_dim + 1, total_dim, pcoords, subsize, from, to, lo_send, state, recv_block_descriptions, recv_from_ranks, myrank);
	}
}

void calculate_send_redistribution_blocks(int current_dim, int total_dim, int* pcoords, int* subsize, int* from, int* to, int* lo_send, redistribution_info* state, block_description* send_block_descriptions, int* send_to_ranks, int myrank)
{
	if (total_dim == current_dim)
	{
		int cart_rank = get_cart_rank(state->recv_dimension, state->recv_pgrid, pcoords);
		if (cart_rank != myrank)
		{
			send_block_descriptions[state->send_count].from = new int[state->send_dimension]();
			send_block_descriptions[state->send_count].to = new int[state->send_dimension]();
			std::memcpy(send_block_descriptions[state->send_count].from, from, sizeof(int)*total_dim);
			std::memcpy(send_block_descriptions[state->send_count].to, to, sizeof(int)*total_dim);

			// These are only used if MODE == datatype
			MPI_Type_create_subarray(total_dim, state->send_array_dimension, subsize, from, MPI_ORDER_C, MPI_INT, &state->send_types[state->send_count]);
	        MPI_Type_commit(&state->send_types[state->send_count]);
	        send_to_ranks[state->send_count] = cart_rank;
	        state->send_count++;
		}
		return;
	}

	int lp = std::max(0, (state->send_coords[current_dim] * state->send_array_dimension[current_dim]) / state->recv_array_dimension[current_dim]);
	int up = std::min(state->recv_pgrid[current_dim], int_ceil((state->send_coords[current_dim] + 1) * state->send_array_dimension[current_dim], state->recv_array_dimension[current_dim]));

	for (auto idx = lp; idx < up; ++idx)
	{
        int actual_idx0 = current_dim;

        int xi = (idx * state->recv_array_dimension[current_dim]) / state->send_array_dimension[current_dim];
        int lambda = idx * state->recv_array_dimension[current_dim] % state->send_array_dimension[current_dim];
        int kappa = int_ceil(state->recv_array_dimension[current_dim] + lambda, state->send_array_dimension[current_dim]);
        int idx_dst = state->send_coords[current_dim] - xi;

        if (idx_dst < 0 || idx_dst >= kappa) continue;
        int lo = (idx_dst == 0 ? lambda : 0);
        int uo = (idx_dst == kappa - 1 ? state->recv_array_dimension[current_dim] + lambda - idx_dst * state->send_array_dimension[current_dim] : state->send_array_dimension[current_dim]);
        subsize[current_dim] = uo - lo;
        from[current_dim] = lo;
        pcoords[current_dim] = idx;
        to[current_dim] = from[current_dim] + subsize[current_dim];
        calculate_send_redistribution_blocks(current_dim+1, total_dim, pcoords, subsize, from, to, lo_send, state, send_block_descriptions, send_to_ranks, myrank);
	}
}

void fill_redistribution_information(redistribution_info* state, int myrank)
{
    int* pcoords = new int[state->recv_dimension]();
    int* subsize = new int[state->recv_dimension]();
    int* from = new int[state->recv_dimension]();
    int* to = new int[state->recv_dimension]();
    int* lo_send = new int[state->send_dimension]();

    calculate_send_redistribution_blocks(0, state->recv_dimension, pcoords, subsize, from, to, lo_send, state, state->send_block_descriptions, state->send_to_ranks, myrank);
    calculate_recv_redistribution_blocks(0, state->send_dimension, pcoords, subsize, from, to, lo_send, state, state->recv_block_descriptions, state->recv_from_ranks, myrank);
    
    return;
}

void redistribute(redistribution_info* state, int* A, int* A_shape_in, int* B, int* B_shape_in, std::string MODE)
{
	redistribute_by_dimension(state, A, A_shape_in, B, B_shape_in, MODE);
}