/*
 * This code is for 2022 DPHPC course project, for topic "tensor redistribution".
 * As it's related to Dace project, some code is reused from / modified by Dace framework.
 */

#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <cstring>

#include "dace_helper.h"
#include "redistribute_by_dimension.hpp"

void parse_array_dimension(std::vector<int>& dimension_vector, std::string input_value)
{
	std::size_t left_param = input_value.find("[");
	std::size_t right_param = input_value.find("]");
	if (left_param == std::string::npos || right_param == std::string::npos || right_param < left_param+2)
	{
		std::cout << "Use format [dim1,dim2,dim3,...] to specify array dimensions.";
		return;
	}
	std::string dimension_str = input_value.substr(left_param+1,right_param-left_param-1);

	size_t delimiter_pos;
	while ((delimiter_pos = dimension_str.find(",")) != std::string::npos)
	{
		std::string new_dimension = dimension_str.substr(0, delimiter_pos);
		try
		{
			int dimension = std::stoi(new_dimension);
			dimension_vector.push_back(dimension);
		}
		catch(...)
		{
			std::cout << "Encounter error when converting dimensions to integer";
		}
		dimension_str.erase(0, delimiter_pos+1);
	}
}

void check_parameter_valid(std::vector<int>& pgrid_send, std::vector<int>& pgrid_recv, std::vector<int>& total_array_size, int total_processor)
{
	if (pgrid_send.empty() || pgrid_recv.empty() || total_array_size.empty())
	{
		std::cout << "Required parameter value missed. Please check all parameters are assigned in cfg file.";
		std::exit(EXIT_FAILURE);
	}

	if (pgrid_send.size() != total_array_size.size() || pgrid_send.size() != pgrid_recv.size())
	{
		std::cout << "processor grid dimension and  array size dimension should be same.";
		std::exit(EXIT_FAILURE);
	}

	int send_size = 1, recv_size = 1;
	for (int i = 0; i < pgrid_send.size(); i++) 
	{
		send_size *= pgrid_send[i];
		recv_size *= pgrid_recv[i];
	}

	if (send_size != total_processor || recv_size != total_processor)
	{
		std::cout << "the product of processor grid should be same as mpi process number." << std::endl;
		std::exit(EXIT_FAILURE);
	}

}

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

	int send_array_dimension = state->total_array_size[current_dim] / state->send_pgrid[current_dim];  // Bx
	int recv_array_dimension = state->total_array_size[current_dim] / state->recv_pgrid[current_dim];  // By
	int xi = state->recv_coords[current_dim] * recv_array_dimension / send_array_dimension;  // (__state->__pgrid_1_coords[0] * int(P*m)) / int(m);
    int lambda = state->recv_coords[current_dim] * recv_array_dimension % send_array_dimension; // __state->__pgrid_1_coords[0] * int(P*m) % int(m);
    int kappa = int_ceil(int(recv_array_dimension) + lambda, int(send_array_dimension)); // int_ceil(int(P*m) + lambda[0], int(m));

    int rem = recv_array_dimension; // By
	for (auto idx = 0; idx < kappa; ++idx)
	{
		pcoords[current_dim] = xi + idx;
        int lo = (idx == 0 ? lambda : 0);
        int uo = std::min(int(send_array_dimension), lo + rem);
        subsize[current_dim] = uo - lo;
        from[current_dim] = recv_array_dimension - rem;
        to[current_dim] = (recv_array_dimension - rem) + (uo - lo);
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
			send_to_ranks[state->send_count] = cart_rank;
			state->send_count++;
		}
		return;
	}

	int send_array_dimension = state->total_array_size[current_dim] / state->send_pgrid[current_dim];  // Bx
	int recv_array_dimension = state->total_array_size[current_dim] / state->recv_pgrid[current_dim];  // By
	state->send_array_dimension[current_dim] = send_array_dimension;
	state->recv_array_dimension[current_dim] = recv_array_dimension;
	state->send_element_count *= send_array_dimension;
	state->recv_element_count *= recv_array_dimension;
	int lp = std::max(0, (state->send_coords[current_dim] * send_array_dimension) / recv_array_dimension);
	int up = std::min(state->recv_pgrid[current_dim], int_ceil((state->send_coords[current_dim] + 1) * send_array_dimension, recv_array_dimension));

	for (auto idx = lp; idx < up; ++idx)
	{
        int actual_idx0 = current_dim;

        int xi = (idx * recv_array_dimension) / send_array_dimension;
        int lambda = idx * recv_array_dimension % send_array_dimension;
        int kappa = int_ceil(recv_array_dimension + lambda, send_array_dimension);
        int idx_dst = state->send_coords[current_dim] - xi;

        if (idx_dst < 0 || idx_dst >= kappa) continue;
        int lo = (idx_dst == 0 ? lambda : 0);
        int uo = (idx_dst == kappa - 1 ? recv_array_dimension + lambda - idx_dst * send_array_dimension : send_array_dimension);
        subsize[current_dim] = uo - lo;
        from[current_dim] = lo;
        pcoords[current_dim] = idx;
        to[current_dim] = from[current_dim] + subsize[current_dim];
        calculate_send_redistribution_blocks(current_dim+1, total_dim, pcoords, subsize, from, to, lo_send, state, send_block_descriptions, send_to_ranks, myrank);
	}
}

void redistribute(redistribution_info* state, int* A, int* A_shape_in, int* B, int* B_shape_in)
{
	redistribute_by_dimension(state, A, A_shape_in, B, B_shape_in);
}


void fill_redistribution_information(redistribution_info* state, int max_sends, int max_recvs, int myrank)
{
	int* pcoords = new int[state->recv_dimension]();
	int* subsize = new int[state->recv_dimension]();
	int* from = new int[state->recv_dimension]();
	int* to = new int[state->recv_dimension]();
	int* lo_send = new int[state->send_dimension]();
	state->recv_block_descriptions = new block_description[max_recvs];
	state->recv_from_ranks = new int[max_recvs];
	state->send_block_descriptions = new block_description[max_sends];
	state->send_to_ranks = new int[max_sends];
	state->send_array_dimension = new int[state->send_dimension]();
	state->recv_array_dimension = new int[state->recv_dimension]();
	state->send_element_count = 1;
	state->recv_element_count = 1;

	calculate_recv_redistribution_blocks(0, state->send_dimension, pcoords, subsize, from, to, lo_send, state, state->recv_block_descriptions, state->recv_from_ranks, myrank);
	calculate_send_redistribution_blocks(0, state->recv_dimension, pcoords, subsize, from, to, lo_send, state, state->send_block_descriptions, state->send_to_ranks, myrank);
	
	return;
}

void print_debug_message(redistribution_info* state, int myrank)
{
	for (int i = 0; i < state->recv_count; i++)
	{
		std::string message = "I am rank " + std::to_string(myrank) + " and I receive from " + std::to_string(state->recv_from_ranks[i]);
		message += " in new layout it's from ";
		for (int j = 0; j < state->recv_dimension; j++)
		{
			message += " " + std::to_string(state->recv_block_descriptions[i].from[j]);
		}

		message += " to ";
		for (int j = 0; j < state->recv_dimension; j++)
		{
			message += " " + std::to_string(state->recv_block_descriptions[i].to[j]);
		}
		std::cout << message << std::endl;
	}

	for (int i = 0; i < state->copy_count; i++)
	{
		std::string message = "I am rank " + std::to_string(myrank) + " and I self-copy, src is ";
		for (int j = 0; j < state->send_dimension; j++)
		{
			message += " " + std::to_string(state->self_src[i*state->copy_count+j]);
		}
		message += " dst is ";
		for (int j = 0; j < state->recv_dimension; j++)
		{
			message += " " + std::to_string(state->self_dst[i*state->copy_count+j]);
		}
		message += " size is ";
		for (int j = 0; j < state->send_dimension; j++)
		{
			message += " " + std::to_string(state->self_size[i*state->copy_count+j]);
		}
		std::cout << message << std::endl;
	}

	for (int i = 0; i < state->send_count; i++)
	{
		std::string message = "I am rank " + std::to_string(myrank) + " and I send to " + std::to_string(state->send_to_ranks[i]);
		message += " in old layout it's from ";
		for (int j = 0; j < state->send_dimension; j++)
		{
			message += " " + std::to_string(state->send_block_descriptions[i].from[j]);
		}

		message += " to ";
		for (int j = 0; j < state->recv_dimension; j++)
		{
			message += " " + std::to_string(state->send_block_descriptions[i].to[j]);
		}
		std::cout << message << std::endl;
	}
}

int main(int argc, char** argv)
{
	int size, rank, received_threads;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &received_threads);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("my rank is %d\n", rank);

	std::vector<int> pgrid_send;
	std::vector<int> pgrid_recv;
	std::vector<int> total_array_size;
	redistribution_info* state = new redistribution_info; 

	std::ifstream cfg_input;
	cfg_input.open("redistribute.cfg");
	for (std::string line; std::getline(cfg_input, line);)
	{
		std::size_t found = line.find("=");
		if (found == std::string::npos)
			continue;
		std::string key = line.substr(0, found);
		std::string value = line.substr(found+1, line.length()-key.length()-1);

		if (key == "processor_grid_send")
		{
			parse_array_dimension(pgrid_send, value);
		}
		else if (key == "processor_grid_recv")
		{
			parse_array_dimension(pgrid_recv, value);
		}
		else if (key == "total_array_size")
		{
			parse_array_dimension(total_array_size, value);
		}
		else
		{
			std::cout << "unknown parameter " << key << " detected and ignore" << std::endl;
		}
	}
	cfg_input.close();

	state->send_pgrid = pgrid_send.data();
	state->recv_pgrid = pgrid_recv.data();
	state->total_array_size = total_array_size.data();
	state->send_coords = new int[pgrid_send.size()]();
	state->recv_coords = new int[pgrid_recv.size()]();
	int* pgrid_period_send =  new int[pgrid_send.size()]();
	int* pgrid_period_recv =  new int[pgrid_recv.size()]();
	state->send_dimension = pgrid_send.size();
	state->recv_dimension = pgrid_recv.size();
	state->send_count = 0;
	state->recv_count = 0;
	state->copy_count = 0;

	int max_sends = 1;
	int max_recvs = 1;
	for (int i = 0; i < state->send_dimension; i++)
	{
		int Bx = total_array_size[i] / pgrid_send[i];
		int By = total_array_size[i] / pgrid_recv[i];
		max_sends *= int_ceil(Bx + By - 1, By);
    	max_recvs *= int_ceil(By - 1, Bx) + 1;
    }
    state->self_src = new int[max_sends * state->send_dimension];
    state->self_dst = new int[max_recvs * state->recv_dimension];
    state->self_size = new int[max_sends * state->send_dimension];
    std::cout << "max_sends = " << max_sends << " max_recvs = " << max_recvs << std::endl;

	check_parameter_valid(pgrid_send, pgrid_recv, total_array_size, size);

	MPI_Cart_create(MPI_COMM_WORLD, state->send_dimension, state->send_pgrid, pgrid_period_send, 0, &state->send_comm);
	MPI_Cart_create(MPI_COMM_WORLD, state->recv_dimension, state->recv_pgrid, pgrid_period_recv, 0, &state->recv_comm);
	MPI_Barrier(MPI_COMM_WORLD);

	if (state->send_comm == MPI_COMM_NULL)
	{
		std::cout << "MPI_Cart_create error" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	MPI_Comm_group(state->send_comm, &state->send_group);
    MPI_Comm_rank(state->send_comm, &state->send_rank);
    MPI_Comm_size(state->send_comm, &state->send_size);
    MPI_Cart_coords(state->send_comm, state->send_rank, state->send_dimension, state->send_coords);

	if (state->recv_comm == MPI_COMM_NULL)
	{
		std::cout << "MPI_Cart_create error" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	MPI_Comm_group(state->recv_comm, &state->recv_group);
    MPI_Comm_rank(state->recv_comm, &state->recv_rank);
    MPI_Comm_size(state->recv_comm, &state->recv_size);
    MPI_Cart_coords(state->recv_comm, state->recv_rank, state->recv_dimension, state->recv_coords);


    fill_redistribution_information(state, max_sends, max_recvs, rank);

    print_debug_message(state, rank);
    int* A = new int[state->send_element_count];
    int* B = new int[state->recv_element_count];

    redistribute(state, A, state->send_array_dimension, B, state->recv_array_dimension);
    MPI_Finalize();
	return 0;
}
