/*
 * This code is for 2022 DPHPC course project, for topic "tensor redistribution".
 * As it's related to Dace project, some code is reused from / modified by Dace framework.
 */

#include "redistribute.hpp"
#define DEBUG false

int main(int argc, char** argv)
{
	int size, rank, received_threads;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &received_threads);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::vector<int> pgrid_send;
	std::vector<int> pgrid_recv;
	std::vector<int> total_array_size;
	redistribution_info* state = new redistribution_info; 
    read_cfg_file(pgrid_send, pgrid_recv, total_array_size);
    fill_state_information(pgrid_send, pgrid_recv, total_array_size, state);

	check_parameter_valid(pgrid_send, pgrid_recv, total_array_size, size);

	int* pgrid_period_send =  new int[pgrid_send.size()]();
	int* pgrid_period_recv =  new int[pgrid_recv.size()]();
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

    fill_redistribution_information(state, rank);

    if (DEBUG)
    	print_debug_message(state, rank);

    int* A = new int[state->send_element_count];
    int* B = new int[state->recv_element_count];

    for (int i = 0; i < state->send_element_count; i++)
    {
    	A[i] = i;
    }

    for (int i = 0; i < state->recv_element_count; i++)
    {
    	B[i] = 0;
   	}

	std::string MODE = argv[1];
    redistribute(state, A, state->send_array_dimension, B, state->recv_array_dimension, MODE);

    MPI_Finalize();
	return 0;
}
