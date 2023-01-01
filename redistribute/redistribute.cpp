/*
 * This code is for 2022 DPHPC course project, for topic "tensor redistribution".
 * As it's related to Dace project, some code is reused from / modified by Dace framework.
 */

#include "redistribute.hpp"

#define DEBUG false

#define RUN 1000
#define WARMUP 10
#define SYNC 10

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
    read_cfg_file(pgrid_send, pgrid_recv, total_array_size, "redistribute.cfg");
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

	for (int threads = 4; threads <= 4; ++threads)
	{
        omp_set_num_threads(threads);

		for (int num_chunks = 1; num_chunks <= 10; ++num_chunks)
		{
			const char* modes[3];
			modes[0] = "manual";
			modes[1] = threads > 1 || num_chunks > 1 ? NULL : "datatype";
			modes[2] = NULL;
			for (int i=0; modes[i] != NULL; i++)
			{
				std::string MODE = modes[i];
				// Clear received data before redistributing again
				for (int i = 0; i < state->recv_element_count; i++)
				{
					B[i] = 0;
				}

				// LSB timing setup
				std::string name = "5d_redistribute_" + MODE + "_t" + std::to_string(omp_get_max_threads());
				LSB_Init(name.c_str(), 0);
				LSB_Set_Rparam_int("rank", rank);
				LSB_Set_Rparam_string("mode", MODE.c_str());
				LSB_Set_Rparam_string("type", "sync");
				LSB_Set_Rparam_int("threads", omp_get_max_threads());
				LSB_Set_Rparam_int("chunks", num_chunks);
				LSB_Set_Rparam_double("err", 0); // meaningless here
				double win;
				double* recorded_values = new double[1000];
				bool finished = false;
				LSB_Rec_disable();

				for (auto run_idx = 0; run_idx < WARMUP + SYNC + RUN && !finished; ++run_idx)
				{
					configure_LSB_and_sync(run_idx, WARMUP, SYNC, &win);
					LSB_Res();

					redistribute(state, A, state->send_array_dimension, B, state->recv_array_dimension, MODE, num_chunks);

					int num_recorded_values = run_idx - WARMUP - SYNC + 1;
					LSB_Rec(std::max(num_recorded_values, 0));
					if (num_recorded_values >= 1)
					{
						aggregate_CIs(num_recorded_values, recorded_values, size, &finished);
					}
				}

				LSB_Finalize();
			}
		}
	}

    MPI_Finalize();
	return 0;
}
