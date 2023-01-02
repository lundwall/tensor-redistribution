#include <mpi.h>
#include <liblsb.h>
#include <math.h>
#include <vector>
#include <fstream>
#include "dace_helper.h"

struct block_description
{
	int* from;
	int* to;
};

struct redistribution_info 
{
	MPI_Comm send_comm;
    MPI_Group send_group;
    int* send_coords;
    int* send_pgrid;
    int send_rank;
    int send_size;
    int send_dimension;
    block_description* send_block_descriptions;
    int* send_to_ranks;
    int* send_array_dimension;
    int send_count;
    int send_element_count;
    MPI_Datatype* send_types;
    MPI_Request* send_req;

    MPI_Comm recv_comm;
    MPI_Group recv_group;
    int* recv_coords;
    int* recv_pgrid;
    int recv_rank;
    int recv_size;
    int recv_dimension;
    block_description* recv_block_descriptions;
    int* recv_from_ranks;
    int* recv_array_dimension;
    int recv_count;
    int recv_element_count;
    MPI_Datatype* recv_types;

    int* total_array_size;

    int* self_src;
    int* self_dst;
    int* self_size;
    int copy_count;
};

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

void read_cfg_file(std::vector<int>& pgrid_send, std::vector<int>& pgrid_recv, std::vector<int>& total_array_size, std::string cfg_filename)
{
    std::ifstream cfg_input;
    cfg_input.open(cfg_filename);
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
}

void fill_state_information(std::vector<int>& pgrid_send, std::vector<int>& pgrid_recv, std::vector<int>& total_array_size, redistribution_info* state)
{
    state->send_pgrid = pgrid_send.data();
    state->recv_pgrid = pgrid_recv.data();
    state->total_array_size = total_array_size.data();
    state->send_coords = new int[pgrid_send.size()]();
    state->recv_coords = new int[pgrid_recv.size()]();

    state->send_dimension = pgrid_send.size();
    state->recv_dimension = pgrid_recv.size();
    state->send_count = 0;
    state->recv_count = 0;
    state->copy_count = 0;
    state->send_array_dimension = new int[state->send_dimension]();
    state->recv_array_dimension = new int[state->recv_dimension]();
    state->send_element_count = 1;
    state->recv_element_count = 1;

    int max_sends = 1;
    int max_recvs = 1;
    for (int i = 0; i < state->send_dimension; i++)
    {
        int Bx = total_array_size[i] / pgrid_send[i];
        int By = total_array_size[i] / pgrid_recv[i];
        max_sends *= int_ceil(Bx + By - 1, By);
        max_recvs *= int_ceil(By - 1, Bx) + 1;
        state->send_array_dimension[i] = state->total_array_size[i] / state->send_pgrid[i];
        state->recv_array_dimension[i] = state->total_array_size[i] / state->recv_pgrid[i];
        state->send_element_count *= state->send_array_dimension[i];
        state->recv_element_count *= state->recv_array_dimension[i];
    }
    state->self_src = new int[max_sends * state->send_dimension];
    state->self_dst = new int[max_recvs * state->recv_dimension];
    state->self_size = new int[max_sends * state->send_dimension];
    state->send_types = new MPI_Datatype[max_sends];
    state->recv_types = new MPI_Datatype[max_recvs];
    state->send_req = new MPI_Request[10*max_sends];
    state->recv_block_descriptions = new block_description[max_recvs];
    state->recv_from_ranks = new int[max_recvs];
    state->send_block_descriptions = new block_description[max_sends];
    state->send_to_ranks = new int[max_sends];

}

void delete_state_information(redistribution_info* state) {
    MPI_Group_free(&state->send_group);
    MPI_Comm_free(&state->send_comm);

    MPI_Group_free(&state->recv_group);
    MPI_Comm_free(&state->recv_comm);

    delete[] state->send_coords;
    delete[] state->recv_coords;
    delete[] state->send_array_dimension;
    delete[] state->recv_array_dimension;
    delete[] state->self_src;
    delete[] state->self_dst;
    delete[] state->self_size;
    delete[] state->send_types;
    delete[] state->recv_types;
    delete[] state->send_req;

    for (int i = 0; i < state->send_count; i++) {
        delete[] state->send_block_descriptions[i].from;
        delete[] state->send_block_descriptions[i].to;
    }
    delete[] state->send_block_descriptions;
    delete[] state->send_to_ranks;

    for (int i = 0; i < state->recv_count; i++) {
        delete[] state->recv_block_descriptions[i].from;
        delete[] state->recv_block_descriptions[i].to;
    }
    delete[] state->recv_block_descriptions;
    delete[] state->recv_from_ranks;

    delete state;
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

void aggregate_CIs(int num_recorded_values, double* recorded_values, int size, bool* finished)
{
    double current_time;
    LSB_Fold(num_recorded_values, LSB_MAX, &current_time);
    double times_across_ranks[size];
    MPI_Allgather(&current_time, 1, MPI_DOUBLE, times_across_ranks, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    double max_time = 0;
    for (int i = 0; i < size; ++i) {
        max_time = std::max(max_time, times_across_ranks[i]);
    }

    int i = num_recorded_values - 2;
    while (i >= 0 && recorded_values[i] > max_time)
    {
        recorded_values[i + 1] = recorded_values[i];
        --i;
    }
    recorded_values[i + 1] = max_time;

    if (num_recorded_values >= 10) {
        int lower_ci_index = (int) floor(((double)num_recorded_values - 1.96*sqrt((double)num_recorded_values))/2.0);
        int upper_ci_index = (int) ceil(1 + ((double)num_recorded_values + 1.96*sqrt((double)num_recorded_values))/2.0);
        double lower_ci = recorded_values[lower_ci_index];
        double upper_ci = recorded_values[upper_ci_index];
        double median;
        if (num_recorded_values % 2 != 0)
        {
            median = recorded_values[num_recorded_values/2];
        }
        else
        {
            median = (recorded_values[(num_recorded_values-1)/2] + recorded_values[num_recorded_values/2])/2.0;
        }
        double diff = ((median - lower_ci) / median) * 100;
        *finished = diff < 5.0;
    }
}

void configure_LSB_and_sync(int run_idx, int num_warmup, int num_sync, double* win)
{
    if (run_idx == num_warmup)
    {
        LSB_Rec_enable();
    }
    else if (run_idx == num_warmup + num_sync)
    {
        LSB_Fold(0, LSB_MAX, win);
        *win *= 4;
        LSB_Sync_init(MPI_COMM_WORLD, *win);
        LSB_Set_Rparam_string("type", "running");
    }

    if (run_idx < num_warmup + num_sync)
    {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        double err = LSB_Sync();
        LSB_Set_Rparam_double("err", err);
        if (err > 0.0)
        {
            *win += err;
            LSB_Sync_reset(*win);
        }
    }
}