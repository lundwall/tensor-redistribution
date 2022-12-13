#include <mpi.h>

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

void read_cfg_file(std::vector<int>& pgrid_send, std::vector<int>& pgrid_recv, std::vector<int>& total_array_size)
{
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
    state->send_req = new MPI_Request[max_sends];
    state->recv_block_descriptions = new block_description[max_recvs];
    state->recv_from_ranks = new int[max_recvs];
    state->send_block_descriptions = new block_description[max_sends];
    state->send_to_ranks = new int[max_sends];

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