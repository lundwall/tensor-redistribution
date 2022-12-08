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

    int* total_array_size;

    int* self_src;
    int* self_dst;
    int* self_size;
    int copy_count;
};
