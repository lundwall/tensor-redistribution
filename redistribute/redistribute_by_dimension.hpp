#include <iostream>
#include <algorithm>
#include <send_recv.hpp>
#include <validation.hpp>
#include <liblsb.h>
#include <omp.h>
#include <math.h>
#include "utils.hpp"
#include "send_recv_5d.hpp"

template<std::size_t... I, typename U>
auto as_tuple(const U &arr, std::index_sequence<I...>) {
    return std::make_tuple(arr[I]...);
}

template<typename T, std::size_t N>
auto as_tuple(const T (&arr)[N]) {
    return as_tuple(arr, std::make_index_sequence<N>{});
}

template<typename T, std::size_t N>
auto as_tuple(const std::array<T, N> &arr) {
    return as_tuple(arr, std::make_index_sequence<N>{});
}

template <size_t N>
auto array_to_tuple(int* input_array)
{
	int array_copy[N];
	std::memcpy(array_copy, input_array, N*sizeof(int));
	return as_tuple(array_copy);
}

template<size_t N, size_t dim, typename... Ts>
void CopyNDDynamicHelper(int copy_idx, int* src_strides, int* rcv_strides, int* self_size, Ts... args)
{
	if constexpr (dim == N)
	{
		CopyNDDynamic<int, 1, false, N>::Dynamic::Copy(args...);
	}
	else
	{
		CopyNDDynamicHelper<N, dim+1>(copy_idx, src_strides, rcv_strides, self_size, args..., self_size[copy_idx*N+dim], src_strides[dim], rcv_strides[dim]);
	}
}

template <size_t N>
void redistribute_by_dimension_template(redistribution_info* state, int* A, int* A_shape_in, int* B, int* B_shape_in, std::string MODE, int num_chunks)
{
	int* _inp_buffer = &A[0];
    int* _out_buffer = &B[0];
    int* chunk = new int[N];

    for (int i = 0; i < N; i++)
    {
    	chunk[i] = 1;
    }

    auto A_shape_tup = array_to_tuple<N>(A_shape_in);
    auto B_shape_tup = array_to_tuple<N>(B_shape_in);
    auto chunk_num_tup = array_to_tuple<N>(chunk);

	size_t allocate_number = std::max(state->send_count, 1);
	int** send_buffers = new int*[allocate_number];

	NdIndices<N>* from_tup_send = new NdIndices<N>[state->send_count];
	NdIndices<N>* to_tup_send = new NdIndices<N>[state->send_count];
	NdIndices<N>* from_tup_recv = new NdIndices<N>[state->recv_count]; 
	NdIndices<N>* to_tup_recv = new NdIndices<N>[state->recv_count];

	for (auto idx = 0; idx < state->send_count; ++idx)
	{
		from_tup_send[idx] = array_to_tuple<N>(state->send_block_descriptions[idx].from);
		to_tup_send[idx] = array_to_tuple<N>(state->send_block_descriptions[idx].to);
		NdIndices<N> range = to_tup_send[idx] - from_tup_send[idx];
		int sending_total = get_product<N>(range);
		send_buffers[idx] = new int[sending_total];
	}
	for (auto idx = 0; idx < state->send_count; ++idx)
	{
		if (MODE == "manual")
		{
			send_5d(_inp_buffer, state->send_to_ranks[idx], A_shape_in, state->send_block_descriptions[idx].from, state->send_block_descriptions[idx].to, num_chunks, &state->send_req[num_chunks*idx], send_buffers[idx]);
		}
		else if (MODE == "datatype")
		{
			MPI_Isend(_inp_buffer, 1, state->send_types[idx], state->send_to_ranks[idx], 0, MPI_COMM_WORLD, &state->send_req[idx]);
		}
	}

	for (auto idx = 0; idx < state->recv_count; ++idx) 
	{
		if (MODE == "manual")
		{
			recv_5d(_out_buffer, state->recv_from_ranks[idx], B_shape_in, state->recv_block_descriptions[idx].from, state->recv_block_descriptions[idx].to, num_chunks);
		}
		else if (MODE == "datatype")
		{
			MPI_Recv(_out_buffer, 1, state->recv_types[idx], state->recv_from_ranks[idx], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	MPI_Waitall(state->send_count * num_chunks, state->send_req, MPI_STATUSES_IGNORE);

	for (auto idx = 0; idx < state->send_count; ++idx)
	{
		delete[] send_buffers[idx];
	}

	delete[] send_buffers;
	send_buffers = nullptr;

	delete[] chunk;

	delete[] from_tup_send;
	delete[] to_tup_send;
	delete[] from_tup_recv;
	delete[] to_tup_recv;

	int* copy_source = _inp_buffer;
	int* copy_dest = _out_buffer;
	int src_strides[N];
	int rcv_strides[N];
	for (auto idx = 0; idx < state->copy_count; ++idx) 
	{
		src_strides[N-1] = 1;
		rcv_strides[N-1] = 1;
		int src_stride = 1;
		int rcv_stride = 1;
		for (auto dim_idx = state->send_dimension-1; dim_idx >= 0; dim_idx--)
		{
			copy_source += state->self_src[idx * state->send_dimension + dim_idx] * src_stride;
			src_strides[dim_idx] = src_stride;
			src_stride *= A_shape_in[dim_idx];
		}

		for (auto dim_idx = state->recv_dimension=1; dim_idx >= 0; dim_idx--)
		{
			copy_dest += state->self_dst[idx * state->recv_dimension + dim_idx] * rcv_stride;
			rcv_strides[dim_idx] = rcv_stride;
			rcv_stride *= B_shape_in[dim_idx];
		}

		CopyNDDynamicHelper<N, 1>(idx, src_strides, rcv_strides, state->self_size, copy_source, copy_dest, state->self_size[0], src_strides[0], rcv_strides[0]);
	}
}

void redistribute_by_dimension(redistribution_info* state, int* A, int* A_shape_in, int* B, int* B_shape_in, std::string MODE, int num_chunks)
{
	switch(state->send_dimension)
	{
		case 2:
		{
			redistribute_by_dimension_template<2>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		case 3:
		{
			redistribute_by_dimension_template<3>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		case 4:
		{
			redistribute_by_dimension_template<4>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		case 5:
		{
			redistribute_by_dimension_template<5>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		case 6:
		{
			redistribute_by_dimension_template<6>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		case 7:
		{
			redistribute_by_dimension_template<7>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		case 8:
		{
			redistribute_by_dimension_template<8>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		case 9:
		{
			redistribute_by_dimension_template<9>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		case 10:
		{
			redistribute_by_dimension_template<10>(state, A, A_shape_in, B, B_shape_in, MODE, num_chunks);
			break;
		}
		default:
			break;
	}
}
