#include <iostream>
#include <send_recv.hpp>
#include <validation.hpp>
#include <liblsb.h>
#include "utils.hpp"

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

void redistribute_by_dimension(redistribution_info* state, int* A, int* A_shape_in, int* B, int* B_shape_in)
{
	switch(state->send_dimension)
	{
		case 2:
		{
			int* _inp_buffer = &A[0];
		    int* _out_buffer = &B[0];
		    NdIndices<2> chunk_num = {1, 1};
		    int A_shape_copy[2]; int B_shape_copy[2];
		    std::memcpy(A_shape_copy, A_shape_in, 2*sizeof(int));
		    std::memcpy(B_shape_copy, B_shape_in, 2*sizeof(int));
		    auto A_shape_tup = as_tuple(A_shape_copy);
		    auto B_shape_tup = as_tuple(B_shape_copy);

		    LSB_Init("2d_redistribute", 0);
    		set_lsb_chunk_size<2>(chunk_num);
		    for (auto idx = 0; idx < state->send_count; ++idx)
		    {
		    	NdIndices<2> from; NdIndices<2> to;
		    	int from_int[2]; int to_int[2];
		    	std::memcpy(from_int, state->send_block_descriptions[idx].from, 2*sizeof(int));
		    	std::memcpy(to_int, state->send_block_descriptions[idx].to, 2*sizeof(int));
		    	auto from_tup = as_tuple(from_int);
		    	auto to_tup = as_tuple(to_int);
		    	LSB_Res();
		    	send<int, 2>(_inp_buffer, state->send_to_ranks[idx], A_shape_tup, from_tup, to_tup, chunk_num);
		    	LSB_Rec(0);
		    }

		    for (auto idx = 0; idx < state->recv_count; ++idx) 
		    {
		    	int from_int[2]; int to_int[2];
		    	std::memcpy(from_int, state->recv_block_descriptions[idx].from, 2*sizeof(int));
		    	std::memcpy(to_int, state->recv_block_descriptions[idx].to, 2*sizeof(int));
		    	auto from_tup = as_tuple(from_int);
		    	auto to_tup = as_tuple(to_int);
		    	LSB_Res();
		    	recv<int, 2>(_out_buffer, state->recv_from_ranks[idx], B_shape_tup, from_tup, to_tup, chunk_num);
		    	LSB_Rec(0);
		    }

		    int* copy_source = _inp_buffer;
		    int* copy_dest = _out_buffer;
		    int src_strides[2];
		    int rcv_strides[2];
			for (auto idx = 0; idx < state->copy_count; ++idx) 
			{
		        src_strides[1] = 1;
		        rcv_strides[1] = 1;
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

		        CopyNDDynamic<int, 1, false, 2>::Dynamic::Copy(
		            copy_source, copy_dest, state->self_size[idx * 2 + 0], src_strides[0], rcv_strides[0], state->self_size[idx * 2 + 1], src_strides[1], rcv_strides[1]
		        );
		    }
		    LSB_Finalize();
    		LSB_chunk_dim_cstr_free_all<2>();
			break;
		}
		default:
			break;
	}
}