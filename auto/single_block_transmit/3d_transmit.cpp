#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <omp.h>
#include <string>
#include <tuple>
#include <cassert>
#include <send_recv.hpp>

int main(int argc, char** argv){
    constexpr size_t RUNS = 1;
    constexpr size_t N = 3;
    using T = int;
    
    std::string file_name = std::to_string(N) + std::string("d_transmit_auto_") + std::string(std::getenv("OMP_NUM_THREADS"));

    NdIndices<N> chunk_num = {2, 2, 3};

    NdIndices<N> current_size = {5, 6, 9};
    NdIndices<N> new_size = {9, 9, 9};

    NdIndices<N> from = {1, 2, 3};
    NdIndices<N> to = {3, 6, 9};

    T* current_array;
    T* new_array;
    int rank;
    rank = init<T, N>(argc, argv, file_name, chunk_num, current_array, new_array, current_size, new_size);

    for (int i = 0; i < get_product<N>(current_size); i++){
        current_array[i] = i;
    }
    for (int i = 0; i < get_product<N>(new_size); i++){
        new_array[i] = 0;
    }

    if (rank == 0) {
        for (int k = 0; k < RUNS; ++k) {
            LSB_Res();
            send<T, N>(current_array, 1, current_size, from, to, chunk_num);
            LSB_Rec(k);
        }
    }

    if(rank == 1){
        for (int k = 0; k < RUNS; ++k){
            LSB_Res();
            recv<T, N>(new_array, 0, new_size, from, to, chunk_num);
            LSB_Rec(k);
        }
        for (int i = 0; i < std::get<0>(new_size); i++){
            for (int j = 0; j < std::get<1>(new_size); j++){
                for (int k = 0; k < std::get<2>(new_size); k++){
                    std::cout << new_array[i * std::get<1>(new_size) * std::get<2>(new_size) + j * std::get<2>(new_size) + k] << "\t";
                }
                std::cout << "\n";
            }
            std::cout << "\n\n";
        }
    }

    term(current_array, new_array);
}
