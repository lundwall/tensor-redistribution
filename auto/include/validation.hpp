#include "send_recv.hpp"

#ifndef VALIDATION_HPP
#define VALIDATION_HPP

template<typename T, size_t N>
void validate_init(T* array, NdIndices<N> size){
    size_t total = get_product<N>(size);

    for (int i = 0; i < total; i++){
        array[i] = i + 1;
    }
}

template<typename T, size_t N>
bool validate_check(T* orig_array, T* new_array, NdIndices<N> orig_size, NdIndices<N> new_size, NdIndices<N> from, NdIndices<N> to){
    size_t new_total = get_product<N>(new_size);
    bool* checked = new bool[new_total];
    memset(checked, 0, sizeof(bool) * new_total);
    for (int i = 0; i < new_total; i++) {
        checked[i] = new_array[i] == i + 1;
    }

    NdIndices<N> orig_size_postfix = get_tuple_range<1, N + 1>(get_postfix_product<N>(orig_size));
    NdIndices<N> new_size_postfix = get_tuple_range<1, N + 1>(get_postfix_product<N>(new_size));
    
    auto checker = [&](auto... indices) {
        NdIndices<N> ind_tuple = std::make_tuple(indices...);
        size_t orig_ind = sum_tuple<N>(ind_tuple * orig_size_postfix);
        size_t new_ind = sum_tuple<N>(ind_tuple * new_size_postfix);
        checked[new_ind] |= orig_array[orig_ind] == new_array[new_ind];
    };

    n_for<N>(from, to, checker);

    size_t correct_num = 0;
    for (int i = 0; i < new_total; i++) {
        correct_num += checked[i];
    }
    return correct_num == new_total;
}

#endif