#include <mpi.h>
#include <stdlib.h>
#include <cstring>
#include <liblsb.h>
#include <omp.h>
#include <string>
#include <tuple>
#include <cassert>

#ifndef SEND_RECV_HPP
#define SEND_RECV_HPP

using std::size_t;

template <typename T, size_t N>
struct tuple_of {
    template <size_t>
    using Capture = T;

    template <size_t... Is>
    static auto helper(std::index_sequence<Is...>) {
        return std::tuple<Capture<Is>...>{};
    }

    using type = decltype(helper(std::make_index_sequence<N>{}));
};

template <typename T, size_t N>
using tuple_of_t = typename tuple_of<T, N>::type;

namespace tuple_arith_impl{
template <typename Op, typename T, std::size_t... Is>
inline auto helper2(const Op& op, const T& l, const T& r, const std::index_sequence<Is...>&) {
    return std::make_tuple(op(std::get<Is>(l), std::get<Is>(r))...);
}

template <typename Op, typename T>
inline auto helper1(const Op& op, const T& l, const T& r) {
    return helper2(op, l, r, std::make_index_sequence<std::tuple_size<T>{}>{});
}
}  // namespace NdIndices_impl

template <typename ... Ts>
inline auto operator+(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
    return tuple_arith_impl::helper1(std::plus<>{}, l, r);
}

template <typename ... Ts>
inline auto operator-(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
    return tuple_arith_impl::helper1(std::minus<>{}, l, r);
}

template <typename ... Ts>
inline auto operator*(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
    return tuple_arith_impl::helper1(std::multiplies<>{}, l, r);
}

template <typename ... Ts>
inline auto operator/(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
    return tuple_arith_impl::helper1(std::divides<>{}, l, r);
}

template <typename ... Ts>
inline auto operator%(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
    return tuple_arith_impl::helper1(std::modulus<>{}, l, r);
}

template <typename ... Ts>
inline auto operator<(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
    return tuple_arith_impl::helper1(std::less<>{}, l, r);
}

template <size_t N, size_t... Is>
constexpr auto make_index_range_helper(std::index_sequence<Is...>) { 
    return std::index_sequence<N + Is...>{}; 
}

template <size_t S, size_t E>
using make_index_range = decltype(make_index_range_helper<S>(std::make_index_sequence<E - S>{}));

template <typename... Ts, size_t... Is>
inline auto get_tuple_helper(std::tuple<Ts...> inp, std::index_sequence<Is...> const&){
    return std::make_tuple(std::get<Is>(inp)...);
}

template <size_t S, size_t E, typename... Ts>
inline auto get_tuple_range(std::tuple<Ts...> inp){
    make_index_range<S, E> range;
    return get_tuple_helper(inp, range);
}

template <typename T>
struct to_MPI_type_helper {
    template <typename Tl>
    static constexpr MPI_Datatype helper(typename std::enable_if<std::is_same_v<Tl, int>, int>::type = 0){
        return MPI_INT;
    };

    static constexpr MPI_Datatype value = helper<T>();
};

template <typename T>
constexpr MPI_Datatype to_MPI_type = to_MPI_type_helper<T>::value;

template <size_t N>
using NdIndices = tuple_of_t<size_t, N>;

template <size_t N>
struct create_zero_NdIndices {
    template <size_t>
    static constexpr int capture = 0;

    template <size_t... Is>
    static constexpr auto helper(std::index_sequence<Is...>) {
        return NdIndices<N>{capture<Is>...};
    }

    static constexpr NdIndices<N> value = helper(std::make_index_sequence<N>{});
};

template <size_t N>
constexpr auto zero_NdIndices = create_zero_NdIndices<N>::value;

template <size_t N, size_t... Is>
inline size_t get_product_by_index(NdIndices<N> inp, const std::index_sequence<Is...>&) { 
    return (std::get<Is>(inp) * ...); 
}

template <size_t N, size_t S = 0, size_t E = N>
inline size_t get_product(NdIndices<N> inp){
    make_index_range<S, E> range;
    return get_product_by_index<N>(inp, range);
}

template <size_t N,  size_t... Is>
inline NdIndices<N + 1> get_postfix_product_helper(NdIndices<N> inp, const std::index_sequence<Is...>&){
    return std::make_tuple(get_product<N, Is, N>(inp)..., 1);
}

template <size_t N>
inline NdIndices<N + 1> get_postfix_product(NdIndices<N> inp){
    auto n_ind = std::make_index_sequence<N>();
    return get_postfix_product_helper<N>(inp, n_ind);
}

template <size_t N,  size_t... Is>
inline size_t sum_tuple_helper(NdIndices<N> inp, const std::index_sequence<Is...>&){
    return (std::get<Is>(inp) + ...);
}

template <size_t N>
inline size_t sum_tuple(NdIndices<N> inp){
    auto n_ind = std::make_index_sequence<N>();
    return sum_tuple_helper<N>(inp, n_ind);
}

template <size_t N,  size_t... Is>
inline void check_divisible_helper(NdIndices<N> top, NdIndices<N> bot, const std::index_sequence<Is...>&){
    if((std::get<Is>(top % bot) + ...)){
        throw std::runtime_error("Chunk sizes cannot evenly divide");
    }
}

template <size_t N>
inline void check_divisible(NdIndices<N> top, NdIndices<N> bot){
    auto n_ind = std::make_index_sequence<N>();
    return check_divisible_helper<N>(top, bot, n_ind);
}

template <size_t I>
struct LSB_chunk_dim_cstr{
    inline static char* static_str;

    static char* set(const std::string& str){
        static_str = strdup(str.c_str());
        return static_str;
    }

    static char* get(){
        return static_str;
    }
};

template <size_t... Is>
inline void LSB_chunk_dim_cstr_free_all_helper(const std::index_sequence<Is...>&){
    (free(LSB_chunk_dim_cstr<Is>::get()), ...);
}

template <ssize_t N>
inline void LSB_chunk_dim_cstr_free_all(){
    auto n_ind = std::make_index_sequence<N>();
    return LSB_chunk_dim_cstr_free_all_helper(n_ind);
}

template <size_t N,  size_t... Is>
inline void set_lsb_chunk_size_helper(NdIndices<N> inp, const std::index_sequence<Is...>&){
    (LSB_Set_Rparam_long(LSB_chunk_dim_cstr<Is>::set(std::string("chunk_dim_") + std::to_string(Is)), std::get<Is>(inp)), ...);
}

template <size_t N>
inline void set_lsb_chunk_size(NdIndices<N> inp){
    auto n_ind = std::make_index_sequence<N>();
    return set_lsb_chunk_size_helper<N>(inp, n_ind);
}

template<size_t N, size_t I, typename Callable>
inline constexpr void n_for_helper(NdIndices<N> begin, NdIndices<N> end, Callable&& c) {
    for(size_t i = std::get<I>(begin); i != std::get<I>(end); ++i){
        if constexpr(I == N - 1){
            c(i);
        } else {
            auto bind_an_argument = [i, &c](auto... args) {
                c(i, args...);
            };
            n_for_helper<N, I + 1>(begin, end, bind_an_argument);
        }
    }
}

template<size_t N, typename Callable>
inline constexpr void n_for(NdIndices<N> begin, NdIndices<N> end, Callable&& c) {
    n_for_helper<N, 0>(begin, end, c);
}

template<size_t N, size_t I, typename Callable>
inline constexpr void n_for_task_helper(NdIndices<N> begin, NdIndices<N> end, Callable&& c) {
    for(size_t i = std::get<I>(begin); i != std::get<I>(end); ++i){
        if constexpr(I == N - 1){
#pragma omp task
            c(i);
        } else {
            auto bind_an_argument = [i, &c](auto... args) {
                c(i, args...);
            };
            n_for_task_helper<N, I + 1>(begin, end, bind_an_argument);
        }
    }
}

template<size_t N, typename Callable>
inline constexpr void n_task_for(NdIndices<N> begin, NdIndices<N> end, Callable&& c) {
#pragma omp parallel
#pragma omp single
    n_for_task_helper<N, 0>(begin, end, c);
}

template <typename T, size_t N>
void send(T* source, int other_rank, NdIndices<N> current_size, NdIndices<N> from, NdIndices<N> to, NdIndices<N> chunk_num){
    NdIndices<N> range = to - from;
    check_divisible<N>(range, chunk_num);

    NdIndices<N> chunk_size = range / chunk_num;
    NdIndices<N> size_postfix = get_tuple_range<1, N + 1>(get_postfix_product<N>(current_size));
    NdIndices<N> chunk_num_postfix = get_tuple_range<1, N + 1>(get_postfix_product<N>(chunk_num));
    NdIndices<N> chunk_size_postfix = get_tuple_range<1, N + 1>(get_postfix_product<N>(chunk_size));
    NdIndices<N - 1> chunk_size_1 = get_tuple_range<0, N - 1>(chunk_size);

    size_t chunk_total = get_product<N>(chunk_size);
    size_t sending_total = get_product<N>(range);
    size_t row_length = std::get<N - 1>(chunk_size);
    
    MPI_Datatype datatype = to_MPI_type<T>;
    T* buffer = new T[sending_total];
    MPI_Request send_req;
    bool is_first = true;

    // indices in chunk_num
    auto chunk_iter = [&](auto... indices) {
        NdIndices<N> chunk_ind_tuple = std::make_tuple(indices...);
        NdIndices<N> chunk_start = chunk_ind_tuple * chunk_size;
        size_t chunk_id = sum_tuple<N>(chunk_ind_tuple * chunk_num_postfix);
        // indices in chunk_size[:-1]
        auto ele_iter = [&](auto... indices) {
            NdIndices<N> ele_ind_tuple = std::make_tuple(indices..., 0);
            size_t source_start = sum_tuple<N>((from + chunk_start + ele_ind_tuple) * size_postfix);
            size_t buffer_start = chunk_total * chunk_id + sum_tuple<N>(ele_ind_tuple * chunk_size_postfix);
            memcpy(buffer + buffer_start, source + source_start, sizeof(T) * row_length);
        };

        n_for<N - 1>(zero_NdIndices<N - 1>, chunk_size_1, ele_iter);

        if (!is_first) {
            MPI_Waitall(1, &send_req, MPI_STATUSES_IGNORE);
        }else{
            is_first = false;
        }
        MPI_Isend(&(buffer[chunk_total * chunk_id]), chunk_total, datatype, other_rank, chunk_id, MPI_COMM_WORLD, &send_req);
    };

    n_for<N>(zero_NdIndices<N>, chunk_num, chunk_iter);
    MPI_Waitall(1, &send_req, MPI_STATUSES_IGNORE);

    delete[] buffer;
}

template <typename T, size_t N>
void recv(T* target, int other_rank, NdIndices<N> new_size, NdIndices<N> from, NdIndices<N> to, NdIndices<N> chunk_num){
    NdIndices<N> range = to - from;
    check_divisible<N>(range, chunk_num);

    NdIndices<N> chunk_size = range / chunk_num;
    NdIndices<N> size_postfix = get_tuple_range<1, N + 1>(get_postfix_product<N>(new_size));
    NdIndices<N> chunk_num_postfix = get_tuple_range<1, N + 1>(get_postfix_product<N>(chunk_num));
    NdIndices<N> chunk_size_postfix = get_tuple_range<1, N + 1>(get_postfix_product<N>(chunk_size));
    NdIndices<N - 1> chunk_size_1 = get_tuple_range<0, N - 1>(chunk_size);

    size_t chunk_total = get_product<N>(chunk_size);
    size_t receiving_total = get_product<N>(range);
    size_t row_length = std::get<N - 1>(chunk_size);
    
    MPI_Datatype datatype = to_MPI_type<T>;
    T* buffer = new T[receiving_total];
    MPI_Request recv_req;
    bool is_first = true;
    NdIndices<N> prev_chunk_ind_tuple;

    auto transmit_from_buffer = [&](NdIndices<N> chunk_ind_tuple) {
        NdIndices<N> chunk_start = chunk_ind_tuple * chunk_size;
        size_t chunk_id = sum_tuple<N>(chunk_ind_tuple * chunk_num_postfix);
        MPI_Waitall(1, &recv_req, MPI_STATUSES_IGNORE);
        auto ele_iter = [&](auto... indices) {
            NdIndices<N> ele_ind_tuple = std::make_tuple(indices..., 0);
            size_t target_start = sum_tuple<N>((from + chunk_start + ele_ind_tuple) * size_postfix);
            size_t buffer_start = chunk_total * chunk_id + sum_tuple<N>(ele_ind_tuple * chunk_size_postfix);
            memcpy(target + target_start, buffer + buffer_start, sizeof(T) * row_length);
        };

        n_for<N - 1>(zero_NdIndices<N - 1>, chunk_size_1, ele_iter);
    };

    // indices in chunk_num
    auto chunk_iter = [&](auto... indices) {
        NdIndices<N> chunk_ind_tuple = std::make_tuple(indices...);
        size_t chunk_id = sum_tuple<N>(chunk_ind_tuple * chunk_num_postfix);

        MPI_Irecv(&(buffer[chunk_total * chunk_id]), chunk_total, datatype, other_rank, chunk_id, MPI_COMM_WORLD, &recv_req);

        if (!is_first) {
            transmit_from_buffer(prev_chunk_ind_tuple);
        }else{
            is_first = false;
        }
        prev_chunk_ind_tuple = chunk_ind_tuple;
    };

    n_for<N>(zero_NdIndices<N>, chunk_num, chunk_iter);
    transmit_from_buffer(prev_chunk_ind_tuple);

    delete[] buffer;
}

template <typename T, size_t N>
int init(int argc, char** argv, const std::string& name, NdIndices<N> chunk_num, T*& orig_arr, T*& new_arr, NdIndices<N> current_size, NdIndices<N> new_size){
    int size, rank, received_threads;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &received_threads);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    LSB_Init(name.c_str(), 0);
    LSB_Set_Rparam_int("rank", rank);
    set_lsb_chunk_size<N>(chunk_num);


    // size should be 2!
    if (size != 2) {
        throw std::runtime_error("When testing only 1 block transmission, only 2 processors are needed");
    }

    char processor_name[256];
    int len_processor_name = 0;
    MPI_Get_processor_name(processor_name, &len_processor_name);
    std::cout << processor_name << std::endl;


    size_t current_total = get_product<N>(current_size);
    size_t new_total = get_product<N>(new_size);

    orig_arr = new T[current_total];
    new_arr = new T[new_total];

    return rank;
}

template <typename T, size_t N>
void term(T*& orig_arr, T*& new_arr){
    LSB_Finalize();
    MPI_Finalize();

    LSB_chunk_dim_cstr_free_all<N>();

    delete[] orig_arr; 
    orig_arr = nullptr;
    delete[] new_arr;
    new_arr = nullptr;
}

#endif