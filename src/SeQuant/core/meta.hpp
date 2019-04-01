//
// Created by Eduard Valeyev on 3/29/18.
//

#ifndef SEQUANT_META_HPP
#define SEQUANT_META_HPP

#include <complex>
#include <memory>
#include <type_traits>

namespace sequant {
namespace meta {

template<typename T>
struct type_printer;

///////// is_shared_ptr /////////

template<class T>
struct is_shared_ptr
    : std::false_type {
};
template<class T>
struct is_shared_ptr<std::shared_ptr<T> >
    : std::true_type {
};
template<class T> static constexpr bool is_shared_ptr_v = is_shared_ptr<T>::value;

///////// is_complex /////////

template<class T>
struct is_complex
    : std::false_type {
};
template<class T>
struct is_complex<std::complex<T> >
    : std::true_type {
};
template<class T> static constexpr bool is_complex_v = is_complex<T>::value;

///////// is_less_than_comparable /////////

template<typename T, typename = std::void_t<>>
struct is_less_than_comparable : public std::false_type {};

template<typename T>
struct is_less_than_comparable<T, std::void_t<decltype(std::declval<const T &>() < std::declval<const T &>())>>
    : public std::true_type {
};

template<typename T> static constexpr bool is_less_than_comparable_v = is_less_than_comparable<T>::value;

///////// is_initializer_list /////////

template<typename T>
struct is_initializer_list : public std::false_type {};

template<typename T>
struct is_initializer_list<std::initializer_list<T>>
    : public std::true_type {
};

template<typename T> static constexpr bool is_initializer_list_v = is_initializer_list<T>::value;

///////// tuple_has_type /////////

template <typename T, typename Tuple>
struct tuple_has_type;

template <typename T, typename... Us>
struct tuple_has_type<T, std::tuple<Us...>> : std::disjunction<std::is_same<T, Us>...> {};

template<typename T, typename Tuple> static constexpr bool tuple_has_type_v = tuple_has_type<T, Tuple>::value;

///////// tuple_index_of /////////

/// tuple_index_of<T, Tuple>::value evaluates to the position of the first appearance of T in Tuple, or, if not found, to a negative number
template <typename T, typename Tuple>
struct tuple_index_of;

template <typename T, typename U, typename... Us>
struct tuple_index_of<T, std::tuple<U, Us...>> {
  static const long value = std::is_same_v<T,U> ? 0 : 1 + tuple_index_of<T, std::tuple<Us...>>::value;
};

template <typename T>
struct tuple_index_of<T, std::tuple<>> {
  static const long value = std::numeric_limits<long>::lowest();
};

template<typename T, typename Tuple> static constexpr bool tuple_index_of_v = tuple_index_of<T, Tuple>::value;

///////// has_memfn_to_latex /////////

template<typename T, typename = std::void_t<>>
struct has_memfn_to_latex : public std::false_type {};

template<typename T>
struct has_memfn_to_latex<T, std::void_t<decltype(static_cast<std::wstring>(std::declval<const T &>().to_latex()))>>
    : public std::true_type {
};

template<typename T> static constexpr bool has_memfn_to_latex_v = has_memfn_to_latex<T>::value;

///////// has_memfn_to_wolfram /////////

template<typename T, typename = std::void_t<>>
struct has_memfn_to_wolfram : public std::false_type {};

template<typename T>
struct has_memfn_to_wolfram<T, std::void_t<decltype(static_cast<std::wstring>(std::declval<const T &>().to_wolfram()))>>
    : public std::true_type {
};

template<typename T> static constexpr bool has_memfn_to_wolfram_v = has_memfn_to_wolfram<T>::value;


}  // namespace meta
}  // namespace sequant

#endif //SEQUANT_META_HPP
