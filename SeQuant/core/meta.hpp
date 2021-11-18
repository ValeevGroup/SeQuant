//
// Created by Eduard Valeyev on 3/29/18.
//

#ifndef SEQUANT_META_HPP
#define SEQUANT_META_HPP

#include <complex>
#include <memory>
#include <type_traits>
#include <range/v3/range/access.hpp>

#include <range/v3/range/access.hpp>

namespace sequant {
namespace meta {

template <typename T>
struct type_printer;

///////// is_detected ///////////

struct nonesuch {
  ~nonesuch() = delete;
  nonesuch(nonesuch const &) = delete;
  void operator=(nonesuch const &) = delete;
};

namespace detail {
template <class Default, class AlwaysVoid, template <class...> class Op,
          class... Args>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};

template <class Default, template <class...> class Op, class... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> {
  using value_t = std::true_type;
  using type = Op<Args...>;
};

}  // namespace detail

template <template <class...> class Op, class... Args>
using is_detected =
    typename detail::detector<nonesuch, void, Op, Args...>::value_t;

template <template <class...> class Op, class... Args>
using detected_t = typename detail::detector<nonesuch, void, Op, Args...>::type;

template <class Default, template <class...> class Op, class... Args>
using detected_or = detail::detector<Default, void, Op, Args...>;

template <template <class...> class Op, class... Args>
constexpr inline bool is_detected_v = is_detected<Op, Args...>::value;

template <class Default, template <class...> class Op, class... Args>
using detected_or_t = typename detected_or<Default, Op, Args...>::type;

template <class Expected, template <class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

template <class Expected, template <class...> class Op, class... Args>
constexpr inline bool is_detected_exact_v =
    is_detected_exact<Expected, Op, Args...>::value;

template <class To, template <class...> class Op, class... Args>
using is_detected_convertible =
    std::is_convertible<detected_t<Op, Args...>, To>;

template <class To, template <class...> class Op, class... Args>
constexpr inline bool is_detected_convertible_v =
    is_detected_convertible<To, Op, Args...>::value;

///////// is_shared_ptr /////////

template <class T>
struct is_shared_ptr : std::false_type {};
template <class T>
struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {};
template <class T>
static constexpr bool is_shared_ptr_v = is_shared_ptr<T>::value;

///////// is_complex /////////

template <class T>
struct is_complex : std::false_type {};
template <class T>
struct is_complex<std::complex<T>> : std::true_type {};
template <class T>
static constexpr bool is_complex_v = is_complex<T>::value;

///////// is_less_than_comparable /////////

template <typename T, typename = std::void_t<>>
struct is_less_than_comparable : public std::false_type {};

template <typename T>
struct is_less_than_comparable<T,
                               std::void_t<decltype(std::declval<const T &>() <
                                                    std::declval<const T &>())>>
    : public std::true_type {};

template <typename T>
static constexpr bool is_less_than_comparable_v =
    is_less_than_comparable<T>::value;

///////// is_initializer_list /////////

template <typename T>
struct is_initializer_list : public std::false_type {};

template <typename T>
struct is_initializer_list<std::initializer_list<T>> : public std::true_type {};

template <typename T>
static constexpr bool is_initializer_list_v = is_initializer_list<T>::value;

///////// tuple_has_type /////////

template <typename T, typename Tuple>
struct tuple_has_type;

template <typename T, typename... Us>
struct tuple_has_type<T, std::tuple<Us...>>
    : std::disjunction<std::is_same<T, Us>...> {};

template <typename T, typename Tuple>
static constexpr bool tuple_has_type_v = tuple_has_type<T, Tuple>::value;

///////// tuple_index_of /////////

/// tuple_index_of<T, Tuple>::value evaluates to the position of the first appearance of T in Tuple, or, if not found, to a negative number
template <typename T, typename Tuple>
struct tuple_index_of;

template <typename T, typename U, typename... Us>
struct tuple_index_of<T, std::tuple<U, Us...>> {
  static const long value =
      std::is_same_v<T, U> ? 0
                           : 1 + tuple_index_of<T, std::tuple<Us...>>::value;
};

template <typename T>
struct tuple_index_of<T, std::tuple<>> {
  static const long value = std::numeric_limits<long>::lowest();
};

template <typename T, typename Tuple>
static constexpr bool tuple_index_of_v = tuple_index_of<T, Tuple>::value;

///////// has_memfn_to_latex /////////

template <typename T, typename = std::void_t<>>
struct has_memfn_to_latex : public std::false_type {};

template <typename T>
struct has_memfn_to_latex<T, std::void_t<decltype(static_cast<std::wstring>(
                                 std::declval<const T &>().to_latex()))>>
    : public std::true_type {};

template <typename T>
static constexpr bool has_memfn_to_latex_v = has_memfn_to_latex<T>::value;

///////// has_memfn_to_wolfram /////////

template <typename T, typename = std::void_t<>>
struct has_memfn_to_wolfram : public std::false_type {};

template <typename T>
struct has_memfn_to_wolfram<T, std::void_t<decltype(static_cast<std::wstring>(
                                   std::declval<const T &>().to_wolfram()))>>
    : public std::true_type {};

template <typename T>
static constexpr bool has_memfn_to_wolfram_v = has_memfn_to_wolfram<T>::value;

/// is_range

namespace is_range_impl {  // detects presence of std::{begin,end}

template <class T>
using std_begin_t = decltype(std::begin(std::declval<T &>()));

template <class T>
using std_end_t = decltype(std::begin(std::declval<T &>()));

template <class T>
using ranges_begin_t = decltype(ranges::begin(std::declval<T &>()));

template <class T>
using ranges_end_t = decltype(ranges::end(std::declval<T &>()));

}  // namespace is_range_impl

template <typename T>
static constexpr bool is_range_v =
    (is_detected_v<is_range_impl::std_begin_t, T> &&
     is_detected_v<is_range_impl::std_end_t, T>) ||
    (is_detected_v<is_range_impl::ranges_begin_t, T> &&
     is_detected_v<is_range_impl::ranges_end_t, T>);

}  // namespace meta
}  // namespace sequant

#endif //SEQUANT_META_HPP
