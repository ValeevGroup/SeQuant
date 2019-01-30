//
// Created by Eduard Valeyev on 2019-01-29.
//

#ifndef SEQUANT2_HASH_HPP
#define SEQUANT2_HASH_HPP

#include <type_traits>

namespace sequant2 {

namespace detail {
template<typename T, typename Enabler = void>
struct has_hash_value_member_fn_helper : public std::false_type {};
template<typename T>
struct has_hash_value_member_fn_helper<T, std::void_t<decltype(std::declval<const T &>().hash_value())>>
    : public std::true_type {
};
}

template<typename T>
static constexpr bool has_hash_value_member_fn_v = detail::has_hash_value_member_fn_helper<T>::value;

template<typename T>
auto hash_value(const T &obj, std::enable_if_t<has_hash_value_member_fn_v<T>> * = nullptr) {
  return obj.hash_value();
}

}

#endif //SEQUANT2_HASH_HPP
