//
// Created by Eduard Valeyev on 2019-01-29.
//

#ifndef SEQUANT_HASH_HPP
#define SEQUANT_HASH_HPP

#include <boost/container_hash/hash.hpp>

#include <type_traits>

namespace sequant {

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

/// specialization of boost::hash_range(begin,end) that guarantees to hash a _range_ consisting of a single
/// object to the same value as the hash of that object itself
template <typename It>
std::size_t hash_range(It begin, It end) {
  if (begin != end) {
    using boost::hash_value;
    std::size_t seed = hash_value(*begin);
    boost::hash_range(seed, begin+1, end);
    return seed;
  }
  else
    return boost::hash_range(begin, end);
}

/// redirect to boost::hash_range(seed,begin,end)
template <typename It>
void hash_range(size_t& seed, It begin, It end) {
  boost::hash_range(seed, begin, end);
}


}

#endif //SEQUANT_HASH_HPP
