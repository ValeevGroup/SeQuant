//
// Created by Eduard Valeyev on 2019-01-29.
//

#ifndef SEQUANT_HASH_HPP
#define SEQUANT_HASH_HPP

#include <boost/container_hash/hash.hpp>

#include <type_traits>

#include <boost/functional/hash.hpp>
#include <boost/container_hash/hash.hpp>

#include "meta.hpp"

namespace sequant {

class Expr;

namespace detail {
template <typename T, typename Enabler = void>
struct has_hash_value_member_fn_helper : public std::false_type {};
template <typename T>
struct has_hash_value_member_fn_helper<
    T, std::void_t<decltype(std::declval<const T&>().hash_value())>>
    : public std::true_type {};
}  // namespace detail

template <typename T>
static constexpr bool has_hash_value_member_fn_v =
    detail::has_hash_value_member_fn_helper<T>::value;

template <typename T>
auto hash_value(const T& obj,
                std::enable_if_t<has_hash_value_member_fn_v<T>>* = nullptr) {
  return obj.hash_value();
}

namespace detail {

template <typename T, typename = std::void_t<>>
struct has_boost_hash_value : std::false_type {};

template <typename T>
struct has_boost_hash_value<
    T, std::void_t<decltype(boost::hash_value(std::declval<T>()))>>
    : std::true_type {};

template <typename T>
constexpr bool has_boost_hash_value_v = has_boost_hash_value<const T&>::value;

template <typename T, typename = std::void_t<>>
struct has_hash_value : std::false_type {};

template <typename T>
struct has_hash_value<
    T, std::void_t<decltype(hash_value(std::declval<const T&>()))>>
    : std::true_type {};

template <typename T>
constexpr bool has_hash_value_v = has_hash_value<T>::value;

}  // namespace detail

using boost::hash_value;

// clang-format off
// rationale:
// boost::hash_combine is busted ... it dispatches to one of 3 implementations (all line numbers refer to boost 1.72.0):
// - templated: container_hash/hash.hpp:309
// - 32-bit: container_hash/hash.hpp:315
// - 64-bit: container_hash/hash.hpp:336
// n.b. The latter is disabled if int64_t is not available or 32-bit gcc is used
// Somehow Macos clang dispatches to the first version, but everywhere else (Linux) dispatch ends up to the last one.
// This skips the first version. Since boost::hash_range is implemented in terms of boost::hash_range must reimplement
// that also.
// clang-format on

namespace hash {

template <typename T, typename Enabler = void>
struct _;

template <class T>
inline void combine(std::size_t& seed, T const& v) {
//  std::size_t seed_ref = seed;
//  boost::hash_combine(seed_ref, v);
  _<T> hasher;
  if constexpr (sizeof(std::size_t) == sizeof(boost::uint32_t) &&
                sizeof(decltype(hasher(v))) == sizeof(boost::uint32_t)) {
    const boost::uint32_t value = hasher(v);
    return boost::hash_detail::hash_combine_impl(
        reinterpret_cast<boost::uint32_t&>(seed), value);
  } else if constexpr (sizeof(std::size_t) == sizeof(boost::uint64_t) &&
                       sizeof(decltype(hasher(v))) == sizeof(boost::uint64_t)) {
    const boost::uint64_t value = hasher(v);
    return boost::hash_detail::hash_combine_impl(
        reinterpret_cast<boost::uint64_t&>(seed), value);
  } else {
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
//  assert(seed == seed_ref);
}

template <class It>
inline std::size_t range(It first, It last) {
//  const std::size_t seed_ref = boost::hash_range(first, last);
  std::size_t seed = 0;

  for (; first != last; ++first) {
    hash::combine(seed, *first);
  }

//  assert(seed == seed_ref);
  return seed;
}

template <class It>
inline void range(std::size_t& seed, It first, It last) {
//  std::size_t seed_ref = seed;
//  boost::hash_range(seed_ref, first, last);
  for (; first != last; ++first) {
    hash::combine(seed, *first);
  }
//  assert(seed == seed_ref);
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

template <typename T>
struct _<T, std::enable_if_t<!(detail::has_hash_value_v<T>) && meta::is_range_v<T>>> {
  std::size_t operator()(T const& v) const { return range(ranges::begin(v), ranges::end(v)); }
};

template <typename T>
struct _<T, std::enable_if_t<!(!(detail::has_hash_value_v<T>) && meta::is_range_v<T>)>> {
  std::size_t operator()(T const& v) const {
    using boost::hash_value;
    return hash_value(v);
  }
};

template <typename T>
auto value(const T& obj) {
  return _<T>{}(obj);
}

}  // namespace hash

}  // namespace sequant

#endif //SEQUANT_HASH_HPP
