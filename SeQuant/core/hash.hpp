//
// Created by Eduard Valeyev on 2019-01-29.
//

#ifndef SEQUANT_HASH_HPP
#define SEQUANT_HASH_HPP

#include <boost/container_hash/hash.hpp>

#include <type_traits>

#include <boost/container_hash/hash.hpp>
#include <boost/functional/hash.hpp>

#include "meta.hpp"

namespace sequant {

/// @return the version of hashing used by SeQuant, depends on the version of
/// Boost
constexpr int hash_version() {
#if BOOST_VERSION < 108100
  return 1;
#else
  return 2;
#endif
}

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

#if BOOST_VERSION < 108100
template <typename T,
          typename = std::enable_if_t<has_hash_value_member_fn_v<T>>>
auto hash_value(const T& obj) {
  return obj.hash_value();
}
#endif

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
// This skips the first version. Since boost::hash_range is implemented in terms of boost::hash_combine must reimplement
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

#if BOOST_VERSION >= 108100
  boost::hash_combine(seed, hasher(v));
#else  // older boost workarounds
  // in boost 1.78 hash_combine_impl implementation changed
  // https://github.com/boostorg/container_hash/commit/21f2b5e1db1a118c83a3690055c110d0f5637da3
  // probably no longer need these acrobatics
  if constexpr (sizeof(std::size_t) == sizeof(boost::uint32_t) &&
                sizeof(decltype(hasher(v))) == sizeof(boost::uint32_t)) {
    const boost::uint32_t value = hasher(v);
#if BOOST_VERSION >= 107800
    seed = boost::hash_detail::hash_combine_impl<32>::fn(
        static_cast<boost::uint32_t>(seed), value);
#else
    // N.B. seed passed by reference
    boost::hash_detail::hash_combine_impl(
        reinterpret_cast<boost::uint32_t&>(seed), value);
#endif
    return;
  } else if constexpr (sizeof(std::size_t) == sizeof(boost::uint64_t) &&
                       sizeof(decltype(hasher(v))) == sizeof(boost::uint64_t)) {
    const boost::uint64_t value = hasher(v);

#if BOOST_VERSION >= 107800
    seed = boost::hash_detail::hash_combine_impl<64>::fn(
        static_cast<boost::uint64_t>(seed), value);
#else
    // N.B. seed passed by reference
    boost::hash_detail::hash_combine_impl(
        reinterpret_cast<boost::uint64_t&>(seed), value);
#endif
    return;
  } else {
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
#endif  // older boost workarounds

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

/// specialization of boost::hash_range(begin,end) that guarantees to hash a
/// _range_ consisting of a single object to the same value as the hash of that
/// object itself
template <typename It>
std::size_t hash_range(It begin, It end) {
  if (begin != end) {
    std::size_t seed;
    if constexpr (has_hash_value_member_fn_v<std::decay_t<decltype(*begin)>>)
      seed = begin->hash_value();
    else {
      using boost::hash_value;
      std::size_t seed = hash_value(*begin);
    }
    boost::hash_range(seed, begin + 1, end);
    return seed;
  } else
    return boost::hash_range(begin, end);
}

/// redirect to boost::hash_range(seed,begin,end)
template <typename It>
void hash_range(size_t& seed, It begin, It end) {
  boost::hash_range(seed, begin, end);
}

template <typename T>
struct _<
    T, std::enable_if_t<!(detail::has_hash_value_v<T>)&&meta::is_range_v<T>>> {
  std::size_t operator()(T const& v) const { return range(begin(v), end(v)); }
};

template <typename T>
struct _<T, std::enable_if_t<!(
                !(detail::has_hash_value_v<T>)&&meta::is_range_v<T>)>> {
  std::size_t operator()(T const& v) const {
    if constexpr (has_hash_value_member_fn_v<T>)
      return v.hash_value();
    else {
      using boost::hash_value;
      return hash_value(v);
    }
  }
};

template <typename T>
auto value(const T& obj) {
  return _<T>{}(obj);
}

}  // namespace hash

}  // namespace sequant

#endif  // SEQUANT_HASH_HPP
