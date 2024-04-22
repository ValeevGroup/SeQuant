//
// Created by Eduard Valeyev on 2019-01-29.
//

#ifndef SEQUANT_HASH_HPP
#define SEQUANT_HASH_HPP

#ifdef SEQUANT_USE_SYSTEM_BOOST_HASH
#include <boost/version.hpp>
#define SEQUANT_BOOST_VERSION BOOST_VERSION
#include <boost/container_hash/hash.hpp>
namespace sequant_boost = boost;
#else
#include <SeQuant/external/boost/container_hash/hash.hpp>
#endif

#if SEQUANT_BOOST_VERSION < 108100
#error "SeQuant requires Boost 1.81 or later for hashing"
#endif

#include <SeQuant/core/meta.hpp>

namespace sequant {

namespace hash {

/// the hashing versions known to SeQuant (N.B. hashing changed in Boost 1.81)
enum class Impl { BoostPre181 = 1, Boost181OrLater = 2 };

}  // namespace hash

/// @return the version of hashing used by SeQuant, depends on the version of
/// Boost
constexpr hash::Impl hash_version() { return hash::Impl::Boost181OrLater; }

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

namespace detail {

template <typename T, typename = std::void_t<>>
struct has_boost_hash_value : std::false_type {};

template <typename T>
struct has_boost_hash_value<T, std::void_t<decltype(sequant_boost::hash_value(
                                   std::declval<const T&>()))>>
    : std::true_type {};

template <typename T, typename = std::void_t<>>
struct has_hash_value : std::false_type {};

template <typename T>
struct has_hash_value<
    T, std::void_t<decltype(hash_value(std::declval<const T&>()))>>
    : std::true_type {};

}  // namespace detail

template <typename T>
constexpr bool has_boost_hash_value_v =
    detail::has_boost_hash_value<const T&>::value;

template <typename T>
constexpr bool has_hash_value_v = detail::has_hash_value<T>::value;

// hash_value specialization for types that have a hash_value member function
template <typename T,
          std::enable_if_t<has_hash_value_member_fn_v<T>, short> = 0>
auto hash_value(const T& obj) {
  return obj.hash_value();
}

// hash_value specialization that don't have a hash_value member function but
// have an applicable boost::hash_value function
template <typename T, std::enable_if_t<!has_hash_value_member_fn_v<T> &&
                                           has_boost_hash_value_v<T>,
                                       int> = 0>
auto hash_value(const T& obj) {
  return sequant_boost::hash_value(obj);
}

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
  _<T> hasher;

  sequant_boost::hash_combine(seed, hasher(v));

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
    std::size_t seed = hash_value(*begin);
    sequant_boost::hash_range(seed, begin + 1, end);
    return seed;
  } else
    return sequant_boost::hash_range(begin, end);
}

/// redirect to boost::hash_range(seed,begin,end)
template <typename It>
void hash_range(size_t& seed, It begin, It end) {
  sequant_boost::hash_range(seed, begin, end);
}

template <typename T>
struct _<T, std::enable_if_t<!(has_hash_value_v<T>) && meta::is_range_v<T>>> {
  std::size_t operator()(T const& v) const { return range(begin(v), end(v)); }
};

template <typename T>
struct _<T,
         std::enable_if_t<!(!(has_hash_value_v<T>) && meta::is_range_v<T>)>> {
  std::size_t operator()(T const& v) const { return hash_value(v); }
};

template <typename T>
auto value(const T& obj) {
  return _<T>{}(obj);
}

}  // namespace hash

}  // namespace sequant

#endif  // SEQUANT_HASH_HPP
