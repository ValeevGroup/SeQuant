//
// Created by Eduard Valeyev on 3/27/18.
//

#ifndef SEQUANT_CONTAINER_HPP
#define SEQUANT_CONTAINER_HPP

#include <SeQuant/core/hash.hpp>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/range.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>

#include <map>
#include <set>
#include <vector>

namespace sequant {

namespace container {

template <typename T>
using vector = std::vector<T>;
template <typename T, std::size_t N = 8>
using svector = boost::container::small_vector<T, N>;

template <typename Key, typename Compare = std::less<Key>>
using set = boost::container::flat_set<Key, Compare>;
template <typename Key, typename Value, typename Compare = std::less<Key>>
using map = boost::container::flat_map<Key, Value, Compare>;
template <typename Key, typename Value, typename Compare = std::less<Key>>
using multimap = boost::container::flat_multimap<Key, Value, Compare>;
template <typename Value, class Hash = sequant::hash::_<Value>,
          class EqualTo = std::equal_to<Value>>
using unordered_set = boost::unordered::unordered_set<Value, Hash, EqualTo>;
template <typename Key, typename Value, class Hash = sequant::hash::_<Key>,
          class EqualTo = std::equal_to<Key>>
using unordered_map =
    boost::unordered::unordered_map<Key, Value, Hash, EqualTo>;

using boost::begin;
using boost::end;
}  // namespace container

}  // namespace sequant

#endif  // SEQUANT_CONTAINER_HPP
