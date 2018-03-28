//
// Created by Eduard Valeyev on 3/27/18.
//

#ifndef SEQUANT2_VECTOR_HPP
#define SEQUANT2_VECTOR_HPP

#include <boost/container/small_vector.hpp>
#include <boost/range.hpp>

namespace sequant2 {

namespace container {

template<typename T> using vector = std::vector<T>;
//template<typename T, std::size_t N = 8> using svector = boost::container::small_vector<T, N>;
template<typename T> using svector = std::vector<T>;

using boost::begin;
using boost::end;
}  // namespace sequant2::container

}  // namespace sequant2

#endif //SEQUANT2_VECTOR_HPP
