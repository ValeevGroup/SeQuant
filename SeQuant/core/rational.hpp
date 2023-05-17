#ifndef SEQUANT_CORE_RATIONAL_H
#define SEQUANT_CORE_RATIONAL_H

#include <SeQuant/core/hash.hpp>

#include <boost/rational.hpp>

namespace boost {

template <typename T>
inline auto hash_value(const boost::rational<T>& i) {
  auto val = sequant::hash::value(i.numerator());
  sequant::hash::combine(val, i.denominator());
  return val;
}

}  // namespace boost

namespace sequant {
using rational = boost::rational<std::int64_t>;
}

#endif  // SEQUANT_CORE_RATIONAL_H
