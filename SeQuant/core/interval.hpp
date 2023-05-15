#ifndef SEQUANT_EXTERNAL_BOOST_INTERVAL_H
#define SEQUANT_EXTERNAL_BOOST_INTERVAL_H

#include <SeQuant/core/hash.hpp>

// boost/numeric/interval does not know about arm rounding .. on arm64/macos use
// c99 rounding
#if defined(__arm64__) && defined(__APPLE__) && !defined(__USE_ISOC99)
#define __USE_ISOC99 1
#include <boost/numeric/interval.hpp>
#undef __USE_ISOC99
#else
#include <boost/numeric/interval.hpp>
#endif

namespace boost::numeric {

template <typename T>
inline auto hash_value(const boost::numeric::interval<T>& i) {
  auto val = sequant::hash::value(i.lower());
  sequant::hash::combine(val, i.upper());
  return val;
}

}  // namespace boost::numeric

#endif  // SEQUANT_EXTERNAL_BOOST_INTERVAL_H
