#ifndef SEQUANT_CORE_RATIONAL_H
#define SEQUANT_CORE_RATIONAL_H

#include <SeQuant/core/hash.hpp>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace boost {

template <typename T>
inline auto hash_value(const boost::rational<T>& i) {
  auto val = sequant::hash::value(i.numerator());
  sequant::hash::combine(val, i.denominator());
  return val;
}

template <typename T>
std::string to_string(const boost::rational<T>& i) {
  using std::to_string;
  return i.denominator() == 1
             ? to_string(i.numerator())
             : to_string(i.numerator()) + "/" + to_string(i.denominator());
}

template <typename T>
std::wstring to_wstring(const boost::rational<T>& i) {
  using std::to_wstring;
  return i.denominator() == 1
             ? to_wstring(i.numerator())
             : to_wstring(i.numerator()) + L"/" + to_wstring(i.denominator());
}

}  // namespace boost

namespace sequant {

// introducing multiprecision
using namespace boost::multiprecision;

// typedef boost::multiprecision::number<
//     cpp_int_backend<0, 0, signed_magnitude, checked, void>>
//     mp_int_type;

// try high fixed precision
typedef int512_t mp_int_type;

using rational = boost::rational<mp_int_type>;

// using rational = boost::rational<std::int64_t>;

/// shorter=sweeter? sometimes
using ratio = rational;

// clang-format off
///@note: boost::multiprecision::cpp_int only has limited support for constexpr
/// see: https://www.boost.org/doc/libs/1_82_0/libs/multiprecision/doc/html/boost_multiprecision/tut/lits.html
// clang-format on

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
inline rational to_rational(
    T t, T eps = std::sqrt(std::numeric_limits<T>::epsilon()),
    std::size_t max_niter = 1000) {
  if (std::isnan(t) || std::isinf(t)) {
    throw std::invalid_argument(
        "sequant::to_rational: cannot make a rational out of " +
        std::to_string(t));
  }
  // e.g.
  // https://gist.github.com/mikeando/7073d62385a34a61a6f7#file-main2-cpp-L42
  auto sbtree = [max_niter](T f, T tol) {
    using fraction = rational;
    using Int = typename rational::int_type;
    Int base = std::floor(f);
    f -= base;
    if (f < tol) return fraction(base, 1);
    if (1 - tol < f) return fraction(base + 1, 1);

    fraction lower(0, 1);
    fraction upper(1, 1);

    std::size_t niter = 0;
    while (niter < max_niter) {
      fraction middle(lower.numerator() + upper.numerator(),
                      lower.denominator() + upper.denominator());

      if (middle.denominator() * (f + tol) < middle.numerator()) {
        upper = middle;
      } else if (middle.numerator() < middle.denominator() * (f - tol)) {
        lower = middle;
      } else {
        return fraction(middle.denominator() * base + middle.numerator(),
                        middle.denominator());
      }
      ++niter;
    }
    throw std::invalid_argument(
        "sequant::rationalize: could not rationalize " + std::to_string(f) +
        " to eps=" + std::to_string(tol) + " in " + std::to_string(max_niter) +
        " iterations");  // unreachable
  };
  return sbtree(t, eps);
}

template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
rational to_rational(T t) {
  return rational{t};
}

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T> ||
                                                  std::is_integral_v<T>>>
inline ratio to_ratio(T t) {
  return to_rational(t);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_RATIONAL_H
