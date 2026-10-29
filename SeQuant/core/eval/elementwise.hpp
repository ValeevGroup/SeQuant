#ifndef SEQUANT_CORE_EVAL_ELEMENTWISE_HPP
#define SEQUANT_CORE_EVAL_ELEMENTWISE_HPP

#include <SeQuant/core/rational.hpp>

#include <cmath>
#include <complex>  // the complex overloads of std::pow
#include <cstdint>
#include <limits>

namespace sequant::detail {

///
/// \brief The per-element kernel of \c Result::pow: `x -> x^e` for a fixed
///        rational exponent \c e.
///
/// The exponent is decoded once, at construction, so applying the kernel to
/// each element costs no multiprecision arithmetic. An integer exponent is
/// exact: repeated multiplication (by squaring), then the reciprocal for a
/// negative exponent, so `0^-n` is `inf` per IEEE and is not checked. Any other
/// exponent goes through \c std::pow with the exponent converted to \c double
/// (the complex overload for complex elements).
///
class ElementwisePow {
 public:
  explicit ElementwisePow(rational const& exponent) {
    auto const num = numerator(exponent);
    if (denominator(exponent) == 1 &&
        abs(num) <= std::numeric_limits<std::uint64_t>::max()) {
      integral_ = true;
      negative_ = num < 0;
      n_ = abs(num).convert_to<std::uint64_t>();
    } else {
      exponent_ = static_cast<double>(exponent);
    }
  }

  template <typename T>
  [[nodiscard]] T operator()(T const& x) const {
    if (!integral_) return static_cast<T>(std::pow(x, exponent_));
    T result{1};
    T base = x;
    for (auto n = n_; n != 0;) {
      if (n & 1) result *= base;
      n >>= 1;
      if (n != 0) base *= base;
    }
    return negative_ ? T{1} / result : result;
  }

 private:
  bool integral_ = false;
  bool negative_ = false;
  std::uint64_t n_ = 0;
  double exponent_ = 0;
};

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_EVAL_ELEMENTWISE_HPP
