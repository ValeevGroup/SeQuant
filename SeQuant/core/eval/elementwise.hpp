#ifndef SEQUANT_CORE_EVAL_ELEMENTWISE_HPP
#define SEQUANT_CORE_EVAL_ELEMENTWISE_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cmath>
#include <complex>  // the complex overloads of std::pow
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>

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

///
/// \brief The operand mode of each mode of a Result::broadcast result.
///
/// \param operand the operand's mode labels
/// \param result the result's mode labels, a superset of \p operand's
/// \return for each mode of \p result, its position in \p operand, or
///         \c std::nullopt for a mode the broadcast adds
/// \throw Exception if a mode of \p operand is not a mode of \p result
///
template <typename Labels>
[[nodiscard]] container::svector<std::optional<std::size_t>>
broadcast_operand_modes(Labels const& operand, Labels const& result) {
  container::svector<std::optional<std::size_t>> pos;
  pos.reserve(result.size());
  std::size_t nshared = 0;
  for (auto const& lbl : result) {
    pos.push_back(find_position(operand, lbl));
    if (pos.back()) ++nshared;
  }
  if (nshared != operand.size())
    throw Exception(
        "broadcast: a mode of the operand is not a mode of the result");
  return pos;
}

///
/// \brief The per-element kernel of \c Result::broadcast: copy a strided
///        operand into a strided result of \p extents, where a zero operand
///        stride repeats the operand along that mode.
///
/// `out[sum_p r_p out_strides[p]] = in[sum_p r_p in_strides[p]]` for every
/// multi-index `r` within \p extents.
///
template <typename T>
void broadcast_copy(T* out, container::svector<std::size_t> const& extents,
                    container::svector<std::size_t> const& out_strides,
                    T const* in,
                    container::svector<std::size_t> const& in_strides) {
  auto const rank = extents.size();
  SEQUANT_ASSERT(out_strides.size() == rank && in_strides.size() == rank);
  if (rank == 0) {
    *out = *in;
    return;
  }
  for (auto const e : extents)
    if (e == 0) return;
  auto const last = rank - 1;
  container::svector<std::size_t> idx(rank, 0);
  std::size_t o = 0, i = 0;
  while (true) {
    for (std::size_t k = 0; k < extents[last]; ++k)
      out[o + k * out_strides[last]] = in[i + k * in_strides[last]];
    // advance the multi-index over the modes before the last
    std::size_t d = last;
    while (true) {
      if (d == 0) return;
      --d;
      ++idx[d];
      o += out_strides[d];
      i += in_strides[d];
      if (idx[d] < extents[d]) break;
      o -= idx[d] * out_strides[d];
      i -= idx[d] * in_strides[d];
      idx[d] = 0;
    }
  }
}

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_EVAL_ELEMENTWISE_HPP
