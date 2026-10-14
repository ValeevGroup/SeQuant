//
// Created by Eduard Valeyev on 5/17/23.
//

#ifndef SEQUANT_CORE_MATH_HPP
#define SEQUANT_CORE_MATH_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <utility>

namespace sequant {

constexpr std::size_t pow2(std::size_t n) {
  SEQUANT_ASSERT(n <= 63);
  return 1ul << n;
}

/// Returns the factorial of `n`

/// @param n the argument; there is no practical limit on the supported
/// magnitude of `n`
/// @return `n!`
/// @internal uses pretabulated factorials for up to `n=20`, and memoized values
/// for `n>20`.
/// @note this function is reentrant
inline sequant::intmax_t factorial(std::size_t n) {
  using result_t = sequant::intmax_t;
  constexpr std::size_t n_max_precomputed =
      20;  // 20! is the largest factorial that fits into uint64_t
  const result_t values[n_max_precomputed + 1] = {1ull,
                                                  1ull,
                                                  2ull,
                                                  6ull,
                                                  24ull,
                                                  120ull,
                                                  720ull,
                                                  5040ull,
                                                  40320ull,
                                                  362880ull,
                                                  3628800ull,
                                                  39916800ull,
                                                  479001600ull,
                                                  6227020800ull,
                                                  87178291200ull,
                                                  1307674368000ull,
                                                  20922789888000ull,
                                                  355687428096000ull,
                                                  6402373705728000ull,
                                                  121645100408832000ull,
                                                  2432902008176640000ull};
  if (n <= n_max_precomputed)
    return values[n];
  else {
    // this memoizes values for n>n_max_precomputed
#if defined(__cpp_lib_atomic_shared_ptr) && \
    __cpp_lib_atomic_shared_ptr >= 201711L
    using ValHolder =
        std::atomic<std::shared_ptr<std::vector<sequant::intmax_t>>>;
    auto load = [](auto&& holder) { return holder.load(); };
    auto store = [](auto&& holder, ValHolder::value_type arg) {
      holder.store(std::move(arg));
    };
#else
    using ValHolder = std::shared_ptr<std::vector<sequant::intmax_t>>;
    auto load = [](auto&& holder) { return std::atomic_load(&holder); };
    auto store = [](auto&& holder, ValHolder arg) {
      std::atomic_store(&holder, std::move(arg));
    };
#endif
    static ValHolder memvals;
    // used to serialize access to memvals
    static std::mutex memvals_mutex;

    const std::size_t n_in_memvals = n - n_max_precomputed - 1;
    const auto memvals_handle = load(memvals);
    if (memvals_handle && memvals_handle->size() > n_in_memvals) {
      return (*memvals_handle)[n_in_memvals];
    } else {
      std::lock_guard<std::mutex> lock(memvals_mutex);
      // memvals might have been updated while we were waiting for the lock
      // hence reload the state
      const auto memvals_handle = load(memvals);
      if (memvals_handle && memvals_handle->size() > n_in_memvals) {
        return (*memvals_handle)[n_in_memvals];
      } else {
        auto new_memvals = std::make_shared<std::vector<sequant::intmax_t>>(
            memvals_handle
                ? *memvals_handle
                : std::vector<sequant::intmax_t>{});  // copy old memvals
        new_memvals->reserve(n_in_memvals + 1);
        for (auto i = new_memvals->size(); i <= n_in_memvals; ++i)
          new_memvals->emplace_back(
              (i == 0 ? values[n_max_precomputed] : new_memvals->back()) *
              (i + n_max_precomputed + 1));
        store(memvals, new_memvals);
        return (*new_memvals)[n_in_memvals];
      }
    }
  }
}

/// @param[in] n a non-negative integer
/// @param[in] k the order of the root, at least 1
/// @return the floor of the @p k -th root of @p n
inline intmax_t iroot(const intmax_t& n, unsigned k) {
  SEQUANT_ASSERT(n >= 0);
  SEQUANT_ASSERT(k >= 1);
  if (k == 1 || n < 2) return n;
  // Newton's iteration from above converges monotonically to the floor, but
  // only by a factor of 1 - 1/k per step until it is close, so it starts
  // from a floating-point estimate rounded up by more than its error
  const unsigned bits = msb(n) + 1;
  const unsigned shift = bits > 53 ? bits - 53 : 0;
  const double log2_root =
      (std::log2(intmax_t(n >> shift).convert_to<double>()) + shift) / k + 1e-9;
  const auto root_shift = static_cast<unsigned>(log2_root);
  intmax_t x{static_cast<unsigned long long>(
      std::ceil(std::exp2(log2_root - root_shift) * (1ull << 52)))};
  if (root_shift >= 52)
    x <<= root_shift - 52;
  else
    x >>= 52 - root_shift;
  ++x;
  while (true) {
    const intmax_t y = ((k - 1) * x + n / pow(x, k - 1)) / k;
    if (y >= x) break;
    x = y;
  }
  // in case the estimate was below the root after all
  while (pow(x + 1, k) <= n) ++x;
  return x;
}

/// splits positive integers @p numbers into a coprime base: pairwise coprime
/// integers greater than 1 of which every number is a product of powers.
/// Needs only gcds, unlike a prime factorization.
inline container::svector<intmax_t> coprime_base(
    const container::svector<intmax_t>& numbers) {
  container::svector<intmax_t> base;
  for (const auto& n : numbers)
    if (n > 1) base.push_back(n);
  // replacing a, b by a/g, b/g, g for g = gcd(a, b) > 1 lowers the product
  // of the elements, so this terminates
  for (bool refined = true; refined;) {
    refined = false;
    for (std::size_t i = 0; i < base.size(); ++i) {
      for (std::size_t j = i + 1; j < base.size(); ++j) {
        const intmax_t g = gcd(base[i], base[j]);
        if (g == 1) continue;
        base[i] /= g;
        base[j] /= g;
        base.push_back(g);
        refined = true;
      }
    }
    base.erase(std::remove(base.begin(), base.end(), intmax_t(1)), base.end());
  }
  return base;
}

/// the simplest radical form of a positive real number,
/// `rational_part * radicand^(1/root_order)`
struct RadicalForm {
  rational rational_part;
  intmax_t radicand;    ///< positive integer; 1 iff root_order is 1
  intmax_t root_order;  ///< positive integer
};

/// the default bound, in bits, on the size of a scalar built by folding
/// Powers (see simplest_radical_form() and Power::flatten())
inline constexpr unsigned default_max_folded_bits = 1u << 16;

/// @brief computes the simplest radical form of a product of rational powers
/// of positive rationals
///
/// The bases are split by coprime_base(), each element with a non-integer
/// exponent is factored by trial division up to @p max_trial_divisor (one
/// with an integer exponent contributes only to the rational part and is
/// left whole), and the exponents of each prime, and of the cofactor trial
/// division leaves, are summed; a cofactor that is a perfect power is
/// replaced by its root first. The integer part of each sum
/// goes into the rational part, the fractional part into the radicand, and
/// the root order is the lcm of the denominators of the fractional parts. The
/// radicand is thus not divisible by the `root_order`-th power of a prime up
/// to @p max_trial_divisor, nor of a larger prime whose power is the whole
/// cofactor of an element, so the result depends only on the value of the
/// product unless a cofactor has two or more distinct prime factors above
/// @p max_trial_divisor and is not a perfect power.
/// @param[in] powers (base, exponent) pairs; every base must be positive
/// @throw Exception if a base is not positive
/// @param[in] max_trial_divisor the largest divisor tried when factoring
/// @param[in] max_bits the largest size, in bits, of the radicand and of the
/// numerator and the denominator of the rational part; a larger result is
/// refused rather than built (exponents with large coprime denominators
/// produce one, since their lcm multiplies every exponent); a base whose
/// integer part of the exponent alone would exceed it is refused before it
/// is factored
/// @return the form, or null if it would exceed @p max_bits
/// @note the search for the roots of a cofactor costs superlinearly in its
/// size: a base of 20000 bits with no factor below @p max_trial_divisor and
/// a non-integer exponent takes about a second
inline std::optional<RadicalForm> simplest_radical_form(
    const container::svector<std::pair<rational, rational>>& powers,
    unsigned max_trial_divisor = 1000,
    unsigned max_bits = default_max_folded_bits) {
  // (positive integer, exponent) terms: a numerator carries the exponent of
  // its base, a denominator its negative
  container::svector<std::pair<intmax_t, rational>> terms;
  container::svector<intmax_t> numbers;
  for (const auto& [base, exponent] : powers) {
    if (base <= 0)
      throw Exception("simplest_radical_form: every base must be positive");
    terms.emplace_back(numerator(base), exponent);
    terms.emplace_back(denominator(base), -exponent);
    numbers.push_back(numerator(base));
    numbers.push_back(denominator(base));
  }

  // prime or cofactor -> total exponent
  std::map<intmax_t, rational> exponents;
  for (const auto& element : coprime_base(numbers)) {
    rational exponent{0};
    for (const auto& [n, e] : terms) {
      intmax_t rest = n;
      unsigned multiplicity = 0;
      for (; rest % element == 0; rest /= element) ++multiplicity;
      exponent += multiplicity * e;
    }
    if (exponent == 0) continue;
    // whatever root of the element is taken, the integer part of the
    // exponent alone puts at least floor(|exponent|) times the width of the
    // element into the rational part, so a hopeless element is refused
    // before the (costly) factoring
    const intmax_t whole = abs(numerator(exponent)) / denominator(exponent);
    if (whole * (msb(element) + 1) > max_bits) return std::nullopt;
    // an integer exponent leaves nothing for the radicand, so the element
    // goes whole into the rational part without being factored
    if (denominator(exponent) == 1) {
      exponents[element] += exponent;
      continue;
    }

    intmax_t cofactor = element;
    intmax_t d = 2;
    for (; d <= max_trial_divisor && d * d <= cofactor; d += (d == 2 ? 1 : 2)) {
      unsigned multiplicity = 0;
      for (; cofactor % d == 0; cofactor /= d) ++multiplicity;
      if (multiplicity != 0) exponents[d] += multiplicity * exponent;
    }
    if (cofactor == 1) continue;
    // if trial division ran out of candidates below the square root the
    // cofactor is a prime; otherwise it has only prime factors above
    // max_trial_divisor, and an exact k-th root exposes its base, and with
    // it the integer part of the exponent. Only prime k need trying, since
    // a composite k is covered by its prime factors, and the root exceeds
    // max_trial_divisor, so msb(cofactor) >= k * msb(max_trial_divisor)
    if (d * d <= cofactor) {
      const unsigned root_bits =
          max_trial_divisor < 2 ? 1 : std::bit_width(max_trial_divisor) - 1;
      auto is_prime = [](unsigned k) {
        for (unsigned f = 2; f * f <= k; ++f)
          if (k % f == 0) return false;
        return true;
      };
      for (unsigned k = 2; k <= msb(cofactor) / root_bits; ++k) {
        if (!is_prime(k)) continue;
        for (intmax_t r = iroot(cofactor, k); pow(r, k) == cofactor;
             r = iroot(cofactor, k)) {
          cofactor = r;
          exponent *= k;
        }
      }
    }
    exponents[cofactor] += exponent;
  }

  RadicalForm result{rational{1}, intmax_t{1}, intmax_t{1}};
  // upper bound on the bits of base^n; keeping each part of the result under
  // max_bits also keeps n within the unsigned exponent of pow
  auto bits = [](const intmax_t& base, const intmax_t& n) {
    return n * (msb(base) + 1);
  };
  intmax_t numerator_bits{0}, denominator_bits{0}, radicand_bits{0};
  container::svector<std::pair<intmax_t, rational>> radicand_terms;
  for (const auto& [base, exponent] : exponents) {
    intmax_t integer_part = numerator(exponent) / denominator(exponent);
    if (exponent < integer_part) --integer_part;  // floor
    if (integer_part > 0) {
      numerator_bits += bits(base, integer_part);
      if (numerator_bits > max_bits) return std::nullopt;
      result.rational_part *= pow(base, integer_part.convert_to<unsigned>());
    } else if (integer_part < 0) {
      denominator_bits += bits(base, -integer_part);
      if (denominator_bits > max_bits) return std::nullopt;
      result.rational_part /= pow(base, (-integer_part).convert_to<unsigned>());
    }
    const rational fractional_part = exponent - integer_part;
    if (fractional_part != 0) {
      radicand_terms.emplace_back(base, fractional_part);
      result.root_order = lcm(result.root_order, denominator(fractional_part));
    }
  }
  for (const auto& [base, exponent] : radicand_terms) {
    const intmax_t m =
        numerator(exponent) * (result.root_order / denominator(exponent));
    radicand_bits += bits(base, m);
    if (radicand_bits > max_bits) return std::nullopt;
    result.radicand *= pow(base, m.convert_to<unsigned>());
  }
  return result;
}

}  // namespace sequant

#endif  // SEQUANT_CORE_MATH_HPP
