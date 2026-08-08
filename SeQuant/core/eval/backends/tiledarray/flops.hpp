#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_FLOPS_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_FLOPS_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

/// \file backends/tiledarray/flops.hpp
/// Measurement of the real floating-point work of a TiledArray backend
/// operation; feeds sequant::eval::FlopCounter.
///
/// Every extent used here is read off a REAL object: a local tile's
/// `Range`, or -- for a tensor of tensors -- a real inner tensor's `Range`.
/// Nominal index sizes are never consulted, which is the whole point: for a
/// CSV/PNO tensor of tensors the per-pair inner dimensions vary, and the dense
/// extent product is exactly the model being tested.
///
/// Everything here is LOCAL to the calling rank (only local, non-zero tiles
/// are visited) and performs no collective.
///
/// Every entry point is called only from behind
/// `if (eval::FlopCounter::enabled())`, so nothing in this file is on the
/// disabled path.

#include <SeQuant/core/eval/flops.hpp>

#include <tiledarray.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::detail {

/// An annotation split into its outer and inner label lists. SeQuant
/// annotations are comma-separated labels with an optional ';' separating the
/// outer (non-proto-indexed) labels from the inner (proto-indexed) ones; see
/// EvalExpr::indices_annot.
struct AnnotLabels {
  std::vector<std::string_view> outer, inner;
};

[[nodiscard]] inline AnnotLabels split_annot(std::string const& annot) {
  AnnotLabels out;
  std::string_view s{annot};
  auto const semi = s.find(';');
  auto push = [](std::vector<std::string_view>& dst, std::string_view part) {
    while (!part.empty()) {
      auto const comma = part.find(',');
      auto const tok = part.substr(0, comma);
      if (!tok.empty()) dst.push_back(tok);
      if (comma == std::string_view::npos) break;
      part.remove_prefix(comma + 1);
    }
  };
  if (semi == std::string_view::npos) {
    push(out.outer, s);
  } else {
    push(out.outer, s.substr(0, semi));
    push(out.inner, s.substr(semi + 1));
  }
  return out;
}

[[nodiscard]] inline bool contains(std::vector<std::string_view> const& v,
                                   std::string_view x) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

/// Real element count of \p arr's LOCAL, non-zero tiles: the sum of each
/// tile's actual `Range` volume, and for a tensor of tensors the sum of each
/// non-empty inner tensor's actual volume. This is the measured footprint,
/// not the dense extent product.
template <typename ArrayT>
[[nodiscard]] double real_element_count(ArrayT const& arr) {
  using value_type = typename ArrayT::value_type;
  double n = 0;
  if (!arr.is_initialized()) return 0;
  for (auto it = arr.begin(); it != arr.end(); ++it) {
    // own the tile by value: TileReference::get() returns by value and binding
    // to a reference would dangle (see tot_inner_rank's note).
    auto const tile = it->get();
    if constexpr (TA::detail::is_tensor_of_tensor_v<value_type>) {
      for (auto const& inner : tile)
        if (!inner.empty()) n += static_cast<double>(inner.range().volume());
    } else {
      n += static_cast<double>(tile.range().volume());
    }
  }
  return n;
}

/// Nominal extent of the mode of \p arr labelled \p label under \p labels,
/// read from the array's own TiledRange (i.e. the real dimension of the real
/// array, not an Index's approximate_size). 0 if the label is absent.
template <typename ArrayT>
[[nodiscard]] double mode_extent(ArrayT const& arr,
                                 std::vector<std::string_view> const& labels,
                                 std::string_view label) {
  for (std::size_t d = 0; d < labels.size(); ++d)
    if (labels[d] == label)
      return static_cast<double>(arr.trange().dim(d).extent());
  return 0;
}

/// Record a whole-array elementwise operation (add, permute, scale): one
/// arithmetic operation (or, for a permute, one move) per element actually
/// written.
template <typename ArrayT>
inline void count_elementwise(eval::FlopCategory cat, ArrayT const& result) {
  double const v = real_element_count(result);
  eval::FlopCounter::record(cat, cat == eval::FlopCategory::Permute ? 0.0 : v,
                            v);
}

/// Record a binary product whose result is a flat array, or whose contracted
/// indices are all flat/outer: `2 * (real result elements) * K`, with K the
/// product of the contracted modes' real extents. A product with no contracted
/// index is a Hadamard product, one multiply per element.
template <typename LeftT, typename RightT, typename ResultT>
inline void count_product_by_extent(LeftT const& l, RightT const& r,
                                    ResultT const& c,
                                    std::string const& lannot,
                                    std::string const& rannot,
                                    std::string const& cannot) {
  auto const la = split_annot(lannot), ra = split_annot(rannot),
             ca = split_annot(cannot);
  // Contracted labels may sit in either the outer or the inner list of an
  // operand (a flat operand contracted against a ToT's inner mode); pool them.
  auto pool = [](AnnotLabels const& a) {
    std::vector<std::string_view> v = a.outer;
    v.insert(v.end(), a.inner.begin(), a.inner.end());
    return v;
  };
  auto const lp = pool(la), rp = pool(ra), cp = pool(ca);
  std::size_t n_contracted = 0;
  // extents come from whichever operand is flat enough to have them in its
  // TiledRange; for a ToT the outer modes are in the TiledRange too.
  double k = 1;
  for (auto lab : lp)
    if (contains(rp, lab) && !contains(cp, lab)) {
      double e = mode_extent(l, la.outer, lab);
      if (e == 0) e = mode_extent(r, ra.outer, lab);
      if (e == 0) {  // an inner mode of a ToT: extent not in any TiledRange
        eval::FlopCounter::record_unsourced();
        return;
      }
      k *= e;
      ++n_contracted;
    }
  double const v = real_element_count(c);
  if (n_contracted == 0) {
    eval::FlopCounter::record(eval::FlopCategory::Contraction, v, v);
    eval::FlopCounter::record_kernels(1, v, 1, 1);
    return;
  }
  eval::FlopCounter::record(eval::FlopCategory::Contraction, 2 * v * k, v);
  eval::FlopCounter::record_kernels(1, v, 1, k);
}

/// Record a tensor-of-tensor * tensor-of-tensor product that keeps its nest
/// (`ToT * ToT -> ToT`) by walking the result's real inner tensors.
///
/// For each outer element c the work is one GEMM of shape (M(c), N(c), K(c)):
///   |inner_C(c)| = M(c)*N(c)*H(c),  |inner_L(c)| = M(c)*K(c)*H(c)
/// with H the inner Hadamard extent, so with
///   M' := prod of the extents (read off inner_C(c)) of the labels shared by
///         inner_L and inner_C   ( = M(c)*H(c) )
///   K(c) = |inner_L(c)| / M'
/// the flop count is `2 * |inner_C(c)| * K(c)` -- exact, and derived entirely
/// from the two real inner tensors at c. This is where the per-pair variation
/// of a CSV/PNO array is captured.
///
/// Requires the outer indices to be free/Hadamard (no outer contraction), so
/// that each result outer element has exactly one matching left outer element.
/// Anything else, or a left operand whose outer element is not local, is
/// recorded as unsourced rather than guessed.
template <typename LeftT, typename ResultT>
inline void count_product_tot(LeftT const& l, ResultT const& c,
                              std::string const& lannot,
                              std::string const& rannot,
                              std::string const& cannot) {
  auto const la = split_annot(lannot), ra = split_annot(rannot),
             ca = split_annot(cannot);
  // outer contraction (a label in both operands' outer lists but not the
  // result's) breaks the one-to-one outer correspondence this walk needs
  for (auto lab : la.outer)
    if (contains(ra.outer, lab) && !contains(ca.outer, lab)) {
      eval::FlopCounter::record_unsourced();
      return;
    }
  // position of each of the left operand's outer labels within the result's
  // outer labels
  std::vector<std::size_t> l_from_c(la.outer.size());
  for (std::size_t d = 0; d < la.outer.size(); ++d) {
    auto const it =
        std::find(ca.outer.begin(), ca.outer.end(), la.outer[d]);
    if (it == ca.outer.end()) {
      eval::FlopCounter::record_unsourced();
      return;
    }
    l_from_c[d] = static_cast<std::size_t>(it - ca.outer.begin());
  }

  double flops = 0, elems = 0, kernels = 0, m_sum = 0, n_sum = 0, k_sum = 0;
  using index1_type = typename ResultT::trange_type::range_type::index1_type;
  std::vector<index1_type> lidx(la.outer.size());
  for (auto it = c.begin(); it != c.end(); ++it) {
    auto const ctile = it->get();
    for (auto const& cidx : ctile.range()) {
      auto const& inner_c = ctile[cidx];
      if (inner_c.empty()) continue;
      double const vol_c = static_cast<double>(inner_c.range().volume());
      elems += vol_c;
      for (std::size_t d = 0; d < l_from_c.size(); ++d)
        lidx[d] = static_cast<index1_type>(cidx[l_from_c[d]]);
      auto const ltile_idx = l.trange().element_to_tile(lidx);
      if (l.is_zero(ltile_idx)) continue;  // zero operand: no work
      if (!l.is_local(ltile_idx)) {
        eval::FlopCounter::record_unsourced();
        return;
      }
      auto const& inner_l = l.find_local(ltile_idx).get()[lidx];
      if (inner_l.empty()) continue;
      // M' = product of the extents (from inner_C) of the labels inner_L and
      // inner_C share; those are the external-left plus inner-Hadamard modes
      double mprime = 1;
      for (std::size_t d = 0; d < ca.inner.size(); ++d)
        if (contains(la.inner, ca.inner[d]))
          mprime *= static_cast<double>(inner_c.range().extent(d));
      if (mprime <= 0) continue;
      double const kk = static_cast<double>(inner_l.range().volume()) / mprime;
      flops += 2 * vol_c * kk;
      ++kernels;
      m_sum += mprime;
      n_sum += vol_c / mprime;
      k_sum += kk;
    }
  }
  eval::FlopCounter::record(eval::FlopCategory::Contraction, flops, elems);
  eval::FlopCounter::record_kernels(kernels, m_sum, n_sum, k_sum);
}

/// Record a `ToT * ToT -> T` (de-nesting) product: every inner mode is
/// contracted, so each outer element contributes one dot product of length
/// `|inner_L|`. Exact when the outer indices are pure Hadamard (identical
/// label sets on both operands and the result), which is the only shape the
/// CSV path emits; otherwise recorded as unsourced.
template <typename LeftT>
inline void count_product_denest(LeftT const& l, std::string const& lannot,
                                 std::string const& rannot,
                                 std::string const& cannot) {
  auto const la = split_annot(lannot), ra = split_annot(rannot),
             ca = split_annot(cannot);
  auto same_set = [](std::vector<std::string_view> a,
                     std::vector<std::string_view> b) {
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return a == b;
  };
  if (!same_set(la.outer, ra.outer) || !same_set(la.outer, ca.outer)) {
    eval::FlopCounter::record_unsourced();
    return;
  }
  double flops = 0, kernels = 0, k_sum = 0;
  for (auto it = l.begin(); it != l.end(); ++it) {
    auto const tile = it->get();
    for (auto const& inner : tile) {
      if (inner.empty()) continue;
      double const v = static_cast<double>(inner.range().volume());
      flops += 2 * v;
      ++kernels;
      k_sum += v;
    }
  }
  eval::FlopCounter::record(eval::FlopCategory::Contraction, flops, kernels);
  eval::FlopCounter::record_kernels(kernels, kernels, kernels, k_sum);
}

}  // namespace sequant::detail

#endif  // SEQUANT_HAS_TILEDARRAY
#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_FLOPS_HPP
