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
#include <limits>
#include <string>
#include <utility>
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
///
/// LIMITATION, and the only one left in this file: the result side is exact
/// (only local, non-zero tiles are counted, at their real volumes), but K is
/// the operand's full dimension. If an operand is block-sparse ALONG A
/// CONTRACTED MODE, the skipped blocks are still charged. This is an upper
/// bound on those products, never an under-count. The tensor-of-tensor walks
/// below have no such gap: they read K off a real inner tensor.
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
        eval::FlopCounter::record_unsourced(
            "contracted inner mode has no TiledRange extent: " + lannot +
            " * " + rannot + " = " + cannot);
        return;
      }
      k *= e;
      ++n_contracted;
    }
  double const v = real_element_count(c);
  // No kernel shape is recorded here: TiledArray splits a flat contraction
  // into many per-tile GEMMs whose individual shapes this layer never sees, so
  // adding a fictitious single (M,N,K) would corrupt the shape statistics that
  // the tensor-of-tensor walks measure exactly.
  if (n_contracted == 0) {
    eval::FlopCounter::record(eval::FlopCategory::Contraction, v, v);
    return;
  }
  eval::FlopCounter::record(eval::FlopCategory::Contraction, 2 * v * k, v);
}

/// Total real inner volume of a tensor of tensors: the sum of every non-empty
/// local inner tensor's actual volume.
template <typename ArrayT>
[[nodiscard]] inline double inner_volume(ArrayT const& arr) {
  return real_element_count(arr);
}

/// Extents, in \p labels order, of a real inner tensor.
template <typename InnerT>
[[nodiscard]] inline double shared_extent_product(
    InnerT const& inner, std::vector<std::string_view> const& labels,
    std::vector<std::string_view> const& against) {
  double p = 1;
  for (std::size_t d = 0; d < labels.size(); ++d)
    if (contains(against, labels[d]))
      p *= static_cast<double>(inner.range().extent(d));
  return p;
}

/// Accumulator for one tensor-of-tensor product's measured work.
struct TotWork {
  double flops = 0, kernels = 0, m_sum = 0, n_sum = 0, k_sum = 0;
  /// One GEMM of shape (M', |inner_C|/M', K): see count_product_tot.
  void gemm(double vol_c, double mprime, double vol_operand) {
    double const kk = vol_operand / mprime;
    flops += 2 * vol_c * kk;
    ++kernels;
    m_sum += mprime;
    n_sum += vol_c / mprime;
    k_sum += kk;
  }
};

/// Walk the operand \p o (left or right -- the roles are symmetric) whose
/// outer labels cover the result's, one GEMM per real outer element.
/// \return false, having recorded the reason, if a needed result element is
///         not local.
template <typename OperandT, typename ResultT>
[[nodiscard]] inline bool count_tot_walk_operand(
    OperandT const& o, ResultT const& c, AnnotLabels const& oa,
    AnnotLabels const& ca, TotWork& w, std::string const& why) {
  std::vector<std::size_t> c_from_o(ca.outer.size());
  for (std::size_t d = 0; d < ca.outer.size(); ++d)
    c_from_o[d] = static_cast<std::size_t>(
        std::find(oa.outer.begin(), oa.outer.end(), ca.outer[d]) -
        oa.outer.begin());
  using index1_type = typename ResultT::trange_type::range_type::index1_type;
  std::vector<index1_type> cidx(ca.outer.size());
  // The contracted outer index is typically the fastest-varying mode, so a
  // one-entry memo on the result element turns the result lookup into a hit
  // for a whole run of operand elements.
  std::vector<index1_type> memo_idx;
  double memo_vol = 0, memo_mprime = 0;
  bool memo_ok = false;
  for (auto it = o.begin(); it != o.end(); ++it) {
    auto const otile = it->get();
    for (auto const& oidx : otile.range()) {
      auto const& inner_o = otile[oidx];
      if (inner_o.empty()) continue;
      for (std::size_t d = 0; d < c_from_o.size(); ++d)
        cidx[d] = static_cast<index1_type>(oidx[c_from_o[d]]);
      if (!memo_ok || memo_idx != cidx) {
        memo_ok = false;
        auto const ctile_idx = c.trange().element_to_tile(cidx);
        if (c.is_zero(ctile_idx)) continue;  // zero result: no work
        if (!c.is_local(ctile_idx)) {
          eval::FlopCounter::record_unsourced("result outer tile is remote: " +
                                              why);
          return false;
        }
        auto const& inner_c = c.find_local(ctile_idx).get()[cidx];
        if (inner_c.empty()) continue;
        memo_idx = cidx;
        memo_vol = static_cast<double>(inner_c.range().volume());
        memo_mprime = shared_extent_product(inner_c, ca.inner, oa.inner);
        memo_ok = memo_mprime > 0;
        if (!memo_ok) continue;
      }
      w.gemm(memo_vol, memo_mprime,
             static_cast<double>(inner_o.range().volume()));
    }
  }
  return true;
}

/// Walk the RESULT's real inner tensors, and for each the full range of the
/// contracted OUTER indices (an odometer over the left operand's own
/// TiledRanges). Handles every outer structure, at a cost proportional to the
/// contracted outer volume -- which is why the operand walks above are
/// preferred when they apply.
template <typename LeftT, typename ResultT>
[[nodiscard]] inline bool count_tot_walk_result(
    LeftT const& l, ResultT const& c, AnnotLabels const& la,
    AnnotLabels const& ra, AnnotLabels const& ca, TotWork& w,
    std::string const& why) {
  using index1_type = typename ResultT::trange_type::range_type::index1_type;
  constexpr auto npos = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> l_from_c(la.outer.size(), npos);
  std::vector<std::size_t> contracted_dims;
  for (std::size_t d = 0; d < la.outer.size(); ++d) {
    auto const it = std::find(ca.outer.begin(), ca.outer.end(), la.outer[d]);
    if (it != ca.outer.end())
      l_from_c[d] = static_cast<std::size_t>(it - ca.outer.begin());
    else if (contains(ra.outer, la.outer[d]))
      contracted_dims.push_back(d);
    else {
      eval::FlopCounter::record_unsourced(
          "left outer label is neither free nor contracted: " + why);
      return false;
    }
  }
  std::vector<std::pair<index1_type, index1_type>> crange;
  crange.reserve(contracted_dims.size());
  for (auto d : contracted_dims) {
    auto const& er = l.trange().dim(d).elements_range();
    crange.emplace_back(static_cast<index1_type>(er.first),
                        static_cast<index1_type>(er.second));
  }
  std::vector<index1_type> lidx(la.outer.size()), odo(crange.size());
  for (auto it = c.begin(); it != c.end(); ++it) {
    auto const ctile = it->get();
    for (auto const& cidx : ctile.range()) {
      auto const& inner_c = ctile[cidx];
      if (inner_c.empty()) continue;
      double const vol_c = static_cast<double>(inner_c.range().volume());
      double const mprime = shared_extent_product(inner_c, ca.inner, la.inner);
      if (mprime <= 0) continue;
      for (std::size_t d = 0; d < la.outer.size(); ++d)
        if (l_from_c[d] != npos)
          lidx[d] = static_cast<index1_type>(cidx[l_from_c[d]]);
      for (std::size_t j = 0; j < crange.size(); ++j) odo[j] = crange[j].first;
      while (true) {
        for (std::size_t j = 0; j < crange.size(); ++j)
          lidx[contracted_dims[j]] = odo[j];
        auto const ltile_idx = l.trange().element_to_tile(lidx);
        if (!l.is_zero(ltile_idx)) {
          if (!l.is_local(ltile_idx)) {
            eval::FlopCounter::record_unsourced(
                "left outer tile is remote: " + why);
            return false;
          }
          auto const& inner_l = l.find_local(ltile_idx).get()[lidx];
          if (!inner_l.empty())
            w.gemm(vol_c, mprime,
                   static_cast<double>(inner_l.range().volume()));
        }
        // advance the odometer; an empty odometer means exactly one pass
        bool done = true;
        for (std::size_t j = crange.size(); j-- > 0;) {
          if (++odo[j] < crange[j].second) {
            done = false;
            break;
          }
          odo[j] = crange[j].first;
        }
        if (done) break;
      }
    }
  }
  return true;
}

/// Record a tensor-of-tensor * tensor-of-tensor product that keeps its nest
/// (`ToT * ToT -> ToT`) by walking real inner tensors.
///
/// For one (result outer element, contracted-outer element) pair the work is
/// one GEMM of shape (M, N, K):
///   |inner_C| = M*N*H,  |inner_L| = M*K*H
/// with H the inner Hadamard extent, so with
///   M' := product of the extents (read off inner_C) of the labels shared by
///         inner_L and inner_C   ( = M*H )
///   K  = |inner_L| / M'
/// the flop count is `2 * |inner_C| * K` -- exact, and derived entirely from
/// the two real inner tensors. This is where the per-pair variation of a
/// CSV/PNO array is captured.
///
/// Three walks, all exact, picked purely for cost. When one operand's outer
/// labels cover the result's, walking THAT operand visits each GEMM exactly
/// once and touches only real data -- essential when the contracted outer
/// index is a large one (the DF Kappa, a PAO index). Otherwise the result is
/// walked with an odometer over the contracted outer range.
template <typename LeftT, typename RightT, typename ResultT>
inline void count_product_tot(LeftT const& l, RightT const& r,
                              ResultT const& c, std::string const& lannot,
                              std::string const& rannot,
                              std::string const& cannot) {
  auto const la = split_annot(lannot), ra = split_annot(rannot),
             ca = split_annot(cannot);
  std::string const why = lannot + " * " + rannot + " = " + cannot;
  bool outer_contraction = false;
  for (auto lab : la.outer)
    if (contains(ra.outer, lab) && !contains(ca.outer, lab))
      outer_contraction = true;
  auto covers = [&ca](AnnotLabels const& a) {
    for (auto lab : ca.outer)
      if (!contains(a.outer, lab)) return false;
    return true;
  };

  TotWork w;
  bool ok = true;
  if (outer_contraction && covers(la))
    ok = count_tot_walk_operand(l, c, la, ca, w, why);
  else if (outer_contraction && covers(ra))
    ok = count_tot_walk_operand(r, c, ra, ca, w, why);
  else
    ok = count_tot_walk_result(l, c, la, ra, ca, w, why);
  if (!ok) return;
  // Elements written is a property of the RESULT, counted once: accumulating
  // it per kernel would multiply every result element by the number of
  // contracted-outer terms summed into it.
  eval::FlopCounter::record(eval::FlopCategory::Contraction, w.flops,
                            real_element_count(c));
  eval::FlopCounter::record_kernels(w.kernels, w.m_sum, w.n_sum, w.k_sum);
}

/// Record a `ToT * ToT -> T` (de-nesting) product: every inner mode is
/// contracted, so each (left outer element, result free index) pair
/// contributes one dot product of length `|inner_L|`.
///
/// The result carries no inner index, so the per-GEMM N is 1 and the whole
/// count collapses to `replication * 2 * (total real inner volume of the left
/// operand)`, where `replication` is the product of the extents of the result
/// outer labels the left operand does not carry (free indices supplied by the
/// right operand; the left's inner volume is independent of them). Every left
/// outer label must be free or contracted; anything else is unsourced.
template <typename LeftT, typename RightT, typename ResultT>
inline void count_product_denest(LeftT const& l, RightT const& r,
                                 ResultT const& c, std::string const& lannot,
                                 std::string const& rannot,
                                 std::string const& cannot) {
  auto const la = split_annot(lannot), ra = split_annot(rannot),
             ca = split_annot(cannot);
  for (auto lab : la.outer)
    if (!contains(ca.outer, lab) && !contains(ra.outer, lab)) {
      eval::FlopCounter::record_unsourced(
          "ToT*ToT->T left outer label is neither free nor contracted: " +
          lannot + " * " + rannot + " = " + cannot);
      return;
    }
  double replication = 1;
  for (auto lab : ca.outer)
    if (!contains(la.outer, lab)) {
      double const e = mode_extent(r, ra.outer, lab);
      if (e == 0) {
        eval::FlopCounter::record_unsourced(
            "ToT*ToT->T result outer label has no extent on either operand: " +
            lannot + " * " + rannot + " = " + cannot);
        return;
      }
      replication *= e;
    }
  double flops = 0;
  for (auto it = l.begin(); it != l.end(); ++it) {
    auto const tile = it->get();
    for (auto const& inner : tile) {
      if (inner.empty()) continue;
      flops += 2 * static_cast<double>(inner.range().volume());
    }
  }
  // No shape census: these are rank-1 dot products (M = N = 1), and there are
  // orders of magnitude more of them than there are real GEMMs, so counting
  // them would bury the GEMM shape statistics the ToT*ToT->ToT walks measure.
  eval::FlopCounter::record(eval::FlopCategory::Contraction,
                            replication * flops, real_element_count(c));
}

}  // namespace sequant::detail

#endif  // SEQUANT_HAS_TILEDARRAY
#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_FLOPS_HPP
