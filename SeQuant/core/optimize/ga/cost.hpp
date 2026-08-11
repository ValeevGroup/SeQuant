#ifndef SEQUANT_CORE_OPTIMIZE_GA_COST_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_COST_HPP

// Mask-based fast path for the per-merge flop count, reproducing
// opt::detail::flops_counter (with empty inner_pow) exactly: the extent
// product over the deduplicated union of the indices of lhs, rhs and result,
// including the proto-indices of any composite.
//
// No unit test calls flops_counter on the same merge, so the agreement is a
// contract, not a checked one -- if you change either, re-derive the other by
// hand. What the [ga] tests DO pin, indirectly: the python-parity fixtures in
// tests/unit/ga_reference.hpp (exact total costs of whole schedules under the
// prototype convention, which is this formula x2) and `RuntimeCoster` in
// tests/unit/test_ga.cpp (an independent extent-product-over-index-closure
// walk of the EMITTED forest, reconciled against the schedule's own totals).
// Both would break on a change to the index closure here; neither would catch
// a divergence from flops_counter that they happen to share.

#include <SeQuant/core/optimize/ga/key_table.hpp>

namespace sequant::opt::ga {

/// Extent product over the union of the open indices of subsets \p a, \p b and
/// \p a|b of \p T. == flops_counter(F(a), F(b), F(a|b)) except for the
/// scalar-contraction rule, which the caller's CostModel applies.
inline double merge_volume(TermTable const& T, NodeMask a, NodeMask b) {
  const NodeMask ab = a | b;
  IndexMask m = T.outer_mask[a] | T.outer_mask[b] | T.outer_mask[ab] |
                T.inner_mask[a] | T.inner_mask[b] | T.inner_mask[ab];
  double v = 1.;
  for (; m; m &= m - 1) v *= T.extent[std::countr_zero(m)];
  return v;
}

/// Cost conventions of the two accounting modes. Native (default) is
/// SeQuant's flops_counter convention: a merge costs its extent-product,
/// scalar contractions are free, an L2 addition costs the face size.
/// Python-parity is the prototype's -- 2x the extent-product, no
/// scalar-contraction rule -- and exists only to validate against the
/// prototype's reference numbers. Both are pure FLOP models; nothing here
/// prices storage, and \ref volatile_weight is a REPLAY COUNT on the flops
/// already counted, not a second metric.
struct CostModel {
  double merge_factor = 1.0;
  bool zero_scalar_merges = true;
  /// Replay weight on volatile (amplitude-dependent) work: only that part of
  /// the DAG is rebuilt per solver iteration, so the honest objective is
  /// `persistent_flops + N * volatile_flops`. Mirrors
  /// CostParams::volatile_weight; 1.0 is the volatility-blind cost.
  double volatile_weight = 1.0;
  /// Widen the amortization class from "shared" to "shared OR persistent":
  /// emission NAMES every amplitude-free keyed array it can (not only the
  /// multiply-used ones), so the runtime evaluates it once per solve and the
  /// objective may charge it once. Off by default = the historical contract,
  /// bit for bit, on every path. Read ONLY through \ref runtime_amortized, so
  /// `Fitness::resolve` and `emit.cpp` cannot disagree about what it means.
  bool amortize_persistent = false;
  /// Hard per-array cap, in elements, on the class \ref amortize_persistent
  /// adds. Every newly named persistent array becomes resident for the whole
  /// solve (mpqc4's `GAImedStore` never evicts), so an uncapped rule can hoist
  /// an arbitrarily large intermediate. A key over the cap is simply not
  /// amortized: it stays inlined at emission AND stays charged
  /// x volatile_weight in the objective -- the two sides move together because
  /// the cap is read inside the shared predicate. Shared (use_count >= 2) keys
  /// are named regardless, exactly as before, so `naming_cap_elems == 0`
  /// degrades to the pre-\ref amortize_persistent contract. Default 2e8
  /// elements (~1.6 GB fp64) is deliberately conservative and meant to be
  /// tuned against measured peak RSS, not trusted.
  double naming_cap_elems = 2e8;

  static CostModel native() { return {1.0, true}; }
  static CostModel python_parity() { return {2.0, false}; }

  /// Whether cluster \p S has to be rebuilt on every replay: true iff it
  /// contains an amplitude leaf. Monotone in S, so no propagation pass is
  /// needed for L1 clusters.
  static bool is_volatile(TermTable const& T, NodeMask S) {
    return (S & T.volatile_mask) != 0;
  }
  /// How many times work of this volatility is paid over the whole solve.
  double replays(bool volatile_work) const {
    return volatile_work ? volatile_weight : 1.0;
  }

  /// **The matched pair.** Whether the array produced by cluster \p S of term
  /// \p T is built ONCE per solve rather than rebuilt on every replay, given
  /// whether it is rendered at more than one site (\p shared).
  ///
  /// This one function is both the objective's replay rule (`Fitness::resolve`
  /// charges `x volatile_weight` exactly where this is false) and emission's
  /// naming rule (`emit.cpp` names a key iff it is `shared` or this says the
  /// runtime would amortize it). They cannot drift because there is nothing to
  /// keep in sync: promising amortization the emission does not deliver makes
  /// the objective lie, and naming what the objective still charges per replay
  /// wastes residency for nothing.
  ///
  /// - volatile work is never amortized: its value changes every iteration;
  /// - `shared` work is amortized today and always was (emission names it);
  /// - persistent single-use work is amortized only under
  ///   \ref amortize_persistent, and only if it can actually be named: a
  ///   scalar-faced cluster has no slots (`define()`'s `!slots.empty()` assert,
  ///   and mpqc4's, are BOTH compiled out in release, so this guard is the only
  ///   thing standing between such a cluster and a silently slotless head), and
  ///   an over-cap cluster is refused residency (\ref naming_cap_elems).
  ///
  /// \note Volatility-blind tables must short-circuit this entirely: callers
  ///       gate on `KeyTable::volatility_aware` (the L2 replay flag), which is
  ///       not visible from here.
  bool runtime_amortized(TermTable const& T, NodeMask S, bool shared) const {
    if (is_volatile(T, S)) return false;
    if (shared) return true;
    return amortize_persistent && !T.canon_face_bits[S].empty() &&
           T.face_size(S) <= naming_cap_elems;
  }

  /// \p volatile_work: whether the merged value depends on the amplitudes.
  double merge(TermTable const& T, NodeMask a, NodeMask b,
               bool volatile_work) const {
    const double v = merge_volume(T, a, b);
    if (zero_scalar_merges && v == 1.) return 0.;
    return replays(volatile_work) * merge_factor * v;
  }
  /// L1 form: the merge produces cluster `a|b`, whose volatility is decided by
  /// its own tensors.
  double merge(TermTable const& T, NodeMask a, NodeMask b) const {
    return merge(T, a, b, is_volatile(T, a | b));
  }
  double addition(TermTable const& T, NodeMask S, bool volatile_work) const {
    return replays(volatile_work) * T.face_size(S);
  }
  double addition(TermTable const& T, NodeMask S) const {
    return addition(T, S, is_volatile(T, S));
  }
};

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_COST_HPP
