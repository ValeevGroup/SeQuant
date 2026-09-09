#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_MODEL_OBJECT_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_MODEL_OBJECT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize/cost_model.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>

#include <cstddef>
#include <functional>
#include <utility>

namespace sequant::eval::dryrun {

/// Per-index extent OVERRIDE table: narrows specific indices (by identity, so
/// it survives reshaping across prod/sum/permute -- the same shared/
/// contracted Index object may occupy different tensor modes at different
/// nodes) to a runtime-realized element count. Populated by
/// Result::slice_mode()/mode_batches() call sites (see result.hpp); empty =>
/// no override, the regime's nominal extent applies. This table -- not a
/// second cost model -- is what lets a zero-data DryRun Result report the
/// REALIZED (possibly runtime-sliced) size rather than always the full
/// regime extent, which is exactly the signal Task 6's replay witnesses.
using ExtentOverrides = container::map<Index, std::size_t>;

///
/// \brief Bundles the optimizer's own cost closures (memsize/flops/roofline)
/// behind one value type so dry-run Results report MODEL size (not an
/// allocated size), and the harness can additionally read FLOPs and
/// projected execution cost per operation.
///
/// This is a thin wrapper: all arithmetic is delegated verbatim to
/// \c sequant::opt::detail::memsize_counter / \c flops_counter / \c
/// roofline_op_cost (see \c core/optimize/single_term_detail.hpp and \c
/// core/optimize/cost_model.hpp) -- no parallel cost model is implemented
/// here. The only thing this class adds is the ExtentOverrides indirection:
/// each query builds a fresh (cheap; no heap allocation beyond the closure
/// itself) index-to-extent callable that consults \p overrides before
/// falling back to the SizeRegime's nominal extent, then hands that callable
/// to the counter.
///
class CostModel {
 public:
  explicit CostModel(SizeRegime regime, RooflineParams roofline = {})
      : regime_{std::move(regime)}, roofline_{roofline} {}

  ///
  /// \brief Bytes for a tensor with these (literal, canon-order) indices,
  ///        honoring any per-index extent override (a runtime slice_mode()/
  ///        mode_batches() narrowing).
  ///
  /// Delegates the extent-product / composite-moment math to \c
  /// memsize_counter, invoked with \p idxset as the sole (`lhs`) operand and
  /// empty `rhs`/`result` -- an empty operand's tot_indices() split
  /// accumulates the starting product of 1.0, which memsize_counter itself
  /// special-cases to contribute zero bytes, so this reproduces exactly the
  /// single-operand byte count \c memsize_counter is designed to report per
  /// operand.
  ///
  [[nodiscard]] std::size_t memsize(
      container::svector<Index> const& idxset,
      ExtentOverrides const& overrides = {}) const {
    auto const ext = make_extent_fn(overrides);
    auto const mc =
        sequant::opt::detail::memsize_counter(ext, regime_.inner_pow_fn());
    double const elems =
        mc(idxset, container::svector<Index>{}, container::svector<Index>{});
    return static_cast<std::size_t>(elems * numeric_size_);
  }

  ///
  /// \brief Multiply-add count for a contraction whose free (result) indices
  ///        are \p out and whose contracted (summed-over) indices are
  ///        \p contracted.
  ///
  /// Delegates to \c flops_counter, which prices the union of its (lhs, rhs,
  /// result) arguments; passing (\p out, \p contracted, {}) makes that union
  /// exactly `out U contracted` -- the full index set touched by the
  /// contraction, since by construction `contracted` holds precisely the
  /// indices present in both operands but absent from the result.
  ///
  [[nodiscard]] double flops(container::svector<Index> const& out,
                             container::svector<Index> const& contracted,
                             ExtentOverrides const& overrides = {}) const {
    auto const ext = make_extent_fn(overrides);
    auto const fc =
        sequant::opt::detail::flops_counter(ext, regime_.inner_pow_fn());
    return fc(out, contracted, container::svector<Index>{});
  }

  ///
  /// \brief Roofline-projected execution cost of one contraction (see
  ///        \c sequant::opt::detail::roofline_op_cost).
  ///
  /// \p left_bytes / \p right_bytes are operand footprints in BYTES (as
  /// reported by \c Result::size_in_bytes()); converted to elements (the
  /// counter's native unit) via \c numeric_size before delegating.
  ///
  [[nodiscard]] double exec_cost(double flops_count, std::size_t left_bytes,
                                 std::size_t right_bytes) const {
    double const traffic_elems =
        static_cast<double>(left_bytes + right_bytes) / numeric_size_;
    return sequant::opt::detail::roofline_op_cost(
        flops_count, traffic_elems, roofline_.machine_balance,
        roofline_.fast_mem_elems, roofline_.block_tiles,
        roofline_.block_prefactor);
  }

  [[nodiscard]] SizeRegime const& regime() const noexcept { return regime_; }

 private:
  // Index-to-extent callable consulting `overrides` first, else the
  // regime's nominal extent. The returned std::function captures `overrides`
  // (and `this`) BY REFERENCE and is only ever used -- never stored --
  // within the (memsize/flops) call that constructs it, so the reference
  // stays valid for its entire lifetime. Explicit (non-deduced) return type
  // so this can be called from memsize()/flops(), which appear earlier in
  // the class body (a deduced `auto` return type would require the
  // definition to precede every use, even within the same class).
  [[nodiscard]] std::function<std::size_t(Index const&)> make_extent_fn(
      ExtentOverrides const& overrides) const {
    return [this, &overrides](Index const& ix) -> std::size_t {
      if (auto it = overrides.find(ix); it != overrides.end())
        return it->second;
      return regime_.extent(ix);
    };
  }

  SizeRegime regime_;
  RooflineParams roofline_;
  // sizeof(double); see doc/dev/plans/2026-07-04-dryrun-eval-backend.md Task 2
  // note on OptimizeOptions::numeric_size (hardcoded here, matching the C60
  // trace's real-only CSV-CCk path; complex CSV-CCk is out of scope, see the
  // plan's carried-minor N4).
  double numeric_size_ = 8.0;
};

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_MODEL_OBJECT_HPP
