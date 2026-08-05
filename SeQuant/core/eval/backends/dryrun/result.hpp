#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <format>

#include <algorithm>
#include <any>
#include <array>
#include <cstddef>
#include <memory>
#include <numeric>
#include <utility>

namespace sequant::eval::dryrun {

///
/// \brief Annotation type DryRun's Result ops decode from the eval engine's
///        std::any [l,r,res] / [pre,post] triples/pairs.
///
/// Unlike \c EvalExprTAPP (opaque \c int64_t index-label hashes -- see
/// \c backends/tapp/eval_expr.hpp), DryRun's annotation IS the plain literal
/// (canon-order) index list itself: \c Result::prod/sum/permute need each
/// index's actual space/extent (via \c CostModel), not just an opaque
/// identity, to compute a modeled size.
///
using annot_t = container::svector<Index>;

/// Per-mode assembled element coverage recorded by write_into_slice(): maps an
/// outer mode position to the contiguous `[lo, hi)` element range filled so far
/// by scattered blocks. Lets a zero-data DryRun destination report the REALIZED
/// (assembled) size along a partitioned mode and detect gaps/overlaps between
/// blocks -- the assemble-side analogue of ExtentOverrides for slice_mode().
using AssembledCoverage =
    container::map<std::size_t, std::pair<std::size_t, std::size_t>>;

class ResultDryRun;
class ResultDryRunNested;

///
/// \brief Builds whichever concrete DryRun Result type matches \p idx's
///        content: a nested \c ResultDryRunNested if any index in \p idx is
///        proto-indexed (a CSV/PNO composite leg, e.g. a CSV amplitude's PNO
///        domain leg `a_1<i_1,i_2>`), otherwise a flat \c ResultDryRun.
///
/// Dispatch is by CONTENT of the decoded result annotation, not by either
/// operand's concrete type -- exactly mirroring how the real eval engine
/// itself decides tensor-of-tensor-ness (\c EvalExpr::tot(), from the same
/// proto-indexed-leg criterion). This is what lets \c prod()/sum() freely
/// combine a flat operand (e.g. a bare 3-center DF integral) with a nested
/// one (e.g. a CSV/PNO coefficient), exactly as real CSV-CCSD terms do,
/// without either side needing to know the other's concrete type.
///
[[nodiscard]] inline ResultPtr make_dryrun_result(
    container::svector<Index> idx, std::shared_ptr<CostModel const> cm,
    ExtentOverrides overrides = {});

namespace detail {

[[nodiscard]] inline bool has_proto(container::svector<Index> const& idx) {
  return std::any_of(idx.begin(), idx.end(),
                     [](Index const& ix) { return ix.has_proto_indices(); });
}

[[nodiscard]] inline ExtentOverrides merge_overrides(ExtentOverrides const& a,
                                                     ExtentOverrides const& b) {
  ExtentOverrides out = a;
  for (auto const& [ix, n] : b) out[ix] = n;
  return out;
}

// Uniform read access to a DryRun Result's (index list, overrides, cost
// model) regardless of which concrete DryRun type `r` is. `is<T>()`/`as<T>()`
// are public Result methods, so no friendship is needed; declared here (and
// defined below, after both concrete classes) purely because their bodies
// need the concrete classes' definitions.
[[nodiscard]] container::svector<Index> indices_of(Result const& r);
[[nodiscard]] ExtentOverrides overrides_of(Result const& r);

///
/// \brief Shared op bodies for the two DryRun Result concrete types.
///
/// Both \c ResultDryRun and \c ResultDryRunNested carry exactly an (index
/// list, ExtentOverrides, CostModel) triple and differ only in what they
/// additionally expose (\c ResultDryRunNested splits its index list into
/// outer()/inner() views for CSV-composite-aware inspection/testing).
/// Implemented once here so the two classes' prod/sum/permute/slice_mode/
/// mode_batches bodies are one-line forwards, not near-duplicated logic.
///
struct DryRunOps {
  [[nodiscard]] static ResultPtr sum(container::svector<Index> const& idx,
                                     ExtentOverrides const& ov,
                                     std::shared_ptr<CostModel const> const& cm,
                                     Result const& other,
                                     std::array<std::any, 3> const& annot) {
    auto const a = Annot<annot_t>{annot};
    auto merged = merge_overrides(ov, overrides_of(other));
    return make_dryrun_result(
        container::svector<Index>(a.this_annot.begin(), a.this_annot.end()), cm,
        std::move(merged));
  }

  [[nodiscard]] static ResultPtr prod(
      container::svector<Index> const& idx, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm, Result const& other,
      std::array<std::any, 3> const& annot) {
    if (other.is<ResultScalar<double>>()) {
      // Scalar * tensor: shape (and any accumulated slicing) unchanged.
      return make_dryrun_result(idx, cm, ov);
    }
    auto const a = Annot<annot_t>{annot};
    auto merged = merge_overrides(ov, overrides_of(other));

    // Emit the cost model's OWN flops / roofline exec_cost for THIS op into the
    // eval trace (gated on the eval log level), interleaved right before the
    // generic engine's `Eval | Product` line for the same op. This lets trace
    // post-processing weight avoidable recomputation by MODELLED TIME without
    // re-deriving the cost downstream (which would silently drift from the
    // model). The (out, contracted) index sets and the realized (sliced) extent
    // overrides feed the same CostModel closures the static cost_profile() walk
    // uses, so per-op costs are consistent with the whole-forest totals.
    if (Logger::instance().eval.level > 0) {
      container::svector<Index> const rhs = indices_of(other);
      container::svector<Index> out(a.this_annot.begin(), a.this_annot.end());
      container::svector<Index> contracted;
      for (auto const& ix : idx)
        if (std::find(rhs.begin(), rhs.end(), ix) != rhs.end())
          contracted.push_back(ix);
      double const flops = cm->flops(out, contracted, merged);
      double const exec = cm->exec_cost(flops, cm->memsize(idx, ov), 4096);
      write_log(Logger::instance(), "OpCost", std::format(" | {}", flops),
                std::format(" | {}", exec), '\n');
      // Fold this op's SLICED-extent cost into the replay cost sink, if one is
      // attached (cost_profile()'s recompute-aware tally). merged/ov carry the
      // runtime slicing, so a contraction re-executed once per occ block is
      // charged once per block at its sliced size -- the same numbers already
      // logged above, now summed. No-op (byte-identical) when unattached.
      cm->tally_op(flops, exec);
      // Per-node AVOIDABLE-recompute tally, keyed by the LABEL signature (the
      // same one the runtime Build event and IR node carry, so the visualizer
      // joins by hash->sig). Avoidable is measured in FLOPs against the
      // batching-free ideal of building each value ONCE at full extent:
      //  - `flops`      is THIS build's actual (sliced) FLOPs, accumulated;
      //  - `full_flops` is the FLOPs to build the value once at FULL extent (no
      //    overrides). FLOPs is linear in extents, so disjoint per-block slices
      //    that tile the full value sum to exactly full_flops (0 avoidable),
      //    while a value rebuilt full per block sums to N*full ((N-1)*full
      //    avoidable). No slice-context is needed in the key.
      double const full_flops = cm->flops(out, contracted, {});
      std::string const sig = cost_op_signature(out, idx, rhs);
      cm->tally_node(sig, flops, full_flops);
    }

    if (a.this_annot.empty()) {
      // Full contraction -> scalar. No real numeric value is ever tracked by
      // this zero-data backend (only sizes/costs), so the placeholder 0.0
      // is never meant to be read as a physical result.
      return eval_result<ResultScalar<double>>(0.0);
    }
    return make_dryrun_result(
        container::svector<Index>(a.this_annot.begin(), a.this_annot.end()), cm,
        std::move(merged));
  }

  [[nodiscard]] static ResultPtr permute(
      container::svector<Index> const& /*idx*/, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm,
      std::array<std::any, 2> const& ann) {
    auto const post = std::any_cast<annot_t>(ann[1]);
    return make_dryrun_result(
        container::svector<Index>(post.begin(), post.end()), cm, ov);
  }

  [[nodiscard]] static ResultPtr slice_mode(
      container::svector<Index> const& idx, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm, std::size_t mode,
      std::size_t elem_lo, std::size_t elem_hi) {
    SEQUANT_ASSERT(mode < idx.size());
    auto merged = ov;
    merged[idx[mode]] = elem_hi - elem_lo;
    return make_dryrun_result(idx, cm, std::move(merged));
  }

  /// Scatter \p block into the `[block_lo, block_hi)` element slice of the
  /// destination's mode \p mode -- the inverse of slice_mode(). Zero-data:
  /// updates only the destination's modelled size and assembled-coverage
  /// bookkeeping. \p ov and \p cov are the destination's (mutated in place).
  static void write_into_slice(container::svector<Index> const& idx,
                               ExtentOverrides& ov, AssembledCoverage& cov,
                               std::shared_ptr<CostModel const> const& cm,
                               Result const& block, std::size_t mode,
                               std::size_t block_lo, std::size_t block_hi) {
    SEQUANT_ASSERT(mode < idx.size());
    SEQUANT_ASSERT(block_lo < block_hi);
    Index const& mix = idx[mode];
    // Tile/width consistency: the block's own modelled extent on the shared
    // mode index must equal the slice width it is being written into.
    auto const bov = overrides_of(block);
    std::size_t const block_extent = [&] {
      if (auto it = bov.find(mix); it != bov.end()) return it->second;
      return cm->regime().extent(mix);
    }();
    SEQUANT_ASSERT(block_extent == block_hi - block_lo);
    // Merge the block's range into the assembled coverage, requiring
    // contiguity: a block that neither appends after nor prepends before the
    // filled range would leave a gap or overlap another block (a
    // double-count). This is what makes disjoint gap-free tiling the only
    // accepted assembly.
    if (auto it = cov.find(mode); it == cov.end()) {
      cov.emplace(mode,
                  std::pair<std::size_t, std::size_t>{block_lo, block_hi});
    } else {
      auto& lohi = it->second;
      bool const append = block_lo == lohi.second;
      bool const prepend = block_hi == lohi.first;
      SEQUANT_ASSERT(append || prepend);
      if (append)
        lohi.second = block_hi;
      else
        lohi.first = block_lo;
    }
    // Reflect the assembled element width (hi - lo, lobound preserved) as the
    // realized extent of the batch mode so size_in_bytes() tracks the
    // reconstructed footprint.
    auto const& lohi = cov.at(mode);
    ov[mix] = lohi.second - lohi.first;
  }

  /// Build a zero-data destination shaped like \p idx but with mode \p mode's
  /// index widened to \p axis_src's FULL (unsliced) extent at \p
  /// axis_src_mode -- the dry-run analogue of \c ResultTensorTA:: and \c
  /// ResultTensorOfTensorTA::pre_sized_zeros_over_mode's outer-\c
  /// TiledRange1 swap. \p axis_src is the unsliced carrier leaf's token (its
  /// OWN recorded override for the axis index, if any, else the CostModel
  /// regime's natural extent, is the "full" extent -- mirroring how the TA
  /// backend reads the swapped-in dimension straight off \p axis_src rather
  /// than assuming a fixed default). The width recorded here is a structural
  /// (shape) fact, queryable via \c size_in_bytes()/overrides() immediately
  /// -- exactly as the TA test \c batched_scratch_tot_presize_scatter checks
  /// \c trange().dim(mode) right after presizing, before any block is
  /// written. Once the scatter loop's \c write_into_slice() calls begin, the
  /// AssembledCoverage bookkeeping there takes over the mode's reported
  /// extent (the REALIZED/assembled width so far); by the time every
  /// disjoint block has been written the assembled width converges back to
  /// this same full extent.
  [[nodiscard]] static ResultPtr pre_sized_zeros_over_mode(
      container::svector<Index> const& idx, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm, std::size_t mode,
      Result const& axis_src, std::size_t axis_src_mode) {
    SEQUANT_ASSERT(mode < idx.size());
    auto const src_idx = indices_of(axis_src);
    SEQUANT_ASSERT(axis_src_mode < src_idx.size());
    Index const& axis_ix = src_idx[axis_src_mode];
    auto const src_ov = overrides_of(axis_src);
    std::size_t const full_extent = [&] {
      if (auto it = src_ov.find(axis_ix); it != src_ov.end()) return it->second;
      return cm->regime().extent(axis_ix);
    }();
    auto merged = ov;
    merged[idx[mode]] = full_extent;
    return make_dryrun_result(idx, cm, std::move(merged));
  }

  [[nodiscard]] static container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(container::svector<Index> const& idx, ExtentOverrides const& ov,
               std::shared_ptr<CostModel const> const& cm, std::size_t mode,
               std::size_t target_batch_size) {
    SEQUANT_ASSERT(mode < idx.size());
    Index const& ix = idx[mode];
    std::size_t extent;
    if (auto it = ov.find(ix); it != ov.end())
      extent = it->second;
    else
      extent = cm->regime().extent(ix);

    container::svector<std::pair<std::size_t, std::size_t>> out;
    if (target_batch_size == 0 || extent == 0) {
      out.push_back({0, extent});
      return out;
    }
    // CALLER-SUPPLIED PARTITION: if the mode's space has a recorded batch
    // partition (SizeRegime::space_slice_extents) AND the mode is at full
    // extent (the slice extents sum to `extent`, i.e. this is the batch axis,
    // not a spectator already narrowed by `ov`), use it directly -- accumulate
    // the slice extents into contiguous [lo,hi) ranges. `target_batch_size` is
    // ignored here: the caller already applied it when it built the partition
    // (e.g. via batch_slice_extents_from_tiles for a tiled backend). The
    // dry-run backend stays model-agnostic -- it just reads the slice extents.
    auto const& slices = cm->regime().slice_extents(ix);
    std::size_t const slice_sum =
        std::accumulate(slices.begin(), slices.end(), std::size_t{0});
    if (!slices.empty() && slice_sum == extent) {
      std::size_t lo = 0;
      for (std::size_t const s : slices) {
        out.push_back({lo, lo + s});
        lo += s;
      }
      return out;
    }
    // Fallback: uniform target_batch_size blocks (no partition recorded).
    for (std::size_t lo = 0; lo < extent; lo += target_batch_size)
      out.push_back({lo, std::min(extent, lo + target_batch_size)});
    return out;
  }
};

}  // namespace detail

///
/// \brief Flat (non-CSV) zero-data tensor token.
///
/// Carries only its own literal outer index list (canon order -- the same
/// order \c EvalExpr::canon_indices()/annot() use, so \c slice_mode()/
/// \c mode_batches()'s positional `mode` argument indexes it correctly), an
/// \c ExtentOverrides table recording any runtime \c slice_mode()/
/// \c mode_batches() narrowing (keyed by Index so it survives reshaping
/// across prod/sum/permute), and a shared \c CostModel. No tensor data is
/// ever allocated or copied; every op is index-set bookkeeping plus a
/// CostModel query. Mirrors \c ResultTensorTAPP's structure
/// (backends/tapp/result.hpp) with every real-tensor line replaced by that
/// bookkeeping.
///
class ResultDryRun final : public Result {
 public:
  using Result::id_t;

  ResultDryRun(container::svector<Index> idxset,
               std::shared_ptr<CostModel const> cm,
               ExtentOverrides overrides = {})
      : Result{Payload{}},
        indices_{std::move(idxset)},
        cm_{std::move(cm)},
        overrides_{std::move(overrides)} {}

  [[nodiscard]] container::svector<Index> const& indices() const noexcept {
    return indices_;
  }
  [[nodiscard]] ExtentOverrides const& overrides() const noexcept {
    return overrides_;
  }

  /// The contiguous `[lo, hi)` element range of outer mode \p mode assembled so
  /// far by write_into_slice() (empty `{0, 0}` if nothing written).
  [[nodiscard]] std::pair<std::size_t, std::size_t> assembled_range(
      std::size_t mode) const {
    if (auto it = assembled_.find(mode); it != assembled_.end())
      return it->second;
    return {0, 0};
  }

 private:
  struct Payload {};

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<ResultDryRun>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    return detail::DryRunOps::sum(indices_, overrides_, cm_, other, annot);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    return detail::DryRunOps::prod(indices_, overrides_, cm_, other, annot);
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann);
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann);
  }

  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    return detail::DryRunOps::slice_mode(indices_, overrides_, cm_, mode,
                                         elem_lo, elem_hi);
  }

  [[nodiscard]] container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(std::size_t mode, std::size_t target_batch_size) const override {
    return detail::DryRunOps::mode_batches(indices_, overrides_, cm_, mode,
                                           target_batch_size);
  }

  void write_into_slice(Result const& block, std::size_t mode,
                        std::size_t block_lo, std::size_t block_hi) override {
    detail::DryRunOps::write_into_slice(indices_, overrides_, assembled_, cm_,
                                        block, mode, block_lo, block_hi);
  }

  [[nodiscard]] ResultPtr pre_sized_zeros_over_mode(
      std::size_t mode, Result const& axis_src,
      std::size_t axis_src_mode) const override {
    return detail::DryRunOps::pre_sized_zeros_over_mode(
        indices_, overrides_, cm_, mode, axis_src, axis_src_mode);
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<ResultDryRun>() || other.is<ResultDryRunNested>());
    overrides_ =
        detail::merge_overrides(overrides_, detail::overrides_of(other));
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_);
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_);
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t /*factor*/) const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_);
  }

  [[nodiscard]] std::size_t size_in_bytes() const final {
    return cm_->memsize(indices_, overrides_);
  }

  container::svector<Index> indices_;
  std::shared_ptr<CostModel const> cm_;
  ExtentOverrides overrides_;
  AssembledCoverage assembled_;
};

///
/// \brief CSV/PNO tensor-of-tensor zero-data token.
///
/// Like \c ResultDryRun, but additionally exposes an outer()/inner() split of
/// its (canon-order) index list -- inner = the proto-indexed (composite)
/// legs, e.g. a CSV amplitude's PNO domain leg `a_1<i_1,i_2>`; outer = every
/// other (plain) leg, e.g. the PAO index `mu~_1`. The split is purely an
/// observability/testing convenience: \c size_in_bytes()'s arithmetic is
/// IDENTICAL to \c ResultDryRun's (\c CostModel::memsize already routes any
/// index list containing a proto-indexed entry through the moment-aware
/// `inner_pow` path internally, via \c tot_indices/inner_aware_volume --
/// content-driven, not type-driven), so tests that want to confirm "this used
/// the k-th moment, not extent^k" can inspect inner() directly.
///
/// Position semantics for \c slice_mode()/\c mode_batches(): the `mode`
/// argument the runtime passes is always resolved against the FULL
/// canon-order list (an optional trailing constructor argument, defaulting to
/// `outer ++ inner` when the caller does not need position accuracy, e.g. a
/// hand-built test instance); the \c DryRunLeafEvaluator (eval_expr.hpp)
/// always supplies the leaf's true \c canon_indices() order there, since only
/// LEAF-constructed instances are ever sliced by the runtime (\c slice_mode()
/// is invoked only inside the batched evaluator's leaf-wrapping closure, never
/// on a prod()/sum()-produced intermediate).
///
class ResultDryRunNested final : public Result {
 public:
  using Result::id_t;

  ResultDryRunNested(container::svector<Index> outer,
                     container::svector<Index> inner,
                     std::shared_ptr<CostModel const> cm,
                     ExtentOverrides overrides = {},
                     container::svector<Index> canon_order = {})
      : Result{Payload{}},
        outer_{std::move(outer)},
        inner_{std::move(inner)},
        indices_{canon_order.empty()
                     ? [this] {
                         container::svector<Index> c = outer_;
                         c.insert(c.end(), inner_.begin(), inner_.end());
                         return c;
                       }()
                     : std::move(canon_order)},
        cm_{std::move(cm)},
        overrides_{std::move(overrides)} {}

  [[nodiscard]] container::svector<Index> const& outer() const noexcept {
    return outer_;
  }
  [[nodiscard]] container::svector<Index> const& inner() const noexcept {
    return inner_;
  }
  [[nodiscard]] container::svector<Index> const& indices() const noexcept {
    return indices_;
  }
  [[nodiscard]] ExtentOverrides const& overrides() const noexcept {
    return overrides_;
  }

  /// The contiguous `[lo, hi)` element range of outer mode \p mode assembled so
  /// far by write_into_slice() (empty `{0, 0}` if nothing written).
  [[nodiscard]] std::pair<std::size_t, std::size_t> assembled_range(
      std::size_t mode) const {
    if (auto it = assembled_.find(mode); it != assembled_.end())
      return it->second;
    return {0, 0};
  }

 private:
  struct Payload {};

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<ResultDryRunNested>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    return detail::DryRunOps::sum(indices_, overrides_, cm_, other, annot);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    return detail::DryRunOps::prod(indices_, overrides_, cm_, other, annot);
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann);
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann);
  }

  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    return detail::DryRunOps::slice_mode(indices_, overrides_, cm_, mode,
                                         elem_lo, elem_hi);
  }

  [[nodiscard]] container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(std::size_t mode, std::size_t target_batch_size) const override {
    return detail::DryRunOps::mode_batches(indices_, overrides_, cm_, mode,
                                           target_batch_size);
  }

  void write_into_slice(Result const& block, std::size_t mode,
                        std::size_t block_lo, std::size_t block_hi) override {
    detail::DryRunOps::write_into_slice(indices_, overrides_, assembled_, cm_,
                                        block, mode, block_lo, block_hi);
  }

  [[nodiscard]] ResultPtr pre_sized_zeros_over_mode(
      std::size_t mode, Result const& axis_src,
      std::size_t axis_src_mode) const override {
    return detail::DryRunOps::pre_sized_zeros_over_mode(
        indices_, overrides_, cm_, mode, axis_src, axis_src_mode);
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<ResultDryRun>() || other.is<ResultDryRunNested>());
    overrides_ =
        detail::merge_overrides(overrides_, detail::overrides_of(other));
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_);
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_);
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t /*factor*/) const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_);
  }

  [[nodiscard]] std::size_t size_in_bytes() const final {
    return cm_->memsize(indices_, overrides_);
  }

  container::svector<Index> outer_;
  container::svector<Index> inner_;
  container::svector<Index> indices_;  // canon order; outer_++inner_ content
  std::shared_ptr<CostModel const> cm_;
  ExtentOverrides overrides_;
  AssembledCoverage assembled_;
};

[[nodiscard]] inline ResultPtr make_dryrun_result(
    container::svector<Index> idx, std::shared_ptr<CostModel const> cm,
    ExtentOverrides overrides) {
  if (!detail::has_proto(idx))
    return eval_result<ResultDryRun>(std::move(idx), std::move(cm),
                                     std::move(overrides));
  container::svector<Index> outer, inner;
  for (auto const& ix : idx)
    (ix.has_proto_indices() ? inner : outer).push_back(ix);
  return eval_result<ResultDryRunNested>(std::move(outer), std::move(inner),
                                         std::move(cm), std::move(overrides),
                                         std::move(idx));
}

namespace detail {

[[nodiscard]] inline container::svector<Index> indices_of(Result const& r) {
  if (r.is<ResultDryRun>()) return r.as<ResultDryRun>().indices();
  SEQUANT_ASSERT(r.is<ResultDryRunNested>());
  return r.as<ResultDryRunNested>().indices();
}

[[nodiscard]] inline ExtentOverrides overrides_of(Result const& r) {
  if (r.is<ResultDryRun>()) return r.as<ResultDryRun>().overrides();
  SEQUANT_ASSERT(r.is<ResultDryRunNested>());
  return r.as<ResultDryRunNested>().overrides();
}

}  // namespace detail

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP
