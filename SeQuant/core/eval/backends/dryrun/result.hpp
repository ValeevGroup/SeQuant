#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backend_array_ops.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <format>
#include <stdexcept>

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
    ExtentOverrides overrides = {}, ExtentOverrides lobounds = {});

namespace detail {

[[nodiscard]] inline bool has_proto(container::svector<Index> const& idx) {
  return std::any_of(idx.begin(), idx.end(),
                     [](Index const& ix) { return ix.has_proto_indices(); });
}

[[nodiscard]] inline ExtentOverrides merge_overrides(ExtentOverrides const& a,
                                                     ExtentOverrides const& b) {
  ExtentOverrides out = a;
  for (auto const& [pos, n] : b) out[pos] = n;
  return out;
}

// Remap positional overrides from one annotation's mode positions to another's
// by matching labels: position p (labeled `from[p]`) -> the position of that
// same label in `to`. A position whose label is absent from `to` (e.g. an index
// this op contracts away) is dropped -- it has no mode in the result. This is
// how a sliced mode's width survives prod/sum/permute now that the value itself
// carries no labels: the op's annotation is the only place labels live, so two
// positional maps from different operands can only be combined after both are
// projected onto a COMMON annotation (the result's).
[[nodiscard]] inline ExtentOverrides remap_overrides_by_annot(
    ExtentOverrides const& ov, annot_t const& from, annot_t const& to) {
  ExtentOverrides out;
  for (auto const& [pos, w] : ov) {
    if (pos >= from.size()) continue;
    auto const it = std::find(to.begin(), to.end(), from[pos]);
    if (it != to.end())
      out.emplace(static_cast<std::size_t>(it - to.begin()), w);
  }
  return out;
}

// Project a value's positional overrides into LABEL space via its annotation:
// position p -> (annot[p] -> width). For the flops call, whose out/contracted
// index sets ARE annotation labels (see CostModel::flops).
[[nodiscard]] inline container::map<Index, std::size_t> extents_by_label(
    ExtentOverrides const& ov, annot_t const& annot) {
  container::map<Index, std::size_t> out;
  for (auto const& [pos, w] : ov)
    if (pos < annot.size()) out.emplace(annot[pos], w);
  return out;
}

// Uniform read access to a DryRun Result's (index list, overrides, cost
// model) regardless of which concrete DryRun type `r` is. `is<T>()`/`as<T>()`
// are public Result methods, so no friendship is needed; declared here (and
// defined below, after both concrete classes) purely because their bodies
// need the concrete classes' definitions.
[[nodiscard]] container::svector<Index> indices_of(Result const& r);
[[nodiscard]] ExtentOverrides overrides_of(Result const& r);
/// Positional lower bounds of SLICED modes (mode -> element lobound); a mode
/// absent here is whole (lobound 0). Kept parallel to \c ExtentOverrides so a
/// slice preserves its ABSOLUTE position, as the TA backend does.
[[nodiscard]] ExtentOverrides lobounds_of(Result const& r);

/// The realized element range of position \p pos of an operand: (lo, extent)
/// -- sliced modes from the overrides, whole modes from the regime.
/// (lo, extent); extent is 0 when unknown (whole mode whose label the regime
/// does not size, e.g. a fixture hyperindex) -- callers compare only known
/// extents.
[[nodiscard]] inline std::pair<std::size_t, std::size_t> range_at(
    ExtentOverrides const& ov, ExtentOverrides const& lob,
    std::shared_ptr<CostModel const> const& cm, annot_t const& annot,
    std::size_t pos) {
  std::size_t lo = 0, ext = 0;
  if (auto it = ov.find(pos); it != ov.end())
    ext = it->second;
  else
    try {
      ext = cm->regime().extent(annot[pos]);
    } catch (std::exception const&) {
      ext = 0;  // unknown
    }
  if (auto it = lob.find(pos); it != lob.end()) lo = it->second;
  return {lo, ext};
}

/// DRY-RUN CONFORMANCE CHECK (the dry-run analogue of TA's einsum index-map
/// merge and is_range_set_congruent): every label shared by the two operand
/// annotations must realize the SAME element range (lobound AND extent) on
/// both sides. A whole operand meeting a sliced partner (extent mismatch) or
/// two operands sliced to DIFFERENT batches of one loop (lobound mismatch,
/// e.g. a stale per-batch cell reused across batches) is a schedule/runtime
/// slicing defect; the wet backend hangs or asserts on it, the dry run used to
/// pass silently. Throws with both ranges named.
inline void check_shared_ranges(char const* op, annot_t const& lannot,
                                ExtentOverrides const& lov,
                                ExtentOverrides const& llob,
                                annot_t const& rannot,
                                ExtentOverrides const& rov,
                                ExtentOverrides const& rlob,
                                std::shared_ptr<CostModel const> const& cm) {
  for (std::size_t lp = 0; lp < lannot.size(); ++lp) {
    auto const rit = std::find(rannot.begin(), rannot.end(), lannot[lp]);
    if (rit == rannot.end()) continue;
    std::size_t const rp = static_cast<std::size_t>(rit - rannot.begin());
    // Both whole on this label => trivially equal (no regime query needed).
    if (!lov.count(lp) && !llob.count(lp) && !rov.count(rp) && !rlob.count(rp))
      continue;
    auto const [llo, lext] = range_at(lov, llob, cm, lannot, lp);
    auto const [rlo, rext] = range_at(rov, rlob, cm, rannot, rp);
    if (llo == rlo && (lext == rext || lext == 0 || rext == 0)) continue;
    auto lbls = [](annot_t const& a) {
      std::string s;
      for (auto const& x : a) s += toUtf8(x.full_label()) + " ";
      return s;
    };
    throw std::runtime_error(std::format(
        "[dryrun] {}: shared label {} realizes DIFFERENT ranges on the two "
        "operands: L[{}]=[{},{}) vs R[{}]=[{},{}) -- a whole operand against "
        "a sliced partner (extent) or two slices of different batches "
        "(lobound); L=({}) R=({})",
        op, toUtf8(lannot[lp].full_label()), lp, llo, llo + lext, rp, rlo,
        rlo + rext, lbls(lannot), lbls(rannot)));
  }
}

/// Accumulation (add_inplace) conformance: the accumulator and the partial
/// must realize the SAME range on every mode (same slicing, same batch) --
/// summing a partial of one batch into a cell of another is the stale-cell
/// signature. Positional (both carry the result's own index order).
inline void check_accumulate_ranges(
    container::svector<Index> const& idx, ExtentOverrides const& ov,
    ExtentOverrides const& lob, Result const& other,
    std::shared_ptr<CostModel const> const& cm) {
  auto const oov = overrides_of(other);
  auto const olob = lobounds_of(other);
  annot_t const a(idx.begin(), idx.end());
  for (std::size_t pos = 0; pos < idx.size(); ++pos) {
    if (!ov.count(pos) && !lob.count(pos) && !oov.count(pos) &&
        !olob.count(pos))
      continue;
    auto const [lo, ext] = range_at(ov, lob, cm, a, pos);
    auto const [olo, oext] = range_at(oov, olob, cm, a, pos);
    if (lo == olo && (ext == oext || ext == 0 || oext == 0)) continue;
    throw std::runtime_error(std::format(
        "[dryrun] add_inplace: mode {} ({}) accumulator range [{},{}) vs "
        "partial range [{},{}) -- accumulating a partial of a different "
        "slicing/batch into this cell",
        pos, toUtf8(idx[pos].full_label()), lo, lo + ext, olo, olo + oext));
  }
}

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
                                     std::array<std::any, 3> const& annot,
                                     ExtentOverrides const& lob = {}) {
    auto const a = Annot<annot_t>{annot};
    auto const other_ov = overrides_of(other);
    auto const other_lob = lobounds_of(other);
    check_shared_ranges("sum", a.lannot, ov, lob, a.rannot, other_ov, other_lob,
                        cm);
    // Both summands share the result's index set (possibly reordered); project
    // each operand's positional slice widths onto the result annotation by
    // label before merging (position k is a different mode in each operand).
    auto merged = merge_overrides(
        remap_overrides_by_annot(ov, a.lannot, a.this_annot),
        remap_overrides_by_annot(other_ov, a.rannot, a.this_annot));
    auto merged_lob = merge_overrides(
        remap_overrides_by_annot(lob, a.lannot, a.this_annot),
        remap_overrides_by_annot(other_lob, a.rannot, a.this_annot));
    return make_dryrun_result(
        container::svector<Index>(a.this_annot.begin(), a.this_annot.end()), cm,
        std::move(merged), std::move(merged_lob));
  }

  [[nodiscard]] static ResultPtr prod(
      container::svector<Index> const& idx, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm, Result const& other,
      std::array<std::any, 3> const& annot, ExtentOverrides const& lob = {}) {
    if (other.is<ResultScalar<double>>()) {
      // Scalar * tensor: shape (and any accumulated slicing) unchanged.
      return make_dryrun_result(idx, cm, ov, lob);
    }
    auto const a = Annot<annot_t>{annot};
    auto const other_ov = overrides_of(other);
    auto const other_lob = lobounds_of(other);
    check_shared_ranges("prod", a.lannot, ov, lob, a.rannot, other_ov,
                        other_lob, cm);
    auto merged_lob = merge_overrides(
        remap_overrides_by_annot(lob, a.lannot, a.this_annot),
        remap_overrides_by_annot(other_lob, a.rannot, a.this_annot));
    // RESULT overrides: each operand's positional slice widths are positional
    // against ITS OWN annotation (== its canon index order); project both onto
    // the result's annotation by label (dropping any contracted-away mode),
    // then merge. Merging the two raw positional maps directly would be wrong
    // -- position k means a different mode in each operand.
    auto merged = merge_overrides(
        remap_overrides_by_annot(ov, a.lannot, a.this_annot),
        remap_overrides_by_annot(other_ov, a.rannot, a.this_annot));

    // Emit the cost model's OWN flops / roofline exec_cost for THIS op into the
    // eval trace (gated on the eval log level), interleaved right before the
    // generic engine's `Eval | Product` line for the same op. This lets trace
    // post-processing weight avoidable recomputation by MODELLED TIME without
    // re-deriving the cost downstream (which would silently drift from the
    // model). The (out, contracted) index sets and the realized (sliced) extent
    // overrides feed the same CostModel closures the static cost_profile() walk
    // uses, so per-op costs are consistent with the whole-forest totals.
    if (Logger::instance().eval.level > 0) {
      // Cost THIS op in the einsum ANNOTATION label space (lannot/rannot/
      // this_annot), NOT the operands' stored `indices_` (idx / indices_of(
      // other)). A DAG VALUE HAS NO INTRINSIC LABELS -- only ops bind labels to
      // it, meaningful only within that op; a value's `indices_` merely holds
      // whatever labels its PRODUCER used, which say nothing about how a
      // CONSUMER binds it. The two diverge for every CSE-shared value used
      // under a different binding: e.g. the (g.C)(g.C) legs are the SAME cached
      // value (its `indices_` reads [i_1 i_4 K_2 a_4] from its producer) yet
      // THIS op binds it as lannot=[i_2 i_3 K_2 a_3] and rannot=[i_2 i_1 K_2
      // a_1]. Deriving `contracted` from the producer-labeled stored indices
      // instead of this op's annotations unions modes from different label
      // contexts and exploded the flops (6.65e16 vs the correct ~1.3e13). out U
      // (lannot & rannot) == lannot U rannot == the real contraction volume.
      container::svector<Index> out(a.this_annot.begin(), a.this_annot.end());
      container::svector<Index> contracted;
      for (auto const& ix : a.lannot)
        if (std::find(a.rannot.begin(), a.rannot.end(), ix) != a.rannot.end())
          contracted.push_back(ix);
      // Slice widths in LABEL space for the flops call: project each operand's
      // positional overrides through its annotation. A batched label shared by
      // both operands (e.g. a contracted, batched aux index) carries the same
      // width from either side, so the first insertion wins harmlessly.
      auto label_extents = extents_by_label(ov, a.lannot);
      for (auto const& [lbl, w] : extents_by_label(other_ov, a.rannot))
        label_extents.emplace(lbl, w);
      double const flops = cm->flops(out, contracted, label_extents);
      sequant::eval::detail::last_op_flops() = flops;  // for the Build event
      double const exec = cm->exec_cost(flops, cm->memsize(idx, ov), 4096);
      sequant::eval::detail::last_op_exec() = exec;  // for the Build event
      write_log(Logger::instance(), "OpCost", std::format(" | {}", flops),
                std::format(" | {}", exec), '\n');
      // Fold this op's SLICED-extent cost into the replay cost sink, if one is
      // attached (cost_profile()'s recompute-aware tally). merged/ov carry the
      // runtime slicing, so a contraction re-executed once per occ block is
      // charged once per block at its sliced size -- the same numbers already
      // logged above, now summed. No-op (byte-identical) when unattached.
      cm->tally_op(flops, exec);
      // last_op_flops (set above) is THIS build's ACTUAL realized-extent cost.
      // The eval loop's build choke point reads it and records it against the
      // node's IDENTITY, at the (value, SLICE) granularity, in the (root)
      // cache's recompute tally -- so avoidable recompute is the actual FLOPs a
      // slice was rebuilt beyond once, with no ill-defined "full extent"
      // denominator (slicing is non-uniform). prod cannot form the node
      // identity here (it has no node, and the include cycle
      // dryrun/eval_expr.hpp -> cost_model_object.hpp keeps the node type out
      // of the CostSink), so the rollup is done there. See
      // CacheManager::tally_build / recompute_tally().
    }

    if (a.this_annot.empty()) {
      // Full contraction -> scalar. No real numeric value is ever tracked by
      // this zero-data backend (only sizes/costs), so the placeholder 0.0
      // is never meant to be read as a physical result.
      return eval_result<ResultScalar<double>>(0.0);
    }
    return make_dryrun_result(
        container::svector<Index>(a.this_annot.begin(), a.this_annot.end()), cm,
        std::move(merged), std::move(merged_lob));
  }

  [[nodiscard]] static ResultPtr permute(
      container::svector<Index> const& idx, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm,
      std::array<std::any, 2> const& ann, ExtentOverrides const& lob = {}) {
    auto const post = std::any_cast<annot_t>(ann[1]);
    // Reordering modes moves each mode's position, so the positional overrides
    // must move with them: `ov` is positional against `idx` (the pre-permute
    // canon order); project it onto `post` by label.
    return make_dryrun_result(
        container::svector<Index>(post.begin(), post.end()), cm,
        remap_overrides_by_annot(ov, idx, post),
        remap_overrides_by_annot(lob, idx, post));
  }

  [[nodiscard]] static ResultPtr slice_mode(
      container::svector<Index> const& idx, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm, std::size_t mode,
      std::size_t elem_lo, std::size_t elem_hi,
      ExtentOverrides const& lob = {}) {
    SEQUANT_ASSERT(mode < idx.size());
    auto merged = ov;
    merged[mode] = elem_hi - elem_lo;  // positional: mode `mode`, any label
    auto merged_lob = lob;
    merged_lob[mode] = elem_lo;  // ABSOLUTE position preserved (as TA does)
    return make_dryrun_result(idx, cm, std::move(merged),
                              std::move(merged_lob));
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
    // mode index must equal the slice width it is being written into. The
    // block's overrides are positional against ITS OWN index list, so locate
    // the shared index there (its mode need not equal the dest's `mode`).
    auto const bov = overrides_of(block);
    auto const bidx = indices_of(block);
    std::size_t const block_extent = [&] {
      auto const it = std::find(bidx.begin(), bidx.end(), mix);
      if (it != bidx.end()) {
        auto const bpos = static_cast<std::size_t>(it - bidx.begin());
        if (auto ov_it = bov.find(bpos); ov_it != bov.end())
          return ov_it->second;
      }
      return cm->regime().extent(mix);
    }();
    if (block_extent != block_hi - block_lo) {
      std::cerr << "[scatter-mismatch] result idx=[";
      for (auto const& x : idx) std::cerr << toUtf8(x.full_label()) << " ";
      std::cerr << "] slice mode=" << mode << " (" << toUtf8(mix.full_label())
                << ") block idx=[";
      for (auto const& x : bidx) std::cerr << toUtf8(x.full_label()) << " ";
      std::cerr << "] block_extent=" << block_extent << " expected slice=["
                << block_lo << "," << block_hi << ")=" << (block_hi - block_lo)
                << " regime_extent(mix)=" << cm->regime().extent(mix)
                << " block_ov={";
      for (auto const& [bp, ex] : bov) std::cerr << bp << ":" << ex << " ";
      std::cerr << "}" << std::endl;
    }
    SEQUANT_ASSERT(block_extent == block_hi - block_lo);
    if (block_extent != block_hi - block_lo)
      throw std::runtime_error(std::format(
          "[dryrun] write_into_slice: block extent {} on mode {} ({}) != "
          "destination slice [{},{})",
          block_extent, mode, toUtf8(mix.full_label()), block_lo, block_hi));
    // The block's ABSOLUTE position on this mode (if it is a slice) must be
    // the destination slice it is written into.
    {
      auto const blob = lobounds_of(block);
      auto const bit = std::find(bidx.begin(), bidx.end(), mix);
      if (bit != bidx.end()) {
        auto const bpos = static_cast<std::size_t>(bit - bidx.begin());
        if (auto lit = blob.find(bpos);
            lit != blob.end() && lit->second != block_lo)
          throw std::runtime_error(std::format(
              "[dryrun] write_into_slice: block lobound {} on mode {} ({}) != "
              "destination slice lobound {}",
              lit->second, mode, toUtf8(mix.full_label()), block_lo));
      }
    }
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
    ov[mode] = lohi.second - lohi.first;  // positional: dest mode `mode`
  }

  [[nodiscard]] static container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(container::svector<Index> const& idx, ExtentOverrides const& ov,
               std::shared_ptr<CostModel const> const& cm, std::size_t mode,
               std::size_t target_batch_size) {
    SEQUANT_ASSERT(mode < idx.size());
    Index const& ix = idx[mode];  // for the regime extent / slice partition
    std::size_t extent;
    if (auto it = ov.find(mode); it != ov.end())  // positional override
      extent = it->second;
    else
      extent = cm->regime().extent(ix);

    container::svector<std::pair<std::size_t, std::size_t>> out;
    if (target_batch_size == 0 || extent == 0) {
      out.push_back({0, extent});
      return out;
    }
    // CALLER-SUPPLIED PARTITION: if the mode's space has a recorded batch
    // partition (SizeRegime::space_slice_extents), the wet backend slices this
    // axis along whole TILES (mode_batches_of_trange1 reads the operand's real
    // TiledRange1), so a batch boundary always falls on a tile edge and a
    // (sub)range spanning N whole tiles yields N batches -- NOT extent/target
    // uniform blocks. The partition slices ARE those target-grouped tile edges
    // (batch_slice_extents_from_tiles applied the target once, at harvest), so
    // we emit the PREFIX of partition slices that sums to `extent`:
    //   - outer call (ov absent, extent == full axis) => all slices, the full
    //     partition (e.g. aux 672 -> [168,168,168,168], 4 batches);
    //   - nested call (ov narrowed the axis to one outer batch, extent < full)
    //     => the prefix reaching that extent, so a single-tile sub-range (168)
    //     is ONE atomic batch, matching the wet backend, instead of being
    //     re-sliced into ceil(168/64)=3 uniform blocks (which then cascade).
    // `target_batch_size` is not re-applied here: the partition already encodes
    // it. If `extent` does not land on a partition boundary (not tile-aligned),
    // fall through to uniform blocks. The dry-run stays model-agnostic -- it
    // only reads slice extents; the caller decided the tiling.
    auto const& slices = cm->regime().slice_extents(ix);
    if (!slices.empty()) {
      std::size_t lo = 0;
      for (std::size_t const s : slices) {
        if (lo >= extent) break;
        out.push_back({lo, lo + s});
        lo += s;
      }
      if (lo == extent) return out;  // extent tile-aligned to the partition
      out.clear();                   // not aligned -> uniform fallback below
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
               ExtentOverrides overrides = {}, ExtentOverrides lobounds = {})
      : Result{Payload{}},
        indices_{std::move(idxset)},
        cm_{std::move(cm)},
        overrides_{std::move(overrides)},
        lobounds_{std::move(lobounds)} {}
  [[nodiscard]] ExtentOverrides const& lobounds() const noexcept {
    return lobounds_;
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
    return id_for_type<ResultDryRun>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    return detail::DryRunOps::sum(indices_, overrides_, cm_, other, annot,
                                  lobounds_);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    return detail::DryRunOps::prod(indices_, overrides_, cm_, other, annot,
                                   lobounds_);
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann,
                                      lobounds_);
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann,
                                      lobounds_);
  }

  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    return detail::DryRunOps::slice_mode(indices_, overrides_, cm_, mode,
                                         elem_lo, elem_hi, lobounds_);
  }

  void write_into_slice(Result const& block, std::size_t mode,
                        std::size_t block_lo, std::size_t block_hi) override {
    detail::DryRunOps::write_into_slice(indices_, overrides_, assembled_, cm_,
                                        block, mode, block_lo, block_hi);
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<ResultDryRun>() || other.is<ResultDryRunNested>());
    detail::check_accumulate_ranges(indices_, overrides_, lobounds_, other,
                                    cm_);
    overrides_ =
        detail::merge_overrides(overrides_, detail::overrides_of(other));
    lobounds_ = detail::merge_overrides(lobounds_, detail::lobounds_of(other));
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_, lobounds_);
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_, lobounds_);
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t /*factor*/) const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_, lobounds_);
  }

  [[nodiscard]] std::size_t size_in_bytes() const final {
    return cm_->memsize(indices_, overrides_);
  }

  container::svector<Index> indices_;
  std::shared_ptr<CostModel const> cm_;
  ExtentOverrides overrides_;
  ExtentOverrides lobounds_;  // positional lobounds of sliced modes
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
                     container::svector<Index> canon_order = {},
                     ExtentOverrides lobounds = {})
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
        overrides_{std::move(overrides)},
        lobounds_{std::move(lobounds)} {}
  [[nodiscard]] ExtentOverrides const& lobounds() const noexcept {
    return lobounds_;
  }

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
    return detail::DryRunOps::sum(indices_, overrides_, cm_, other, annot,
                                  lobounds_);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    return detail::DryRunOps::prod(indices_, overrides_, cm_, other, annot,
                                   lobounds_);
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann,
                                      lobounds_);
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann,
                                      lobounds_);
  }

  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    return detail::DryRunOps::slice_mode(indices_, overrides_, cm_, mode,
                                         elem_lo, elem_hi, lobounds_);
  }

  void write_into_slice(Result const& block, std::size_t mode,
                        std::size_t block_lo, std::size_t block_hi) override {
    detail::DryRunOps::write_into_slice(indices_, overrides_, assembled_, cm_,
                                        block, mode, block_lo, block_hi);
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<ResultDryRun>() || other.is<ResultDryRunNested>());
    detail::check_accumulate_ranges(indices_, overrides_, lobounds_, other,
                                    cm_);
    overrides_ =
        detail::merge_overrides(overrides_, detail::overrides_of(other));
    lobounds_ = detail::merge_overrides(lobounds_, detail::lobounds_of(other));
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_, lobounds_);
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_, lobounds_);
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t /*factor*/) const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_, lobounds_);
  }

  [[nodiscard]] std::size_t size_in_bytes() const final {
    return cm_->memsize(indices_, overrides_);
  }

  container::svector<Index> outer_;
  container::svector<Index> inner_;
  container::svector<Index> indices_;  // canon order; outer_++inner_ content
  std::shared_ptr<CostModel const> cm_;
  ExtentOverrides overrides_;
  ExtentOverrides lobounds_;  // positional lobounds of sliced modes
  AssembledCoverage assembled_;
};

[[nodiscard]] inline ResultPtr make_dryrun_result(
    container::svector<Index> idx, std::shared_ptr<CostModel const> cm,
    ExtentOverrides overrides, ExtentOverrides lobounds) {
  if (!detail::has_proto(idx))
    return eval_result<ResultDryRun>(std::move(idx), std::move(cm),
                                     std::move(overrides), std::move(lobounds));
  container::svector<Index> outer, inner;
  for (auto const& ix : idx)
    (ix.has_proto_indices() ? inner : outer).push_back(ix);
  return eval_result<ResultDryRunNested>(std::move(outer), std::move(inner),
                                         std::move(cm), std::move(overrides),
                                         std::move(idx), std::move(lobounds));
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
[[nodiscard]] inline ExtentOverrides lobounds_of(Result const& r) {
  if (r.is<ResultDryRun>()) return r.as<ResultDryRun>().lobounds();
  SEQUANT_ASSERT(r.is<ResultDryRunNested>());
  return r.as<ResultDryRunNested>().lobounds();
}

}  // namespace detail

/// \brief Build a \c BackendArrayOps for the dry-run cost backend.
///
/// \details The dry-run analogue of \c make_ta_array_ops: \c make_zeros returns
/// a full-extent dry-run token shaped by the descriptor (every mode at its
/// space's natural CostModel extent -- the same all-full result the old
/// \c DryRunOps::pre_sized_zeros_over_mode produced once the scatter axis was
/// widened), and \c axis_batches partitions an axis's FULL space extent exactly
/// as \c DryRunOps::mode_batches does for an unsliced, override-free
/// single-mode token (so the dry run realizes the same batch COUNT the wet
/// backend does, via the shared \c SizeRegime::space_slice_extents partition).
/// Install it on the eval cache (\c CacheManager::set_array_ops) for a dry-run
/// batched eval.
[[nodiscard]] inline BackendArrayOps make_dryrun_array_ops(
    std::shared_ptr<CostModel const> cm) {
  BackendArrayOps aops;
  aops.axis_batches = [cm](Index const& axis, std::size_t target_batch_size) {
    return detail::DryRunOps::mode_batches(container::svector<Index>{axis}, {},
                                           cm, /*mode=*/0, target_batch_size);
  };
  aops.make_zeros =
      [cm](container::vector<Index> const& descriptor) -> ResultPtr {
    return make_dryrun_result(
        container::svector<Index>(descriptor.begin(), descriptor.end()), cm);
  };
  return aops;
}

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP
