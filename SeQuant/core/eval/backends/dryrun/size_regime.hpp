#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <string>

namespace sequant::eval::dryrun {

/// Per-space extents and per-rank CSV moment tables that define one size
/// regime for a dry-run replay. Extents are element counts; CSV moments are
/// power means over occupied pairs (PNO) or singles (OSV).
struct SizeRegime {
  std::map<std::wstring, std::size_t> space_extent;

  /// OPTIONAL per-space BATCH PARTITION: the element extent of each realized
  /// batch slice along the space's batch axis, keyed by space base_key, in
  /// order. Empty (default) => the dry-run batches a mode into UNIFORM
  /// target_batch_size blocks (backend-model-agnostic fallback). When present
  /// for a batch axis, ResultDryRun::mode_batches uses THIS partition directly
  /// (accumulated to [lo,hi) ranges), so the dry-run's batch COUNT -- hence its
  /// recompute -- matches whatever the wet backend realizes, even when a tile
  /// is coarser than target_batch_size.
  ///
  /// The dry-run backend deliberately does NOT know how these were derived: the
  /// CALLER converts its backend's structure into slice extents and supplies
  /// them here, so a new backend model is supported without touching dry-run
  /// eval internals. For a TILE-based wet backend, \c batch_slice_extents_from_
  /// tiles is the ready-made converter. The extents must sum to space_extent.
  std::map<std::wstring, container::svector<std::size_t>> space_slice_extents;

  // csv_pno_moment[k] / csv_osv_moment[k] hold the k-th POWER MEAN
  // M_k = (mean_over_pairs d^k)^(1/k) of the per-pair PNO / per-orbital OSV
  // domain size d, for k in [1,4] (index 0 is unused, set to 1). inner_pow()
  // returns M_k so that inner_aware_volume's per-member product over a
  // k-composite group is M_k^k = mean(d^k), and outer_nocc^N * M_k^k equals
  // the true block-sparse volume Sum_pairs d^k. Do NOT store raw moments
  // mean(d^k) here: that would over-count k-composite groups by a further
  // power of k. For a constant domain d, M_k = d for all k.
  std::array<double, 5> csv_pno_moment{1.0, 1.0, 1.0, 1.0, 1.0};
  std::array<double, 5> csv_osv_moment{1.0, 1.0, 1.0, 1.0, 1.0};

  // Moment tables for CSV cluster ranks >= 3 (CSV-CCSDT triples and beyond),
  // keyed by cluster rank (= number of proto indices). csv_moment_by_rank[r][k]
  // is the k-th power mean of the rank-r cluster domain. A rank not present
  // falls back to csv_pno_moment (the rank-2 table) in inner_pow(), preserving
  // the pre-rank-general behavior where every proto-rank >= 2 used the PNO
  // table. Ranks 1 and 2 are held by csv_osv_moment / csv_pno_moment above and
  // are NOT expected here (an entry for 1 or 2 is ignored by inner_pow()).
  std::map<std::size_t, std::array<double, 5>> csv_moment_by_rank;

  /// \return the flat extent of \p ix's space; throws \c std::out_of_range
  ///         if the space is not present in \c space_extent (fail loud rather
  ///         than silently defaulting to 1).
  [[nodiscard]] std::size_t extent(Index const& ix) const {
    return space_extent.at(std::wstring{ix.space().base_key()});
  }

  /// \return \p ix's space batch-slice-extent sequence, or an empty span if the
  ///         space has no partition recorded (=> the caller falls back to
  ///         uniform target_batch_size blocks). Never throws.
  [[nodiscard]] container::svector<std::size_t> const& slice_extents(
      Index const& ix) const {
    static const container::svector<std::size_t> empty;
    auto const it =
        space_slice_extents.find(std::wstring{ix.space().base_key()});
    return it != space_slice_extents.end() ? it->second : empty;
  }

  /// \return the k-th power-mean moment for a proto-indexed CSV/PNO composite
  ///         index (\p k clamped to 0..4), or \c pow(extent, k) for a plain
  ///         (non-composite) index. Rank is determined by the number of proto
  ///         indices: 1 => OSV (occupied single), 2 => PNO (occupied pair),
  ///         >= 3 => the rank-specific csv_moment_by_rank table if present,
  ///         else the PNO (rank-2) table.
  [[nodiscard]] double inner_pow(Index const& composite, std::size_t k) const {
    if (k > 4) k = 4;
    auto const& protos = composite.proto_indices();
    if (protos.empty())
      return std::pow(static_cast<double>(extent(composite)),
                      static_cast<double>(k));
    auto const rank = protos.size();
    if (rank <= 1) return csv_osv_moment[k];
    if (rank == 2) return csv_pno_moment[k];
    auto const it = csv_moment_by_rank.find(rank);
    return (it != csv_moment_by_rank.end()) ? it->second[k] : csv_pno_moment[k];
  }

  [[nodiscard]] std::function<std::size_t(Index const&)> idx_to_extent() const {
    return [this](Index const& ix) { return extent(ix); };
  }

  [[nodiscard]] std::function<double(Index const&, std::size_t)> inner_pow_fn()
      const {
    return [this](Index const& ix, std::size_t k) { return inner_pow(ix, k); };
  }
};

/// Convert a TILE-extent sequence into a BATCH-slice-extent sequence
/// (SizeRegime::space_slice_extents) by the SAME whole-tile grouping the wet
/// batched evaluator uses (mode_batches_of_trange1, tiledarray/result.hpp):
/// accumulate consecutive tiles into a slice until appending the next would
/// push the slice over \p target_batch_size, then start a new slice; a lone
/// tile larger than the target still forms its own slice. Slice boundaries fall
/// on tile edges. This is a convenience converter for a TILE-based caller; the
/// dry-run backend never calls it -- it only READS the resulting slice extents,
/// so any other backend model can populate space_slice_extents differently
/// without touching dry-run internals.
[[nodiscard]] inline container::svector<std::size_t>
batch_slice_extents_from_tiles(
    container::svector<std::size_t> const& tile_extents,
    std::size_t target_batch_size) {
  container::svector<std::size_t> slices;
  std::size_t const target = std::max<std::size_t>(target_batch_size, 1);
  std::size_t acc = 0;
  for (std::size_t const tsz : tile_extents) {
    if (acc > 0 && acc + tsz > target) {
      slices.push_back(acc);
      acc = 0;
    }
    acc += tsz;
  }
  if (acc > 0) slices.push_back(acc);
  return slices;
}

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP
