#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <string>

namespace sequant::eval::dryrun {

/// Per-basis-key extents and per-rank CSV moment tables that define one size
/// regime for a dry-run replay. Extents are element counts; CSV moments are
/// power means over occupied pairs (PNO) or singles (OSV).
/// Largest group of composites sharing one cluster the moment tables price:
/// a table carries the occupancy-folded per-member factors for k in
/// [1, kMaxCompositeGroup] (index 0 unused). A larger group is an error, not
/// an extrapolation -- the producer measures higher moments cheaply, so a
/// table that stops short is a configuration bug, and silently reusing the
/// last column under-prices the group by p^((k - K)/K).
inline constexpr std::size_t kMaxCompositeGroup = 8;
/// One moment table: index k holds the per-member factor of a k-group.
using MomentTable = std::array<double, kMaxCompositeGroup + 1>;
/// The dense-domain, unit-size table (every entry 1).
inline constexpr MomentTable unit_moment_table() {
  MomentTable t{};
  for (auto& x : t) x = 1.0;
  return t;
}

struct SizeRegime {
  /// element extent per IndexBasis::base_key(); an unnamed basis instance
  /// shares its space's entry (see IndexBasis)
  std::map<std::wstring, std::size_t> space_extent;

  /// Optional per-key batch partition: the element extent of each realized
  /// batch slice along the key's batch axis, keyed by IndexBasis::base_key(),
  /// in order. Empty (default) => the dry-run batches a mode into uniform
  /// target_batch_size blocks (backend-model-agnostic fallback). When present
  /// for a batch axis, ResultDryRun::mode_batches uses this partition directly
  /// (accumulated to [lo,hi) ranges), so the dry-run's batch count -- hence its
  /// recompute -- matches whatever the wet backend realizes, even when a tile
  /// is coarser than target_batch_size.
  ///
  /// The dry-run backend deliberately does not know how these were derived: the
  /// caller converts its backend's structure into slice extents and supplies
  /// them here, so a new backend model is supported without touching dry-run
  /// eval internals. For a tile-based wet backend, \c batch_slice_extents_from_
  /// tiles is the ready-made converter. The extents must sum to space_extent.
  std::map<std::wstring, container::svector<std::size_t>> space_slice_extents;

  // csv_pno_moment[k] / csv_osv_moment[k] hold the k-th power mean
  // M_k = (mean_over_pairs d^k)^(1/k) of the per-pair PNO / per-orbital OSV
  // domain size d, for k in [1, kMaxCompositeGroup] (index 0 is unused, set to
  // 1). inner_pow() returns M_k so that inner_aware_volume's per-member product
  // over a k-composite group is M_k^k = mean(d^k), and outer_nocc^N * M_k^k
  // equals the true block-sparse volume Sum_pairs d^k. Do not store raw moments
  // mean(d^k) here: that would over-count k-composite groups by a further
  // power of k. For a constant domain d, M_k = d for all k.
  MomentTable csv_pno_moment = unit_moment_table();
  MomentTable csv_osv_moment = unit_moment_table();

  // Moment tables for CSV cluster ranks >= 3 (CSV-CCSDT triples and beyond),
  // keyed by cluster rank (= number of proto indices). csv_moment_by_rank[r][k]
  // is the k-th power mean of the rank-r cluster domain. A rank not present
  // falls back to csv_pno_moment (the rank-2 table) in inner_pow(), so a
  // proto-rank >= 2 with no table of its own is sized by the PNO table.
  // Ranks 1 and 2 are held by csv_osv_moment / csv_pno_moment above and
  // are not expected here (an entry for 1 or 2 is ignored by inner_pow()).
  std::map<std::size_t, MomentTable> csv_moment_by_rank;

  /// The CSV moment tables of one basis instance (see csv_moment_by_instance):
  /// the rank-1 (OSV), rank-2 (PNO) and rank >= 3 tables, laid out and
  /// dispatched exactly like csv_osv_moment / csv_pno_moment /
  /// csv_moment_by_rank.
  struct CsvMomentTables {
    MomentTable osv = unit_moment_table();
    MomentTable pno = unit_moment_table();
    std::map<std::size_t, MomentTable> by_rank;

    /// \return the k-th power mean of the rank-\p rank cluster domain: osv
    ///         for rank 1, pno for rank 2, by_rank[rank] for a rank >= 3 that
    ///         has a table, else pno
    [[nodiscard]] double moment(std::size_t rank, std::size_t k) const {
      if (rank <= 1) return osv[k];
      if (rank == 2) return pno[k];
      auto const it = by_rank.find(rank);
      return it != by_rank.end() ? it->second[k] : pno[k];
    }
  };

  /// Moment tables keyed by the basis instance the composite's IndexBasis
  /// carries (IndexBasis::basis_instance()). A composite whose index carries
  /// an instance listed here is sized by that instance's tables; every other
  /// composite -- no instance, or one not listed -- by csv_osv_moment /
  /// csv_pno_moment / csv_moment_by_rank, which thus double as the default.
  /// A consumer whose one evaluation spans several CSV bases (e.g. the
  /// ground-state amplitudes expanded in one PNO basis and the response
  /// amplitudes in another, each family's indices carrying its basis
  /// instance) files every basis here, so no family is priced by another's
  /// domain sizes. Empty (default) => one basis, the default tables.
  std::map<IndexBasis::instance_type, CsvMomentTables> csv_moment_by_instance;

  /// \return the flat extent of the proto-free index \p ix's
  ///         \c IndexBasis::base_key(); throws \c std::out_of_range if the
  ///         key is not present in \c space_extent (fail loud rather than
  ///         silently defaulting to 1), and \c sequant::Exception for a
  ///         proto-indexed composite: a composite has no integer extent (its
  ///         base key names the space it is expanded in, not its size) and is
  ///         sized by \ref inner_pow alone.
  [[nodiscard]] std::size_t extent(Index const& ix) const {
    if (ix.has_proto_indices())
      throw Exception(
          "SizeRegime::extent: asked for the integer extent of the "
          "proto-indexed (CSV/PNO) composite " +
          toUtf8(ix.full_label()) +
          " -- composites are sized by inner_pow, not by base-space extent");
    return space_extent.at(std::wstring{ix.basis().base_key()});
  }

  /// \return \p ix's \c IndexBasis::base_key() batch-slice-extent sequence, or
  /// an
  ///         empty span if the key has no partition recorded (=> the caller
  ///         falls back to uniform target_batch_size blocks). Never throws.
  [[nodiscard]] container::svector<std::size_t> const& slice_extents(
      Index const& ix) const {
    static const container::svector<std::size_t> empty;
    auto const it =
        space_slice_extents.find(std::wstring{ix.basis().base_key()});
    return it != space_slice_extents.end() ? it->second : empty;
  }

  /// \return the k-th power-mean moment for a proto-indexed CSV/PNO composite
  ///         index (\p k in [1, kMaxCompositeGroup], else throws), or \c
  ///         pow(extent, k) for a plain (non-composite) index. The tables are
  ///         the composite's basis instance's if csv_moment_by_instance lists
  ///         it, else the default ones. Rank is determined by the number of
  ///         proto indices: 1 => OSV (occupied single), 2 => PNO (occupied
  ///         pair), >= 3 => the rank-specific table if present, else the PNO
  ///         (rank-2) table.
  [[nodiscard]] double inner_pow(Index const& composite, std::size_t k) const {
    auto const& protos = composite.proto_indices();
    if (!protos.empty() && (k < 1 || k > kMaxCompositeGroup))
      throw Exception(
          "SizeRegime::inner_pow: a group of " + std::to_string(k) +
          " composites on one cluster exceeds the moment table (k in [1, " +
          std::to_string(kMaxCompositeGroup) +
          "]); measure higher moments or treat the term as ill-formed");
    if (protos.empty())
      return std::pow(static_cast<double>(extent(composite)),
                      static_cast<double>(k));
    auto const rank = protos.size();
    if (!csv_moment_by_instance.empty())
      if (auto const& instance = composite.basis().basis_instance())
        if (auto const it = csv_moment_by_instance.find(*instance);
            it != csv_moment_by_instance.end())
          return it->second.moment(rank, k);
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

/// Convert a tile-extent sequence into a batch-slice-extent sequence
/// (SizeRegime::space_slice_extents) by the same whole-tile grouping the wet
/// batched evaluator uses (mode_batches_of_trange1, tiledarray/result.hpp):
/// accumulate consecutive tiles into a slice until appending the next would
/// push the slice over \p target_batch_size, then start a new slice; a lone
/// tile larger than the target still forms its own slice. Slice boundaries fall
/// on tile edges. This is a convenience converter for a tile-based caller; the
/// dry-run backend never calls it -- it only reads the resulting slice extents,
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
