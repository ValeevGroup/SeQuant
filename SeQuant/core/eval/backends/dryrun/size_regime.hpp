#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP

#include <SeQuant/core/index.hpp>

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

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP
