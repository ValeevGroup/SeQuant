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
  std::array<double, 5> csv_pno_moment{1.0, 1.0, 1.0, 1.0, 1.0};
  std::array<double, 5> csv_osv_moment{1.0, 1.0, 1.0, 1.0, 1.0};

  /// \return the flat extent of \p ix's space; throws \c std::out_of_range
  ///         if the space is not present in \c space_extent (fail loud rather
  ///         than silently defaulting to 1).
  [[nodiscard]] std::size_t extent(Index const& ix) const {
    return space_extent.at(std::wstring{ix.space().base_key()});
  }

  /// \return the k-th power-mean moment for a proto-indexed CSV/PNO composite
  ///         index (\p k clamped to 0..4), or \c pow(extent, k) for a plain
  ///         (non-composite) index. Rank is determined by the number of proto
  ///         indices: 2 => PNO (occupied pair), 1 => OSV (occupied single).
  [[nodiscard]] double inner_pow(Index const& composite, std::size_t k) const {
    if (k > 4) k = 4;
    auto const& protos = composite.proto_indices();
    if (protos.empty())
      return std::pow(static_cast<double>(extent(composite)),
                      static_cast<double>(k));
    return (protos.size() >= 2) ? csv_pno_moment[k] : csv_osv_moment[k];
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
