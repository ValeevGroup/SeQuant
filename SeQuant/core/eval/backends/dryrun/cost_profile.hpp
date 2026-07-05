#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_PROFILE_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_PROFILE_HPP

#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>

#include <cstddef>
#include <functional>
#include <vector>

namespace sequant::eval::dryrun {

/// Configuration for a faithful (gated) dry-run cache: the same footprint gate
/// and free-batchable-axis veto the real batched eval loop applies, so a
/// free-batchable giant (a `mu~`/`K`-carrying DF intermediate) is NOT cached
/// whole but recomputed sliced under each consumer's batch trigger.
///
/// The element types mirror the gated \c sequant::cache_manager overload
/// (\c cache_manager.hpp): \c is_volatile is invoked on every \c TreeNode
/// (deduced as \c EvalNodeDryRun for the dry-run backend), \c
/// is_batchable_index on every result \c Index.
///
/// This struct lives here (rather than only in the test) because Task 4's
/// \c cost_profile() entry point consumes it.
struct CacheConfig {
  /// Footprint gate (bytes): a node whose result footprint exceeds this is not
  /// cached. 0 (default) disables the gate.
  double max_footprint = 0.;
  /// Minimum non-persistent repeats to cache an internal node (CSE rule).
  std::size_t min_repeats = 2;
  /// `bool(EvalNodeDryRun const&)`: true if the node is intrinsically volatile
  /// (typically the amplitude leaves). Empty => nothing is volatile.
  std::function<bool(EvalNodeDryRun const&)> is_volatile;
  /// `bool(Index const&)`: true for an index the runtime batched evaluator
  /// slices over (e.g. DF aux `K` / PAO `mu~`). A node whose result carries
  /// such a free index is vetoed from caching. Empty => nothing is batchable.
  std::function<bool(Index const&)> is_batchable_index;
};

/// Builds a gated dry-run cache from an eval-node range, a \p cfg, and a
/// \p regime that supplies the moment-aware node-size model used for the
/// footprint gate.
///
/// The footprint functor sizes a node's result (its \c canon_indices()) with
/// the SAME moment-aware counter the DryRun \c Result uses
/// (\c memsize_counter over \c regime.idx_to_extent()/inner_pow_fn()), scaled
/// to bytes, so the gate compares like-for-like against \c cfg.max_footprint.
///
/// Unlike the SIMPLE \c cache_manager(nodes) factory the ad-hoc dry-run test
/// sites use, this routes through the GATED overload so free-batchable-axis
/// giants are vetoed (matching the real run). Call \c CacheManager::reset() on
/// the returned cache between summands to drop per-term non-persistent scratch
/// while keeping persistent (cross-term) entries.
///
/// \param nodes the evaluation forest (a range of \c EvalNodeDryRun).
/// \param cfg   footprint/repeat/volatility/batchability configuration.
/// \param regime the size regime supplying extents and CSV moment tables.
/// \return a \c CacheManager over \c EvalNodeDryRun.
template <typename NodeRange>
auto build_dryrun_cache(NodeRange const& nodes, CacheConfig const& cfg,
                        SizeRegime const& regime) {
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());

  // Footprint (bytes) of a node's RESULT: canon_indices() fed to the
  // moment-aware counter (as the counter's `result` slot; the empty lhs/rhs
  // contribute nothing) times 8 bytes/element. Same arithmetic as the DryRun
  // Result::size_in_bytes(), so the gate is faithful.
  auto footprint_of =
      [memsize = std::move(memsize)](EvalNodeDryRun const& n) -> double {
    std::vector<Index> const result(n->canon_indices().begin(),
                                    n->canon_indices().end());
    return memsize(std::vector<Index>{}, std::vector<Index>{}, result) * 8.0;
  };

  // Default the predicates so the gated factory never invokes an empty
  // std::function (nothing volatile / nothing batchable leaves those gates
  // inert, matching the factory's own defaults).
  std::function<bool(EvalNodeDryRun const&)> is_volatile =
      cfg.is_volatile ? cfg.is_volatile
                      : std::function<bool(EvalNodeDryRun const&)>(
                            [](EvalNodeDryRun const&) { return false; });
  std::function<bool(Index const&)> is_batchable_index =
      cfg.is_batchable_index ? cfg.is_batchable_index
                             : std::function<bool(Index const&)>(
                                   [](Index const&) { return false; });

  return sequant::cache_manager(nodes, std::move(is_volatile), cfg.min_repeats,
                                std::move(footprint_of), cfg.max_footprint,
                                std::move(is_batchable_index));
}

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_PROFILE_HPP
