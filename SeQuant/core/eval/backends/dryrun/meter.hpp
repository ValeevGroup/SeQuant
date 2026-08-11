#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_METER_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_METER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/peak_monitor.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <string>
#include <unordered_map>
#include <vector>

namespace sequant::eval::dryrun {

///
/// \brief One value's build-vs-home fidelity witness for a \c MeterReport
/// (see \c assemble_report): how many times a distinct value was built (over
/// its whole recompute tally) versus WHERE it is homed and WHERE it is used,
/// read off the matching \c RichSchedule::ValueCell (looked up by hash).
///
struct HomeFidelity {
  std::string label;       ///< value signature "idx:hash" (idx disambiguates a
                           ///< 64-bit hash collision between distinct nodes;
                           ///< see cost_profile.hpp's sig_full_flops)
  std::size_t hash = 0;    ///< the value's EvalExpr::hash_value()
  std::size_t builds = 0;  ///< total builds across slices (recompute-aware)
  std::string home;        ///< dag-scope of home_modes ("" == root {} -- a
                           ///< whole-nest invariant)
  std::string uses;        ///< dag-scope list of the value's occurrences
};

///
/// \brief Summary of one metered dry-run (or wet) replay: the hierarchy-wide
/// peak (from a wired \c PeakMonitor), the persistent/volatile FLOPs and
/// CostModel exec-time split (rolled up from a \c CacheManager's recompute
/// tally, classified by \c compute_volatility), the total build count, and
/// the per-value build-vs-home fidelity list (see \c HomeFidelity).
///
struct MeterReport {
  double peak_bytes = 0;         ///< PeakMonitor high-water (dense on dry
                                 ///< run / sparse on wet run)
  std::size_t peak_op_hash = 0;  ///< location (op hash) of the peak

  double flops_persistent = 0, flops_volatile = 0;  ///< dense model; dry-only
  double cost_persistent = 0, cost_volatile = 0;  ///< CostModel exec; dry-only

  std::size_t builds_total = 0;

  std::vector<HomeFidelity> home_fidelity;  ///< sorted: builds desc

  bool whole_scope = false;  ///< which executor this report describes
};

///
/// \brief Bottom-up memoized volatility over an evaluation \p forest: a node
/// is volatile iff \p is_volatile flags it directly, or (for an internal
/// node) either child is volatile -- the SAME rule the gated
/// \c sequant::cache_manager factory applies while building its NV/V
/// frontier (see cache_manager.hpp's DAG walk). Keyed by
/// \c TreeNode::hash_value() (rather than the node identity itself) so a
/// caller can classify a \c CacheManager::recompute_tally() entry -- keyed by
/// the SAME node identity but not necessarily the SAME node object -- by its
/// hash.
///
/// \param forest the evaluation forest (a range of eval nodes).
/// \param is_volatile `bool(TreeNode const&)`: true if the node is
///        intrinsically volatile. Only its value on leaves matters in
///        practice (volatility propagates up), but it is consulted on every
///        node, matching \c cache_manager's gated factory.
/// \return a map from \c hash_value() to whether that value is volatile.
///
template <class Forest, class IsVolatile>
std::unordered_map<std::size_t, bool> compute_volatility(
    Forest const& forest, IsVolatile const& is_volatile) {
  using Node = std::ranges::range_value_t<Forest>;

  std::unordered_map<std::size_t, bool> volatile_of;

  auto visit = [&](auto&& self, Node const& n) -> bool {
    std::size_t const h = n->hash_value();
    if (auto it = volatile_of.find(h); it != volatile_of.end())
      return it->second;
    bool v;
    if (n.leaf()) {
      v = is_volatile(n);
    } else {
      bool const vl = self(self, n.left());
      bool const vr = self(self, n.right());
      v = is_volatile(n) || vl || vr;
    }
    volatile_of.emplace(h, v);
    return v;
  };

  for (auto const& tree : forest) visit(visit, tree);
  return volatile_of;
}

///
/// \brief Assemble a \c MeterReport from a metered replay: a walked
/// \p cache (its \c recompute_tally() populated by \c CacheManager::
/// tally_build over the replay), the hierarchy-wide \p mon (\c PeakMonitor),
/// the \p rich linearized schedule (\c compute_dag_boulevard over the SAME
/// \p forest, supplying each value's home/use dag-scope), and \p is_volatile
/// (fed to \c compute_volatility to classify each distinct value).
///
/// Per distinct value (one \c cache.recompute_tally() entry): \c builds is
/// the sum, over its slices, of each slice's build count; \c node_flops /
/// \c node_exec are the sum, over its slices, of build-count times that
/// slice's actual (flops, exec). The value is classified persistent/volatile
/// by \c compute_volatility's verdict for its hash and folded into the
/// matching \c MeterReport::flops_*/cost_* accumulator. Its \c HomeFidelity
/// entry's \c home/uses are read off the \p rich cell sharing its hash (empty
/// if the value has no matching cell, e.g. a leaf never realized as its own
/// distinct product build).
///
/// \param cache the (root) cache whose \c recompute_tally() was populated by
///        a \c Trace::On metered replay.
/// \param mon the \c PeakMonitor wired onto \p cache's scope chain during the
///        replay.
/// \param rich the linearized schedule (\c compute_dag_boulevard) over the
///        SAME forest the replay walked.
/// \param forest the evaluation forest (fed to \c compute_volatility).
/// \param is_volatile `bool(TreeNode const&)`: intrinsic volatility
///        predicate, as for \c compute_volatility.
/// \param whole_scope which executor this report describes (stashed verbatim
///        into \c MeterReport::whole_scope).
/// \return the assembled \c MeterReport.
///
template <class Cache, class Forest, class IsVolatile>
MeterReport assemble_report(Cache const& cache, PeakMonitor const& mon,
                            RichSchedule const& rich, Forest const& forest,
                            IsVolatile const& is_volatile, bool whole_scope) {
  MeterReport report;
  report.whole_scope = whole_scope;
  report.peak_bytes = static_cast<double>(mon.hwmark_bytes);
  report.peak_op_hash = mon.peak.op_hash;

  auto const volatility = compute_volatility(forest, is_volatile);

  // hash -> ValueCell* lookup, mirroring make_node_meta's map build
  // (scope_executor.hpp): rich.cells is a flat vector, not keyed by hash.
  std::unordered_map<std::size_t, ValueCell const*> cell_by_hash;
  cell_by_hash.reserve(rich.cells.size());
  for (auto const& cell : rich.cells) cell_by_hash.emplace(cell.hash, &cell);

  // dag-scope formatting: comma-joined IndexSpace base_keys, no trailing
  // comma -- the same convention make_node_meta uses (scope_executor.hpp).
  auto const dag_scope = [](auto const& modes) {
    std::string s;
    for (auto const& m : modes) {
      if (!s.empty()) s += ",";
      s += toUtf8(m.space().base_key());
    }
    return s;
  };

  std::size_t idx = 0;
  for (auto const& [node, tally] : cache.recompute_tally()) {
    std::size_t builds = 0;
    double node_flops = 0.0, node_exec = 0.0;
    for (auto const& [sig, rec] : tally.slices) {
      builds += rec.count;
      node_flops += static_cast<double>(rec.count) * rec.flops;
      node_exec += static_cast<double>(rec.count) * rec.exec;
    }
    report.builds_total += builds;

    std::size_t const hash = node->hash_value();
    bool const is_vol = [&] {
      auto it = volatility.find(hash);
      return it != volatility.end() && it->second;
    }();

    if (is_vol) {
      report.flops_volatile += node_flops;
      report.cost_volatile += node_exec;
    } else {
      report.flops_persistent += node_flops;
      report.cost_persistent += node_exec;
    }

    HomeFidelity hf;
    hf.label = std::to_string(idx++) + ":" + std::to_string(hash);
    hf.hash = hash;
    hf.builds = builds;
    if (auto it = cell_by_hash.find(hash); it != cell_by_hash.end()) {
      auto const* cell = it->second;
      hf.home = dag_scope(cell->home_modes);
      container::svector<Index> uses_modes;
      for (auto const& occ : cell->occurrences)
        for (auto const& [mode, range] : occ.ectx) uses_modes.push_back(mode);
      hf.uses = dag_scope(uses_modes);
    }
    report.home_fidelity.push_back(std::move(hf));
  }

  std::sort(report.home_fidelity.begin(), report.home_fidelity.end(),
            [](HomeFidelity const& a, HomeFidelity const& b) {
              return a.builds > b.builds;
            });

  return report;
}

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_METER_HPP
