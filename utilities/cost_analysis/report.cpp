#include "report.hpp"

#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/version.hpp>

#include <algorithm>
#include <chrono>
#include <format>
#include <map>

using namespace sequant;

namespace {
constexpr double bytes_per_elem = 8.0;

double mb(const AsyCost& c) {
  return c.ops() * bytes_per_elem / (1024.0 * 1024.0);
}

void report_largest(const CellResult& cell, std::size_t top_n,
                    std::ostream& out) {
  std::vector<std::pair<const TreeNode*, const Record*>> rows;
  for (const auto& [node, rec] : cell.catalog)
    if (rec.memory != AsyCost::zero()) rows.emplace_back(&node, &rec);
  // Sort by size (largest first); tie-break on the construction string so the
  // row order is deterministic regardless of unordered_map iteration order.
  std::ranges::sort(rows, [](const auto& a, const auto& b) {
    const double sa = mb(a.second->memory), sb = mb(b.second->memory);
    if (sa != sb) return sb < sa;
    return full_expr(*a.first) < full_expr(*b.first);
  });
  out << "| Rank | Result | Spaces | Memory (MB) | Uses | Construction | "
         "Local FLOPs |\n|---|---|---|---|---|---|---|\n";
  for (std::size_t i = 0; i < rows.size() && i < top_n; ++i) {
    const auto& [node, rec] = rows[i];
    out << std::format("| {} | {} | {} | {:.4f} | {} | {} | {} |\n", i + 1,
                       rec->label, rec->spaces, mb(rec->memory), rec->uses,
                       full_expr(*node), rec->flops.text());
  }
  out << "\n";
}

void report_expensive(const CellResult& cell, std::size_t top_n,
                      std::ostream& out) {
  std::vector<std::pair<const TreeNode*, const Record*>> rows;
  for (const auto& [node, rec] : cell.catalog)
    if (rec.flops != AsyCost::zero()) rows.emplace_back(&node, &rec);
  // Most FLOPs first; tie-break on the construction string for determinism.
  std::ranges::sort(rows, [](const auto& a, const auto& b) {
    const double fa = a.second->flops.ops(), fb = b.second->flops.ops();
    if (fa != fb) return fb < fa;
    return full_expr(*a.first) < full_expr(*b.first);
  });
  out << "| Rank | Result | Spaces | FLOPs | Uses | Construction | "
         "Cost |\n|---|---|---|---|---|---|---|\n";
  for (std::size_t i = 0; i < rows.size() && i < top_n; ++i) {
    const auto& [node, rec] = rows[i];
    out << std::format("| {} | {} | {} | {:.4e} | {} | {} | {} |\n", i + 1,
                       rec->label, rec->spaces, rec->flops.ops(), rec->uses,
                       full_expr(*node), rec->flops.text());
  }
  out << "\n";
}

void report_shape_histogram(const CellResult& cell, std::ostream& out) {
  struct Agg {
    std::size_t count = 0, uses = 0;
    AsyCost mem;
  };
  std::map<std::string, Agg> by_shape;
  for (const auto& [_, rec] : cell.catalog) {
    if (rec.memory == AsyCost::zero()) continue;
    auto& g = by_shape[rec.spaces];
    ++g.count;
    g.uses += rec.uses;
    // Representative size of the shape = its largest instance. A single O/V/X
    // signature can group differing sizes only when multiple aux spaces (Κ, L)
    // are registered; taking the max keeps this deterministic either way.
    if (mb(g.mem) < mb(rec.memory)) g.mem = rec.memory;
  }
  std::vector<std::pair<std::string, Agg>> rows(by_shape.begin(),
                                                by_shape.end());
  std::ranges::sort(rows, [](const auto& a, const auto& b) {
    const double sa = mb(a.second.mem), sb = mb(b.second.mem);
    if (sa != sb) return sb < sa;
    return a.first < b.first;  // tie-break on shape string
  });
  out << "| Shape | Size (MB) | # Distinct | Total uses |\n|---|---|---|---|\n";
  for (const auto& [shape, g] : rows)
    out << std::format("| {} | {:.4f} | {} | {} |\n", shape, mb(g.mem), g.count,
                       g.uses);
  out << "\n";
}

}  // namespace

std::string full_expr(const TreeNode& n) {
  if (n.leaf()) return n->label();
  if (n->op_type() == EvalOp::Adjoint) return full_expr(n.left()) + "+";
  const char* op = n->op_type() == EvalOp::Product ? " * " : " + ";
  return "(" + full_expr(n.left()) + op + full_expr(n.right()) + ")";
}

void write_report(
    const Config& cfg,
    const std::vector<std::pair<std::string, CellResult>>& results,
    const SimResult& sim, std::ostream& out) {
  const auto now = std::chrono::floor<std::chrono::seconds>(
      std::chrono::system_clock::now());
  out << "# SeQuant cost analysis\n\n"
      << std::format("_Generated: {:%Y-%m-%d %H:%M:%S} UTC_  \n", now)
      << std::format("_SeQuant: {}_\n\n", git_revision())
      << std::format("Sizes in MB ({} bytes/element).\n\n", bytes_per_elem);
  for (const auto& [name, cell] : results) {
    out << std::format("## {}\n\n", name);
    out << std::format(
        "Terms: {}; distinct intermediates: {}; reused: {}; largest: {:.4f} "
        "MB; peak storage: {:.4f} MB.\n\n",
        cell.n_terms, cell.n_distinct, cell.n_reused, mb(cell.largest_mem),
        mb(cell.peak_storage));
    out << std::format("Total FLOPs (symbolic): {}\n\n",
                       cell.total_flops.text());
    out << std::format("Total FLOPs (computed, index extents): {:.4e}\n\n",
                       cell.total_flops.ops());
    out << "### Largest intermediates\n\n";
    report_largest(cell, cfg.out.top_n, out);
    out << "### Most expensive contractions\n\n";
    report_expensive(cell, cfg.out.top_n, out);
    out << "### Shape census\n\n";
    report_shape_histogram(cell, out);
  }
  if (cfg.cache.enabled) {
    out << std::format(
        "## Cache (gated, volatile leaf = \"{}\", min_repeats = {})\n\n",
        cfg.optimize.volatile_leaf, cfg.cache.min_repeats);
    out << "| Metric | Value |\n|---|---|\n";
    out << std::format("| Cached | {} |\n", sim.n_cached);
    out << std::format("| Persistent | {} |\n", sim.n_persistent);
    out << std::format("| Persistent footprint (MB) | {:.4f} |\n",
                       mb(sim.persistent_footprint));
    out << std::format("| Total cached footprint (MB) | {:.4f} |\n\n",
                       mb(sim.cached_footprint));
  }
}
