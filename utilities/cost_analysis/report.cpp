#include "report.hpp"

#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/version.hpp>

#include <algorithm>
#include <format>
#include <map>

using namespace sequant;

namespace {
// Bytes per stored scalar (real=8, complex=16), threaded from the field.
double mb(const AsyCost& c, double bytes_per_elem) {
  return c.ops() * bytes_per_elem / (1024.0 * 1024.0);
}

void report_largest(const EquationResult& cell, std::size_t top_n,
                    double bytes_per_elem, std::ostream& out) {
  std::vector<std::pair<const TreeNode*, const Record*>> rows;
  for (const auto& [node, rec] : cell.catalog)
    if (rec.memory != AsyCost::zero()) rows.emplace_back(&node, &rec);
  // Sort by size (largest first); tie-break on the construction string so the
  // row order is deterministic regardless of unordered_map iteration order.
  std::ranges::sort(rows, [bytes_per_elem](const auto& a, const auto& b) {
    const double sa = mb(a.second->memory, bytes_per_elem),
                 sb = mb(b.second->memory, bytes_per_elem);
    if (sa != sb) return sb < sa;
    return full_expr(*a.first) < full_expr(*b.first);
  });
  out << "| Rank | Result | Spaces | Memory (MB) | Uses | Construction | "
         "Local ops |\n|---|---|---|---|---|---|---|\n";
  for (std::size_t i = 0; i < rows.size() && i < top_n; ++i) {
    const auto& [node, rec] = rows[i];
    out << std::format("| {} | {} | {} | {:.4f} | {} | {} | {} |\n", i + 1,
                       rec->label, rec->spaces, mb(rec->memory, bytes_per_elem),
                       rec->uses, full_expr(*node), rec->flops.text());
  }
}

void report_expensive(const EquationResult& cell, std::size_t top_n,
                      std::ostream& out) {
  std::vector<std::pair<const TreeNode*, const Record*>> rows;
  for (const auto& [node, rec] : cell.catalog)
    if (rec.flops != AsyCost::zero()) rows.emplace_back(&node, &rec);
  // Most FLOPs first, tie-broken on the construction string for determinism.
  std::ranges::sort(rows, [](const auto& a, const auto& b) {
    const double fa = a.second->flops.ops(), fb = b.second->flops.ops();
    if (fa != fb) return fb < fa;
    return full_expr(*a.first) < full_expr(*b.first);
  });
  out << "| Rank | Result | Spaces | Uses | Construction | "
         "Operations |\n|---|---|---|---|---|---|\n";
  for (std::size_t i = 0; i < rows.size() && i < top_n; ++i) {
    const auto& [node, rec] = rows[i];
    out << std::format("| {} | {} | {} | {} | {} | {} |\n", i + 1, rec->label,
                       rec->spaces, rec->uses, full_expr(*node),
                       rec->flops.text());
  }
}

void report_shape_histogram(const EquationResult& cell, double bytes_per_elem,
                            std::ostream& out) {
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
    // Representative size of the shape = its largest instance. 'X' is meant
    // for aux spaces (Κ, L), but space_signature() also puts any other
    // non-hole/non-particle space there (e.g. frozen-core), so a single O/V/X
    // signature can group instances of different sizes even with a single aux
    // space (or none) registered; taking the max keeps this deterministic, at
    // the cost of a coarser (not exact) census size.
    if (mb(g.mem, bytes_per_elem) < mb(rec.memory, bytes_per_elem))
      g.mem = rec.memory;
  }
  std::vector<std::pair<std::string, Agg>> rows(by_shape.begin(),
                                                by_shape.end());
  std::ranges::sort(rows, [bytes_per_elem](const auto& a, const auto& b) {
    const double sa = mb(a.second.mem, bytes_per_elem),
                 sb = mb(b.second.mem, bytes_per_elem);
    if (sa != sb) return sb < sa;
    return a.first < b.first;  // tie-break on shape string
  });
  out << "| Shape | Size (MB) | # Distinct | Total uses |\n|---|---|---|---|\n";
  for (const auto& [shape, g] : rows)
    out << std::format("| {} | {:.4f} | {} | {} |\n", shape,
                       mb(g.mem, bytes_per_elem), g.count, g.uses);
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
    const std::vector<std::pair<std::string, EquationResult>>& results,
    const SimResult& sim, std::ostream& out) {
  const double bytes_per_elem = cfg.real_field ? 8.0 : 16.0;
  // Sections are separated by a leading blank line rather than a trailing one,
  // so the report ends with a single newline (matches the committed fixtures).
  out << "# SeQuant cost analysis\n";
  if (!cfg.out.omit_revision)
    out << std::format("\n_SeQuant: {}_\n", git_revision());
  out << std::format("\nSizes in MB ({} bytes/element).\n", bytes_per_elem);
  for (const auto& [name, cell] : results) {
    out << std::format("\n## {}\n", name);
    out << std::format(
        "\nTerms: {}; distinct intermediates: {}; reused: {}; largest: {:.4f} "
        "MB; peak storage: {:.4f} MB.\n",
        cell.n_terms, cell.n_distinct, cell.n_reused,
        mb(cell.largest_mem, bytes_per_elem),
        mb(cell.peak_storage, bytes_per_elem));
    // Op-count only; no concrete FLOP number.
    out << std::format("\nTotal operations (symbolic): {}\n",
                       cell.total_flops.text());
    out << "\n### Largest intermediates\n\n";
    report_largest(cell, cfg.out.top_n, bytes_per_elem, out);
    out << "\n### Most expensive contractions\n\n";
    report_expensive(cell, cfg.out.top_n, out);
    out << "\n### Shape census\n\n";
    report_shape_histogram(cell, bytes_per_elem, out);
  }
  if (cfg.cache.enabled) {
    out << std::format(
        "\n## Cache (gated, volatile leaf = \"{}\", min_repeats = {})\n\n",
        cfg.optimize.volatile_leaf, cfg.cache.min_repeats);
    out << "| Metric | Value |\n|---|---|\n";
    out << std::format("| Cached | {} |\n", sim.n_cached);
    out << std::format("| Persistent | {} |\n", sim.n_persistent);
    out << std::format("| Persistent footprint (MB) | {:.4f} |\n",
                       mb(sim.persistent_footprint, bytes_per_elem));
    out << std::format("| Total cached footprint (MB) | {:.4f} |\n",
                       mb(sim.cached_footprint, bytes_per_elem));
  }
}
