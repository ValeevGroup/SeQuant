#ifndef SEQUANT_UTILITIES_COST_ANALYSIS_COST_ANALYSIS_HPP
#define SEQUANT_UTILITIES_COST_ANALYSIS_COST_ANALYSIS_HPP

/// \file
/// Data types shared between the pipeline (`cost_analysis.cpp`) and the
/// reporter (`report.cpp`): the parsed JSON `Config` and the per-equation
/// analysis results. Logic lives in the `.cpp` files.

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>

#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

/// `optimize` driver block; fields map onto `sequant::OptimizeOptions` (see
/// there for their semantics). `volatile_leaf` is the amplitude label solved
/// for (empty => no volatile weighting).
struct OptSpec {
  std::string objective = "dense_flops";
  bool reorder = true, cse_subnet = false;
  std::string volatile_leaf = "R";
  double volatile_weight = 10.0, machine_balance = 0.0, fast_mem_elems = 0.0;
};

/// `cache` driver block: parameters of the gated cache-reuse simulation.
struct CacheSpec {
  bool enabled = true;
  std::size_t min_repeats = 2;
  double max_footprint = 0.0;
};

/// `output` driver block: report path and verbosity.
struct OutSpec {
  std::string path = "cost_analysis.md";
  std::size_t top_n = 20;  ///< rows in the largest/expensive tables
  bool dump_tree = false;  ///< also write each result's binarized tree
};

/// One `results[]` entry: a named equation file to analyze.
struct ResultSpec {
  std::string name, equation_file;
};

/// The parsed JSON driver.
struct Config {
  bool spinor = true;
  bool real_field = true;  // true -> Field::Real (the utility's default); the
                           // field is a single global choice for a computation
  std::string convention = "sr";  // standard registry: min_sr|sr|mr|f12
  std::vector<std::string>
      aux;  // extra factorization spaces to register so
            // DF/THC equations parse: df (Κ) and/or thc (L)
  std::map<std::string, std::size_t> sizes;  // space label -> approximate size,
                                             // overriding the standard registry
  OptSpec optimize;
  CacheSpec cache;
  OutSpec out;
  std::vector<ResultSpec> results;
};

/// A binarized evaluation-tree node, keyed in the catalog by structural
/// hash/equality so equal intermediates across terms collapse to one entry.
using TreeNode = sequant::EvalNode<sequant::EvalExpr>;
using Hasher = sequant::TreeNodeHasher<TreeNode>;
using Comparator = sequant::TreeNodeEqualityComparator<TreeNode>;

/// Catalog entry for one distinct intermediate (eval-tree node).
struct Record {
  std::size_t uses = 0;            ///< how many terms build this node
  sequant::AsyCost memory, flops;  ///< result footprint and local op count
  std::string label, spaces;       ///< tensor label and O/V/X signature
};

/// Cost analysis of one result equation.
struct CellResult {
  std::size_t n_terms = 0, n_distinct = 0, n_reused = 0;
  sequant::AsyCost largest_mem, peak_storage, total_flops;
  std::unordered_map<TreeNode, Record, Hasher, Comparator> catalog;
  /// Binarized tree per summand, built once and reused by the catalog scan,
  /// the cache simulation, and --dump_tree.
  std::vector<TreeNode> trees;
};

/// Cache-simulation totals across all analyzed results.
struct SimResult {
  std::size_t n_cached = 0, n_persistent = 0;
  sequant::AsyCost cached_footprint, persistent_footprint;
};

#endif  // SEQUANT_UTILITIES_COST_ANALYSIS_COST_ANALYSIS_HPP
