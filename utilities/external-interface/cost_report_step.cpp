#include "cost_report_step.hpp"
#include "filter_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"
#include "utils.hpp"

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/version.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstddef>
#include <cwctype>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <ostream>
#include <ranges>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(CostReportStep);

namespace {

/// A binarized evaluation-tree node, keyed in the catalog by structural
/// hash/equality so equal intermediates across terms collapse to one entry.
using TreeNode = EvalNode<EvalExpr>;
using Hasher = TreeNodeHasher<TreeNode>;
using Comparator = TreeNodeEqualityComparator<TreeNode>;

/// Catalog entry for one distinct intermediate (eval-tree node).
struct Record {
  std::size_t uses = 0;       ///< how often the trees build this node
  AsyCost memory, flops;      ///< result footprint and local op count
  std::string label, spaces;  ///< tensor label and space signature
};

/// Cost analysis of one result expression.
struct ExpressionCost {
  std::string name;
  std::size_t n_terms = 0, n_distinct = 0, n_reused = 0;
  AsyCost largest_mem, peak_storage, total_flops;
  std::unordered_map<TreeNode, Record, Hasher, Comparator> catalog;
  /// Binarized tree per summand
  std::vector<TreeNode> trees;
};

/// Cache-simulation totals across all analyzed expressions.
struct CacheCost {
  std::size_t n_cached = 0, n_persistent = 0;
  AsyCost cached_footprint, persistent_footprint;
};

std::string space_signature(const EvalExpr &ev,
                            const CostReportStep::Options &options) {
  if (!ev.is_tensor()) return "scalar";
  std::string result;
  for (const auto &idx : ev.as_tensor().const_braketaux_indices()) {
    auto it = options.space_symbols.find(
        mbpt::spinannotation_remove(idx.space().base_key()));
    result += it != options.space_symbols.end() ? std::string_view(it->second)
                                                : std::string_view("X");
  }
  return result;
}

AsyCost result_memory(const EvalExpr &ev) {
  if (!ev.is_tensor()) return AsyCost::zero();
  AsyCost::ExponentMap exponents;
  // Aux indices count towards the size, consistent with space_signature
  for (const auto &idx : ev.as_tensor().const_braketaux_indices())
    ++exponents[idx.basis()];
  return AsyCost{std::move(exponents)};
}

// Like min_storage() in eval_node.hpp, but maxes numerically (by .ops(), the
// size at the registered index-space sizes) rather than by
// AsyCost::operator<, which orders by symbolic polynomial degree: at actual
// sizes a lower-degree node can be the larger one.
AsyCost numeric_peak_storage(const TreeNode &tree) {
  AsyCost result = AsyCost::zero();
  tree.visit([&](const TreeNode &n) {
    AsyCost cost = AsyCost::zero();
    if (n.leaf()) {
      if (n->is_tensor()) cost = result_memory(*n);
    } else {
      if (n.left()->is_tensor()) cost = cost + result_memory(*n.left());
      if (n.right()->is_tensor()) cost = cost + result_memory(*n.right());
      if (n->is_tensor()) cost = cost + result_memory(*n);
    }
    if (result.ops() < cost.ops()) result = cost;
  });
  return result;
}

/// Calls @p f on every term of @p expr. Summands of nested sums are terms too;
/// products are left as they are, since their nesting is the factorization
/// chosen by the optimizer.
template <typename F>
void for_each_term(const ExprPtr &expr, F &&f) {
  if (expr->is<Sum>()) {
    for (const ExprPtr &summand : expr->as<Sum>().summands())
      for_each_term(summand, f);
  } else {
    f(expr);
  }
}

ExpressionCost analyze(std::string name, const ResultExpr &expr,
                       const CostReportStep::Options &options) {
  ExpressionCost res;
  res.name = std::move(name);

  ExprPtr rhs = expr.expression()->clone();
  pop_symmetrizer(rhs);

  auto scan = [&](const TreeNode &tree) {
    tree.visit_internal([&](const TreeNode &node) {
      // A term's root is its result, not an intermediate
      if (node.root()) return;
      auto it = res.catalog.find(node);
      if (it == res.catalog.end()) {
        Record rec;
        rec.uses = 1;
        rec.memory = result_memory(*node);
        rec.flops = Flops{}(node);
        rec.label = node->label();
        rec.spaces = space_signature(*node, options);
        res.catalog.emplace(node, std::move(rec));
      } else {
        ++it->second.uses;
      }
    });
  };

  auto handle = [&](const ExprPtr &summand) {
    ResultExpr term = expr;
    term.expression() = summand;
    auto tree = binarize<EvalExpr>(term);
    scan(tree);
    // Peak storage is the heaviest single operation over all terms; total
    // operations are summed over them
    const AsyCost ps = numeric_peak_storage(tree);
    if (res.peak_storage.ops() < ps.ops()) res.peak_storage = ps;
    res.total_flops = res.total_flops + asy_cost(tree);
    res.trees.push_back(std::move(tree));
  };

  for_each_term(rhs, handle);
  res.n_terms = res.trees.size();

  for (const auto &[_, rec] : res.catalog) {
    if (rec.memory == AsyCost::zero()) continue;
    ++res.n_distinct;
    if (res.largest_mem.ops() < rec.memory.ops()) res.largest_mem = rec.memory;
    if (rec.uses >= 2) ++res.n_reused;
  }
  return res;
}

CacheCost simulate_cache(const CostReportStep::CacheOptions &options,
                         const std::vector<ExpressionCost> &costs) {
  // A view rather than a copy: copying a tree copies all its nodes
  auto forest =
      costs | std::views::transform(&ExpressionCost::trees) | std::views::join;

  auto is_volatile = [&](const TreeNode &n) {
    return n.leaf() && n->is_tensor() &&
           options.is_volatile->matches(n->as_tensor());
  };
  auto cm = cache_manager(
      forest, is_volatile, options.min_repeats,
      [](const TreeNode &n) { return result_memory(*n).ops(); },
      options.max_footprint);

  CacheCost r;
  cm.for_each_key([&](const TreeNode &n) {
    const AsyCost mem = result_memory(*n);
    r.cached_footprint += mem;
    ++r.n_cached;
    if (cm.persistent(n)) {
      r.persistent_footprint += mem;
      ++r.n_persistent;
    }
  });
  return r;
}

/// Serializes a binarized tree back to a parenthesized string
std::string full_expr(const TreeNode &n) {
  if (n.leaf()) return n->label();
  if (n->op_type() == EvalOp::Adjoint) return full_expr(n.left()) + "+";
  if (n->op_type() == EvalOp::Broadcast) return full_expr(n.left());
  if (n->op_type() == EvalOp::Power)
    return full_expr(n.left()) + "^(" + sequant::to_string(power_exponent(n)) +
           ")";
  const char *op = n->op_type() == EvalOp::Product ? " * " : " + ";
  return "(" + full_expr(n.left()) + op + full_expr(n.right()) + ")";
}

double mb(const AsyCost &c, double bytes_per_elem) {
  return c.ops() * bytes_per_elem / (1024.0 * 1024.0);
}

/// Bytes per stored element: 16 if any index of the analyzed expressions lives
/// in a complex-valued space, 8 otherwise
double bytes_per_element(const std::vector<ExpressionCost> &costs) {
  for (const auto &cost : costs)
    for (const auto &tree : cost.trees) {
      bool complex = false;
      tree.visit([&](const TreeNode &n) {
        if (complex || !n->is_tensor()) return;
        for (const auto &idx : n->as_tensor().const_braketaux_indices())
          if (idx.basis().field() == Field::Complex) complex = true;
      });
      if (complex) return 16.0;
    }
  return 8.0;
}

void report_largest(const ExpressionCost &cost, std::size_t top_n,
                    double bytes_per_elem, std::ostream &out) {
  std::vector<std::pair<const TreeNode *, const Record *>> rows;
  for (const auto &[node, rec] : cost.catalog)
    if (rec.memory != AsyCost::zero()) rows.emplace_back(&node, &rec);
  // Largest first; ties broken on the construction string so that the order
  // does not depend on unordered_map iteration
  std::ranges::sort(rows, [bytes_per_elem](const auto &a, const auto &b) {
    const double sa = mb(a.second->memory, bytes_per_elem),
                 sb = mb(b.second->memory, bytes_per_elem);
    if (sa != sb) return sb < sa;
    return full_expr(*a.first) < full_expr(*b.first);
  });
  out << "| Rank | Result | Spaces | Memory (MB) | Uses | Construction | "
         "Local ops |\n|---|---|---|---|---|---|---|\n";
  for (std::size_t i = 0; i < rows.size() && i < top_n; ++i) {
    const auto &[node, rec] = rows[i];
    out << std::format("| {} | {} | {} | {:.4f} | {} | {} | {} |\n", i + 1,
                       rec->label, rec->spaces, mb(rec->memory, bytes_per_elem),
                       rec->uses, full_expr(*node), rec->flops.text());
  }
}

void report_expensive(const ExpressionCost &cost, std::size_t top_n,
                      std::ostream &out) {
  std::vector<std::pair<const TreeNode *, const Record *>> rows;
  for (const auto &[node, rec] : cost.catalog)
    if (rec.flops != AsyCost::zero()) rows.emplace_back(&node, &rec);
  std::ranges::sort(rows, [](const auto &a, const auto &b) {
    const double fa = a.second->flops.ops(), fb = b.second->flops.ops();
    if (fa != fb) return fb < fa;
    return full_expr(*a.first) < full_expr(*b.first);
  });
  out << "| Rank | Result | Spaces | Uses | Construction | "
         "Operations |\n|---|---|---|---|---|---|\n";
  for (std::size_t i = 0; i < rows.size() && i < top_n; ++i) {
    const auto &[node, rec] = rows[i];
    out << std::format("| {} | {} | {} | {} | {} | {} |\n", i + 1, rec->label,
                       rec->spaces, rec->uses, full_expr(*node),
                       rec->flops.text());
  }
}

void report_shape_census(const ExpressionCost &cost, double bytes_per_elem,
                         std::ostream &out) {
  struct Agg {
    std::size_t count = 0, uses = 0;
    AsyCost mem;
  };
  std::map<std::string, Agg> by_shape;
  for (const auto &[_, rec] : cost.catalog) {
    if (rec.memory == AsyCost::zero()) continue;
    auto &g = by_shape[rec.spaces];
    ++g.count;
    g.uses += rec.uses;
    // One signature can cover spaces of different sizes (every unlisted space
    // is 'X'), so the shape's size is that of its largest instance
    if (g.mem.ops() < rec.memory.ops()) g.mem = rec.memory;
  }
  std::vector<std::pair<std::string, Agg>> rows(by_shape.begin(),
                                                by_shape.end());
  std::ranges::sort(rows, [](const auto &a, const auto &b) {
    const double sa = a.second.mem.ops(), sb = b.second.mem.ops();
    if (sa != sb) return sb < sa;
    return a.first < b.first;
  });
  out << "| Shape | Size (MB) | # Distinct | Total uses |\n|---|---|---|---|\n";
  for (const auto &[shape, g] : rows)
    out << std::format("| {} | {:.4f} | {} | {} |\n", shape,
                       mb(g.mem, bytes_per_elem), g.count, g.uses);
}

void write_report(const CostReportStep::Options &options,
                  const std::vector<ExpressionCost> &costs,
                  const CacheCost &cache, std::ostream &out) {
  const double bytes_per_elem = bytes_per_element(costs);
  out << "# SeQuant cost analysis\n";
  if (!options.omit_revision)
    out << std::format("\n_SeQuant: {}_\n", git_revision());
  out << std::format("\nSizes in MB ({} bytes/element).\n", bytes_per_elem);
  for (const auto &cost : costs) {
    out << std::format("\n## {}\n", cost.name);
    out << std::format(
        "\nTerms: {}; distinct intermediates: {}; reused: {}; largest: {:.4f} "
        "MB; peak storage: {:.4f} MB.\n",
        cost.n_terms, cost.n_distinct, cost.n_reused,
        mb(cost.largest_mem, bytes_per_elem),
        mb(cost.peak_storage, bytes_per_elem));
    out << std::format("\nTotal operations (symbolic): {}\n",
                       cost.total_flops.text());
    out << "\n### Largest intermediates\n\n";
    report_largest(cost, options.top_n, bytes_per_elem, out);
    out << "\n### Most expensive contractions\n\n";
    report_expensive(cost, options.top_n, out);
    out << "\n### Shape census\n\n";
    report_shape_census(cost, bytes_per_elem, out);
  }
  if (options.cache) {
    out << std::format(
        "\n## Cache (gated, volatile leaf = {}, min_repeats = {})\n\n",
        options.cache->volatile_spec, options.cache->min_repeats);
    out << "| Metric | Value |\n|---|---|\n";
    out << std::format("| Cached | {} |\n", cache.n_cached);
    out << std::format("| Persistent | {} |\n", cache.n_persistent);
    out << std::format("| Persistent footprint (MB) | {:.4f} |\n",
                       mb(cache.persistent_footprint, bytes_per_elem));
    out << std::format("| Total cached footprint (MB) | {:.4f} |\n",
                       mb(cache.cached_footprint, bytes_per_elem));
  }
}

/// Section name of an expression without a result label
constexpr std::string_view unnamed = "Unnamed";

std::string expression_name(const ResultExpr &expr) {
  return expr.has_label() ? toUtf8(expr.label()) : std::string(unnamed);
}

/// A section name; `labeled` is false for the placeholder of an unlabeled
/// expression
struct SectionName {
  std::string_view name;
  bool labeled;
};

/// Throws if two of @p names coincide. With @p dump_tree the names are also
/// file names, so they must not contain a path separator and must differ in
/// more than case (as file systems may ignore it).
void check_names(const std::vector<SectionName> &names, bool dump_tree,
                 std::string_view kind) {
  // Folded name -> whether its first owner is labeled
  std::map<std::wstring, bool> keys;
  for (const auto &[name, labeled] : names) {
    std::wstring key = toUtf16(name);
    if (dump_tree) {
      if (key.find_first_of(L"/\\") != std::wstring::npos) {
        throw Exception(std::string(kind) + " cannot name a tree dump after '" +
                        std::string(name) +
                        "', which contains a path separator");
      }
      for (wchar_t &c : key) c = std::towlower(c);
    }
    const auto [it, inserted] = keys.emplace(std::move(key), labeled);
    if (inserted) continue;
    if (!labeled && !it->second) {
      throw Exception(std::string(kind) +
                      " got several expressions without a result label; "
                      "report them in separate steps");
    }
    throw Exception(std::string(kind) +
                    " got several expressions with result '" +
                    std::string(name) +
                    (dump_tree ? "' (ignoring case, as tree dumps are named "
                                 "after it)"
                               : "'") +
                    "; report them in separate steps");
  }
}

bool get_bool(const std::string &name, const nlohmann::json &value) {
  if (!value.is_boolean()) {
    throw Exception("Value for " + name + " must be a boolean");
  }
  return value.get<bool>();
}

}  // namespace

std::string CostReportStep::kind() const { return "cost_report"; }

bool CostReportStep::accepts_options() const { return true; }

bool CostReportStep::requires_options() const { return false; }

void CostReportStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    const std::string name = kind() + " option '" + key + "'";
    if (key == "output") {
      if (!value.is_string() || value.get<std::string>().empty()) {
        throw Exception("Value for " + name + " must be a non-empty string");
      }
      options_.output = value.get<std::string>();
    } else if (key == "top_n") {
      options_.top_n = parse_count(value, "Value for " + name);
    } else if (key == "dump_tree") {
      options_.dump_tree = get_bool(name, value);
    } else if (key == "omit_revision") {
      options_.omit_revision = get_bool(name, value);
    } else if (key == "space_symbols") {
      if (!value.is_object()) {
        throw Exception("Value for " + name +
                        " must be an object mapping space labels to "
                        "characters");
      }
      options_.space_symbols.clear();
      for (const auto &[label, symbol] : value.items()) {
        // One code point, which may take several bytes in UTF-8: exactly one
        // byte that is not a continuation byte (10xxxxxx)
        if (!symbol.is_string() ||
            std::ranges::count_if(symbol.get<std::string>(), [](char c) {
              return (static_cast<unsigned char>(c) & 0xC0) != 0x80;
            }) != 1) {
          throw Exception("Value for " + name + " entry '" + label +
                          "' must be a single character");
        }
        options_.space_symbols[toUtf16(label)] = symbol.get<std::string>();
      }
    } else if (key == "cache") {
      if (!value.is_object()) {
        throw Exception(name + " expects a JSON object");
      }
      CacheOptions cache;
      for (const auto &[ckey, cvalue] : value.items()) {
        const std::string cname = kind() + ".cache option '" + ckey + "'";
        if (ckey == "volatile") {
          cache.is_volatile =
              parse_tensor_filter(cvalue, kind() + ".cache.volatile");
          cache.volatile_spec = cvalue.dump();
        } else if (ckey == "min_repeats") {
          cache.min_repeats = parse_count(cvalue, "Value for " + cname);
        } else if (ckey == "max_footprint") {
          if (!cvalue.is_number() || cvalue.get<double>() < 0.0) {
            throw Exception("Value for " + cname +
                            " must be a non-negative number");
          }
          cache.max_footprint = cvalue.get<double>();
        } else {
          throw Exception("Unknown option key for " + kind() + ".cache: '" +
                          ckey + "'");
        }
      }
      if (!cache.is_volatile) {
        throw Exception("The 'volatile' option for " + kind() +
                        ".cache is mandatory");
      }
      options_.cache = std::move(cache);
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

void CostReportStep::check_space_symbols() const {
  const auto known = registered_space_labels();
  for (const auto &[label, _] : options_.space_symbols) {
    if (!known.contains(label)) {
      throw Exception("Unknown index space '" + toUtf8(label) + "' in " +
                      kind() + ".space_symbols");
    }
    // Shape signatures look up spaces with their spin annotation removed, so
    // an annotated key would never match
    if (mbpt::spinannotation_remove(label) != label) {
      throw Exception("Index space '" + toUtf8(label) + "' in " + kind() +
                      ".space_symbols carries a spin annotation; use its "
                      "spin-free label");
    }
  }
}

std::size_t CostReportStep::run(std::string_view step_id, ExecutionContext &ctx,
                                const std::vector<std::string_view> &inputs) {
  check_space_symbols();

  struct NamedExpr {
    const ResultExpr *expr;
    std::string name;
    std::string input_id;
  };
  // Gathered before any analysis, so that naming errors surface at once
  std::vector<NamedExpr> exprs;
  // Data reachable through several inputs is analyzed once
  std::set<std::size_t> seen;
  for (std::string_view input : inputs) {
    for (const ExecutionContext::Data<ProcessingData> &current :
         ctx.get_data(input)) {
      SEQUANT_ASSERT(!current.associated_ids.empty());
      if (!seen.insert(current.index).second) {
        continue;
      }
      const std::string input_id(current.associated_ids.front());
      try {
        const auto &data = convert_data<ExpressionData>(current.data.get());
        for (const ResultExpr &expr : data.expressions) {
          exprs.push_back({&expr, expression_name(expr), input_id});
        }
      } catch (const std::exception &e) {
        throw Exception("Error in " + kind() + " on input " + input_id + ": " +
                        e.what());
      }
    }
  }
  {
    std::vector<SectionName> names;
    names.reserve(exprs.size());
    for (const NamedExpr &e : exprs)
      names.push_back({e.name, e.expr->has_label()});
    check_names(names, options_.dump_tree, kind());
  }

  std::vector<ExpressionCost> costs;
  costs.reserve(exprs.size());
  for (const NamedExpr &e : exprs) {
    try {
      costs.push_back(analyze(e.name, *e.expr, options_));
    } catch (const std::exception &ex) {
      throw Exception("Error in " + kind() + " on input " + e.input_id + ": " +
                      ex.what());
    }
  }

  const std::filesystem::path output(
      options_.output.empty() ? std::string(step_id) + ".md" : options_.output);
  create_parent_directories(output, kind());

  // The report comes first, so that a run that fails writing it leaves no
  // tree dumps behind
  const CacheCost cache =
      options_.cache ? simulate_cache(*options_.cache, costs) : CacheCost{};

  std::ofstream out(output);
  if (!out)
    throw Exception("Cannot open output file '" + output.string() + "'");
  write_report(options_, costs, cache, out);
  out.close();
  if (!out) {
    throw Exception("Failed writing output file '" + output.string() + "'");
  }

  if (options_.dump_tree) {
    for (const ExpressionCost &cost : costs) {
      const std::filesystem::path path =
          output.parent_path() /
          (output.stem().string() + "." + cost.name + ".tree.txt");
      std::ofstream to(path);
      if (!to) throw Exception("Cannot open dump file '" + path.string() + "'");
      for (const auto &tree : cost.trees) to << full_expr(tree) << "\n";
      to.close();
      if (!to) {
        throw Exception("Failed writing dump file '" + path.string() + "'");
      }
    }
  }

  return 0;
}

}  // namespace sequant::util::extint
