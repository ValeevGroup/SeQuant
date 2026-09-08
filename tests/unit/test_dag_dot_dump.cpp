//
// Hidden dump utility ([.][dot-dump]): writes two Graphviz files for the
// serialized CSV-CCSD residual equations (energy, singles, doubles), optimized
// UNBATCHED and binarized one tree per equation exactly as the application
// does -- the FOREST (one node per tree node) and the value DAG (one node per
// distinct value, as compute_dag_boulevard folds the forest). In both, a value
// with more than one occurrence in the forest gets its own colour (the same
// colour in both pictures); a value that occurs once, every leaf, and every
// edge are black.
//
//   SEQUANT_DOT_FOREST=<path> SEQUANT_DOT_DAG=<path> [SEQUANT_DOT_LABELS=1]
//   [SEQUANT_DOT_COLLAPSE_SUMS=1] [SEQUANT_DOT_LATEX=<path>]
//   unit_tests-sequant "[dot-dump]"
//
// SEQUANT_DOT_LATEX writes the forest as LaTeX: one display equation per
// root, the head tensor equal to the optimized (factorized) expression, in
// breqn's dmath* environments so long right-hand sides break automatically.
//
// With SEQUANT_DOT_LABELS set, nodes carry their value id as a label;
// otherwise they are unlabeled dots.
//

#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::string dotdump_slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

}  // namespace

TEST_CASE(
    "dot dump: the residual forest and its value DAG, shared values coloured",
    "[.][dot-dump]") {
  using namespace sequant;
  using namespace sequant::eval::dryrun;
  namespace ev = sequant::eval;

  char const* const forest_path = std::getenv("SEQUANT_DOT_FOREST");
  char const* const dag_path = std::getenv("SEQUANT_DOT_DAG");
  REQUIRE(forest_path != nullptr);
  REQUIRE(dag_path != nullptr);
  bool const labels = std::getenv("SEQUANT_DOT_LABELS") != nullptr;
  char const* const latex_path = std::getenv("SEQUANT_DOT_LATEX");
  // SEQUANT_DOT_COLLAPSE_SUMS: a chain of Sum nodes (an equation's binarized
  // summation spine, every link a one-time node) is drawn as its top node
  // only, with the term trees attached to it directly.
  bool const collapse_sums =
      std::getenv("SEQUANT_DOT_COLLAPSE_SUMS") != nullptr;

  auto ctx0 = get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx0));

  // The reference size regime (the same numbers the dry-run fixtures use), so
  // the unbatched factorization is the one a real run of that size picks.
  SizeRegime regime;
  regime.space_extent = {{L"i", 80u}, {L"μ̃", 896u}, {L"Κ", 1682u}};
  regime.csv_pno_moment = {1.0, 23.175775480059084, 25.865548281212597,
                           28.171416142614103, 30.03848680550367};
  regime.csv_osv_moment = {1.0, 58.987499999999997, 59.289227520688783,
                           59.584437469011633, 59.872014818179686};
  auto idx_to_extent = [&regime](Index const& idx) -> std::size_t {
    if (idx.has_proto_indices())
      return static_cast<std::size_t>(
          std::ceil(std::max(2.0, regime.inner_pow(idx, 1))));
    return regime.extent(idx);
  };

  OptimizeOptions opts;  // no batch policy: unbatched
  opts.objective_function = ObjectiveFunction::DenseTimeSpace;
  opts.reorder = ReorderSum::Reorder;
  opts.idx_to_extent = idx_to_extent;
  opts.inner_pow = regime.inner_pow_fn();
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;

  auto make_head = [](std::size_t rank) -> Tensor {
    std::vector<Index> occ, vir;
    for (std::size_t k = 1; k <= rank; ++k)
      occ.emplace_back(L"i_" + std::to_wstring(k));
    for (std::size_t k = 1; k <= rank; ++k)
      vir.emplace_back(L"a_" + std::to_wstring(k), occ);
    return Tensor(rank == 0 ? L"E" : L"R", bra(vir), ket(occ),
                  Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                  ColumnSymmetry::Symm);
  };

  struct Eqn {
    char const* file;
    std::size_t rank;
  };
  std::vector<EvalNodeDryRun> forest;
  for (Eqn const& eq : {Eqn{"/data/csv_ccsd_energy_df.txt", 0u},
                        Eqn{"/data/csv_ccsd_singles_residual_df.txt", 1u},
                        Eqn{"/data/csv_ccsd_doubles_residual_df.txt", 2u}}) {
    auto const body =
        dotdump_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) + eq.file);
    REQUIRE(!body.empty());
    std::string line = body;
    if (auto nl = line.find('\n'); nl != std::string::npos)
      line = line.substr(0, nl);
    auto sum = deserialize<ExprPtr>(line);
    REQUIRE(sum);
    auto res = optimize_result(sum, opts);
    REQUIRE(res.expr);
    ResultExpr rexpr{make_head(eq.rank), res.expr};
    if (latex_path) {
      static std::ofstream tex;
      if (!tex.is_open()) {
        tex.open(latex_path);
        REQUIRE(tex.good());
        tex << "% forest of the optimized (factorized) equations; needs "
               "\\usepackage{breqn}\n";
      }
      auto const w2s = [](std::wstring const& w) {
        std::string out;
        for (wchar_t c : w) {
          if (c < 0x80) {
            out += static_cast<char>(c);
          } else {  // UTF-8 encode
            unsigned int const u = static_cast<unsigned int>(c);
            if (u < 0x800) {
              out += static_cast<char>(0xC0 | (u >> 6));
              out += static_cast<char>(0x80 | (u & 0x3F));
            } else if (u < 0x10000) {
              out += static_cast<char>(0xE0 | (u >> 12));
              out += static_cast<char>(0x80 | ((u >> 6) & 0x3F));
              out += static_cast<char>(0x80 | (u & 0x3F));
            } else {
              out += static_cast<char>(0xF0 | (u >> 18));
              out += static_cast<char>(0x80 | ((u >> 12) & 0x3F));
              out += static_cast<char>(0x80 | ((u >> 6) & 0x3F));
              out += static_cast<char>(0x80 | (u & 0x3F));
            }
          }
        }
        return out;
      };
      tex << "\\begin{dmath*}\n"
          << w2s(make_head(eq.rank).to_latex()) << " = "
          << w2s(res.expr->to_latex()) << "\n\\end{dmath*}\n\n";
    }
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    forest.push_back(binarize<EvalExprDryRun>(rexpr, BinarizationOptions{}));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }
  REQUIRE(forest.size() == 3);

  // The value DAG: one cell per distinct value, occurrences = every tree node
  // that is that value.
  CostModel const cm{regime};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };
  auto rich = ev::compute_dag_boulevard(forest, cm, block_of);
  REQUIRE(!rich.cells.empty());
  auto const g = ev::detail::ordered_schedule_dep_graph(rich);

  // Colour per shared (multi-occurrence, non-leaf) value: hues spread evenly.
  std::unordered_map<std::size_t, std::string> colour_of_hash;  // hash -> HSV
  {
    std::vector<std::size_t> shared;
    for (auto const& c : rich.cells)
      if (!c.is_leaf && c.occurrences.size() > 1) shared.push_back(c.value_id);
    std::size_t const n = shared.size();
    for (std::size_t k = 0; k < n; ++k) {
      // interleave hues so neighbours in id order are far apart in hue
      double const h =
          std::fmod(0.618033988749895 * static_cast<double>(k), 1.0);
      double const s = 0.75 + 0.25 * static_cast<double>(k % 2);
      double const v = 0.70 + 0.25 * static_cast<double>((k / 2) % 2);
      std::ostringstream os;
      os << std::fixed << std::setprecision(4) << h << " " << s << " " << v;
      colour_of_hash[rich.cells[shared[k]].hash] = os.str();
    }
    std::cerr << "[dot-dump] values=" << rich.cells.size() << " shared=" << n
              << "\n";
  }
  std::unordered_map<std::size_t, std::size_t> vid_of_hash;
  for (auto const& c : rich.cells) vid_of_hash.emplace(c.hash, c.value_id);

  // Roots (the forest's trees, in order): drawn as labelled double circles on
  // one rank in both pictures.
  std::vector<std::pair<std::size_t, std::string>> root_labels;  // hash, label
  for (std::size_t t = 0; t < forest.size(); ++t)
    root_labels.push_back(
        {forest[t]->hash_value(),
         t == 0 ? std::string("E") : "R" + std::to_string(t)});
  auto root_label_of = [&](std::size_t hash) -> std::string {
    for (auto const& [h, l] : root_labels)
      if (h == hash) return l;
    return {};
  };
  auto node_attrs = [&](std::size_t hash, bool leaf) -> std::string {
    std::string attrs;
    if (auto const rl = root_label_of(hash); !rl.empty())
      return "shape=doublecircle, width=0.34, height=0.34, fixedsize=true, "
             "style=filled, fillcolor=white, color=black, label=\"" +
             rl + "\", fontsize=11, fontname=\"Helvetica-Bold\"";
    auto const it = colour_of_hash.find(hash);
    if (!leaf && it != colour_of_hash.end())
      attrs = "style=filled, fillcolor=\"" + it->second + "\", color=\"" +
              it->second + "\"";
    else
      attrs = std::string("style=filled, fillcolor=black, color=black") +
              (leaf ? ", shape=box, width=0.12, height=0.12" : "");
    if (labels) {
      auto const vit = vid_of_hash.find(hash);
      attrs += ", label=\"" +
               (vit != vid_of_hash.end() ? std::to_string(vit->second)
                                         : std::string("?")) +
               "\", fontsize=7, fontcolor=white";
    }
    return attrs;
  };

  // --- forest ---
  {
    std::ofstream out(forest_path);
    REQUIRE(out.good());
    out << "digraph forest {\n  rankdir=TB; nodesep=0.06; ranksep=0.25;\n"
        << "  node [shape=circle, width=0.14, height=0.14, fixedsize=true, "
           "label=\"\"];\n  edge [color=black, arrowhead=none, "
           "penwidth=0.5];\n";
    std::size_t counter = 0;
    auto is_sum = [](EvalNodeDryRun const& n) {
      return !n.leaf() && n->op_type() == EvalOp::Sum;
    };
    // emit(n, attach): draws n's subtree; when n is a Sum inside a Sum chain
    // and collapsing is on, n itself is skipped and its children attach to
    // `attach` (the chain's top node) instead.
    // emit_under(child, parent, parent_is_sum): draws child's subtree under
    // parent; a Sum child of a Sum parent is a link of a summation spine and,
    // with collapsing on, is skipped -- its own children attach to parent.
    std::function<void(EvalNodeDryRun const&, std::size_t, bool)> emit_under;
    std::function<std::size_t(EvalNodeDryRun const&)> emit =
        [&](EvalNodeDryRun const& n) -> std::size_t {
      std::size_t const id = counter++;
      out << "  n" << id << " [" << node_attrs(n->hash_value(), n.leaf())
          << "];\n";
      if (!n.leaf()) {
        emit_under(n.left(), id, is_sum(n));
        emit_under(n.right(), id, is_sum(n));
      }
      return id;
    };
    emit_under = [&](EvalNodeDryRun const& child, std::size_t parent,
                     bool parent_is_sum) {
      if (collapse_sums && is_sum(child) && parent_is_sum) {
        emit_under(child.left(), parent, true);
        emit_under(child.right(), parent, true);
        return;
      }
      auto const cid = emit(child);
      out << "  n" << parent << " -> n" << cid << ";\n";
    };
    std::vector<std::size_t> root_ids;
    for (auto const& tree : forest) root_ids.push_back(emit(tree));
    out << "  {rank=same;";
    for (std::size_t r : root_ids) out << " n" << r;
    out << "}\n}\n";
    std::cerr << "[dot-dump] forest nodes=" << counter << " -> " << forest_path
              << "\n";
  }

  // --- value DAG ---
  {
    std::ofstream out(dag_path);
    REQUIRE(out.good());
    out << "digraph dag {\n  rankdir=TB; nodesep=0.06; ranksep=0.25;\n"
        << "  node [shape=circle, width=0.14, height=0.14, fixedsize=true, "
           "label=\"\"];\n  edge [color=black, arrowhead=none, "
           "penwidth=0.5];\n";
    auto const vmap = ev::build_value_node_map(forest);
    auto cell_is_sum = [&](std::size_t vid) {
      auto const it = vmap.find(rich.cells[vid].hash);
      return it != vmap.end() && !it->second.leaf() &&
             it->second->op_type() == EvalOp::Sum;
    };
    // a Sum cell whose every consumer is a Sum is a chain link: skipped, its
    // operands re-attached to the nearest non-skipped ancestor
    std::vector<bool> skipped(rich.cells.size(), false);
    if (collapse_sums)
      for (auto const& c : rich.cells) {
        if (!cell_is_sum(c.value_id)) continue;
        auto const cit = g.consumers_of.find(c.value_id);
        if (cit == g.consumers_of.end() || cit->second.empty()) continue;
        bool all_sum = true;
        for (std::size_t u : cit->second) all_sum = all_sum && cell_is_sum(u);
        skipped[c.value_id] = all_sum;
      }
    std::function<std::size_t(std::size_t)> attach_point = [&](std::size_t v) {
      if (!skipped[v]) return v;
      return attach_point(g.consumers_of.at(v).front());
    };
    for (auto const& c : rich.cells)
      if (!skipped[c.value_id])
        out << "  v" << c.value_id << " [" << node_attrs(c.hash, c.is_leaf)
            << "];\n";
    std::size_t edges = 0;
    for (auto const& [parent, ops] : g.depends_on)
      for (std::size_t o : ops) {
        if (skipped[o]) continue;  // its operands attach higher up
        out << "  v" << attach_point(parent) << " -> v" << o << ";\n";
        ++edges;
      }
    out << "  {rank=same;";
    for (auto const& [h, l] : root_labels)
      if (auto const it = vid_of_hash.find(h); it != vid_of_hash.end())
        out << " v" << it->second;
    out << "}\n}\n";
    std::cerr << "[dot-dump] dag nodes=" << rich.cells.size()
              << " edges=" << edges << " -> " << dag_path << "\n";
  }
}
