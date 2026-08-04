// SeQuant expression cost/factorization analysis.
//
// Reads a JSON driver (space sizes + equation files), reads each equation in
// as-is, then optimizes -> binarizes it and writes a Markdown report on the
// largest intermediates, shape census, and (gated) cache footprint. Shared
// data types live in cost_analysis.hpp; report formatting lives in
// report.{hpp,cpp}.

#include "cost_analysis.hpp"
#include "report.hpp"

#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/io/shorthands.hpp>  // deserialize
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/runtime.hpp>  // set_locale
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/string.hpp>     // toUtf16
#include <SeQuant/domain/mbpt/convention.hpp>  // make_min_sr_spaces

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using nlohmann::json;
using namespace sequant;
using namespace sequant::mbpt;

namespace {

Config load_config(const json& d) {
  Config c;

  if (d.contains("context")) {
    const auto& ctx = d.at("context");
    const auto spbasis = ctx.value("spbasis", std::string("spinor"));
    if (spbasis != "spinor" && spbasis != "spinfree")
      throw Exception("unknown spbasis (want spinor|spinfree): " + spbasis);
    c.spinor = spbasis == "spinor";
    const auto field = ctx.value("field", std::string("real"));
    if (field != "real" && field != "complex")
      throw std::runtime_error("unknown field (want real|complex): " + field);
    c.real_field = field == "real";
    c.convention = ctx.value("convention", std::string("sr"));
    for (const auto& a : ctx.value("aux", json::array())) {
      const auto s = a.get<std::string>();
      if (s != "df" && s != "thc")
        throw std::runtime_error("unknown aux (want df|thc): " + s);
      c.aux.push_back(s);
    }
  }

  for (const auto& [label, size] : d.at("sizes").items()) {
    const auto sz = size.get<long long>();
    if (sz <= 0)
      throw std::runtime_error("index-space size must be > 0: " + label);
    c.sizes[label] = static_cast<std::size_t>(sz);
  }

  if (d.contains("optimize")) {
    const auto& o = d.at("optimize");
    c.optimize.objective = o.value("objective", std::string("dense_flops"));
    c.optimize.reorder = o.value("reorder", true);
    c.optimize.cse_subnet = o.value("cse_subnet", false);
    c.optimize.volatile_leaf = o.value("volatile_leaf", std::string("R"));
    c.optimize.volatile_weight = o.value("volatile_weight", 10.0);
    // A volatile-contraction weight is conceptually a replay count, so it must
    // be positive
    if (c.optimize.volatile_weight <= 0.0)
      throw Exception("optimize.volatile_weight must be > 0");
    c.optimize.machine_balance = o.value("machine_balance", 0.0);
    c.optimize.fast_mem_elems = o.value("fast_mem_elems", 0.0);
  }
  if (d.contains("cache")) {
    const auto& ca = d.at("cache");
    c.cache.enabled = ca.value("enabled", true);
    // A negative min_repeats wraps to a huge size_t (JSON ints are signed) and
    // silently disables caching, so reject it.
    const long mr = ca.value("min_repeats", 2L);
    if (mr < 0) throw std::runtime_error("cache.min_repeats must be >= 0");
    c.cache.min_repeats = static_cast<std::size_t>(mr);
    c.cache.max_footprint = ca.value("max_footprint", 0.0);
    if (c.cache.max_footprint < 0.0)
      throw std::runtime_error("cache.max_footprint must be >= 0");
  }
  if (d.contains("output")) {
    const auto& ou = d.at("output");
    c.out.path = ou.value("path", std::string("cost_analysis.md"));
    // Read signed then guard: a negative top_n would wrap to a huge size_t.
    const long tn = ou.value("top_n", 20L);
    if (tn < 0) throw std::runtime_error("output.top_n must be >= 0");
    c.out.top_n = static_cast<std::size_t>(tn);
    c.out.dump_tree = ou.value("dump_tree", false);
  }

  for (const auto& r : d.at("equations")) {
    EquationSpec rs;
    rs.name = r.at("name").get<std::string>();
    rs.equation_file = r.at("equation_file").get<std::string>();
    c.equations.push_back(std::move(rs));
  }
  return c;
}

std::shared_ptr<IndexSpaceRegistry> make_registry(const std::string& conv) {
  if (conv == "min_sr") return make_min_sr_spaces();
  if (conv == "sr") return make_sr_spaces();
  if (conv == "mr") return make_mr_spaces();
  if (conv == "f12") return make_F12_sr_spaces();
  throw std::runtime_error("unknown convention (want min_sr|sr|mr|f12): " +
                           conv);
}

// Start from one of SeQuant's standard registries, apply the config's field
// to every registered space, then override each named space's approximate
// size from the config.
void setup_context(const Config& cfg) {
  auto isr = make_registry(cfg.convention);
  // Register auxiliary factorization spaces (Κ for df, L for thc) so
  // pre-factorized DF/THC equations parse; validated to df|thc in load_config.
  for (const auto& a : cfg.aux) {
    if (a == "df")
      add_df_spaces(isr);
    else
      add_thc_spaces(isr);
  }

  // Field is a single global choice for the computation, so apply it to every
  // registry space (not just the sized ones) the way MPQC's
  // scoped_sequant_field does: bra<->ket symmetry must be consistent for all
  // used spaces, else a used-but-unsized space keeps the registry's default
  // (complex) and its tensors canonicalize wrongly over a real field. Collect
  // keys first (iteration is const) then set via mutable lookup.
  const Field field = cfg.real_field ? Field::Real : Field::Complex;
  std::vector<std::wstring> keys;
  for (const auto& sp : *isr) keys.push_back(sp.base_key());
  for (const auto& key : keys)
    if (IndexSpace* sp = isr->retrieve_ptr(key)) sp->field(field);

  for (const auto& [label, size] : cfg.sizes) {
    IndexSpace* sp = isr->retrieve_ptr(toUtf16(label));
    if (!sp) throw std::runtime_error("unknown index space: " + label);
    sp->approximate_size(size);
  }

  auto ctx = get_default_context();
  auto copts = CanonicalizeOptions::default_options().copy_and_set(
      CanonicalizationMethod::Complete);
  ctx.set(Vacuum::SingleProduct)
      .set(cfg.spinor ? SPBasis::Spinor : SPBasis::Spinfree)
      .set(isr)
      .set(BraKetSlotTypesetting::Naive)
      .set(copts);
  set_default_context(ctx);

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
}

ResultExpr parse_equation(const std::string& text) {
  return deserialize<ResultExpr>(std::string_view{text});
}

/// Cost Analysis

ObjectiveFunction objective_of(const std::string& s) {
  if (s == "dense_size") return ObjectiveFunction::DenseSize;
  if (s == "dense_peak_size") return ObjectiveFunction::DensePeakSize;
  if (s == "dense_flops") return ObjectiveFunction::DenseFLOPs;
  throw std::runtime_error(
      "unknown objective (want dense_size|dense_peak_size|dense_flops): " + s);
}

// Classify each index as Occupied / Virtual / other (aux) by matching its space
// Type against the registry's designated hole/particle types.
std::string space_signature(const EvalExpr& ev) {
  if (!ev.is_tensor()) return "scalar";
  const auto& t = ev.as_tensor();
  const auto& reg = *get_default_context().index_space_registry();
  const auto hole = reg.hole_space(/*nulltype_ok=*/true);
  const auto particle = reg.particle_space(/*nulltype_ok=*/true);
  std::string result;
  for (const auto& idx : t.const_braketaux_indices()) {
    const auto ty = idx.space().type();
    if (ty == hole)
      result += 'O';
    else if (ty == particle)
      result += 'V';
    else
      result += 'X';  // neither hole nor particle (e.g. df/thc aux)
  }
  return result;
}

AsyCost result_memory(const EvalExpr& ev) {
  if (!ev.is_tensor()) return AsyCost::zero();
  const auto& t = ev.as_tensor();
  AsyCost::ExponentMap exponents;
  // Include aux (braketaux, not braket) so auxiliary dimensions are counted in
  // the intermediate's memory/cost, consistent with space_signature.
  for (const auto& idx : t.const_braketaux_indices()) ++exponents[idx.space()];
  return AsyCost{std::move(exponents)};
}

// Like SeQuant's eval_node.hpp min_storage(), but maxes numerically (by
// .ops(), i.e. actual size at this driver's configured index-space sizes)
// instead of via AsyCost::operator<, which orders by symbolic polynomial
// degree. A lower-degree node (e.g. O^4) can be numerically larger than a
// higher-degree one (e.g. O*V*aux) at the sizes this tool is run with, so the
// symbolic max used by the library helper can pick the wrong peak.
AsyCost numeric_peak_storage(const TreeNode& tree) {
  AsyCost result = AsyCost::zero();
  tree.visit([&](const TreeNode& n) {
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

OptimizeOptions make_opts(const Config& cfg) {
  OptimizeOptions o;
  o.objective_function = objective_of(cfg.optimize.objective);
  o.reorder =
      cfg.optimize.reorder ? ReorderSum::Reorder : ReorderSum::NoReorder;
  o.CSE.subnet = cfg.optimize.cse_subnet;
  // volatile leaf (label match), MPQC-style volatile weighting.
  if (!cfg.optimize.volatile_leaf.empty()) {
    const std::wstring vl = toUtf16(cfg.optimize.volatile_leaf);
    o.batch_policy.is_volatile_leaf = [vl](const Tensor& t) {
      return t.label() == vl;
    };
    // MPQC: volatile_weight pinned to 1.0 when caching is off.
    o.volatile_weight = cfg.cache.enabled ? cfg.optimize.volatile_weight : 1.0;
    if (cfg.cache.enabled) {
      o.roofline.machine_balance = cfg.optimize.machine_balance;
      o.roofline.fast_mem_elems = cfg.optimize.fast_mem_elems;
    }
  }
  return o;
}

CellResult process(const Config& cfg, const ExprPtr& rhs, const Tensor& head) {
  CellResult res;

  ExprPtr input = rhs->clone();
  auto popped = pop_tensor(input, reserved::symm_label());
  if (!popped.has_value()) pop_tensor(input, reserved::antisymm_label());

  flatten(input);

  const auto opts = make_opts(cfg);

  auto scan = [&](const TreeNode& tree) {
    tree.visit_internal([&](const TreeNode& node) {
      if (node.root())
        return;  // per-term root is the output, not an intermediate
      auto it = res.catalog.find(node);
      if (it == res.catalog.end()) {
        Record rec;
        rec.uses = 1;
        rec.memory = result_memory(*node);
        rec.flops = Flops{}(node);
        rec.label = node->label();
        rec.spaces = space_signature(*node);
        res.catalog.emplace(node, std::move(rec));
      } else {
        ++it->second.uses;
      }
    });
  };

  auto handle = [&](const ExprPtr& summand) {
    ExprPtr os = optimize(summand, opts);
    auto tree = binarize<EvalExpr>(ResultExpr{head, os});
    scan(tree);
    // Peak working set = heaviest single-contraction over the schedule; report
    // the max across the term-by-term summand trees. Total FLOPs = their sum.
    const AsyCost ps = numeric_peak_storage(tree);
    if (res.peak_storage.ops() < ps.ops()) res.peak_storage = ps;
    res.total_flops = res.total_flops + asy_cost(tree);
    res.trees.push_back(std::move(tree));
  };

  if (input->is<Sum>()) {
    const auto& sum = input->as<Sum>();
    res.n_terms = sum.size();
    for (std::size_t i = 0; i < sum.size(); ++i) handle(sum.summand(i));
  } else {
    res.n_terms = 1;
    handle(input);
  }

  for (const auto& [_, rec] : res.catalog) {
    if (rec.memory == AsyCost::zero()) continue;
    ++res.n_distinct;
    if (res.largest_mem.ops() < rec.memory.ops()) res.largest_mem = rec.memory;
    if (rec.uses >= 2) ++res.n_reused;
  }
  return res;
}

/// Cache Simulation

bool is_leaf_labeled(const TreeNode& n, const std::wstring& label) {
  return n.leaf() && n->is_tensor() && n->as_tensor().label() == label;
}

SimResult simulate_cache(
    const Config& cfg,
    const std::vector<std::pair<std::string, CellResult>>& results) {
  std::vector<TreeNode> forest;
  for (const auto& [name, cell] : results)
    for (const auto& tree : cell.trees) forest.push_back(tree);

  const std::wstring vl = toUtf16(cfg.optimize.volatile_leaf);
  auto is_volatile = [vl](const TreeNode& n) { return is_leaf_labeled(n, vl); };
  auto cm = cache_manager(
      forest, is_volatile, cfg.cache.min_repeats,
      [](const TreeNode& n) { return result_memory(*n).ops(); },
      cfg.cache.max_footprint);

  SimResult r;
  cm.for_each_key([&](const TreeNode& n) {
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

std::string read_file(const std::filesystem::path& p) {
  std::ifstream in(p);
  if (!in) throw std::runtime_error("cannot open equation file: " + p.string());
  return std::string(std::istreambuf_iterator<char>(in), {});
}

}  // namespace

int main(int argc, char** argv) {
  set_locale();
  sequant::detail::OpIdRegistrar op_id_registrar;

  CLI::App app("SeQuant expression cost/factorization analysis");
  argv = app.ensure_utf8(argv);
  std::filesystem::path driver;
  app.add_option("--driver", driver, "Path to the JSON driver file")
      ->required();
  bool omit_revision = false;
  app.add_flag("--omit-revision", omit_revision,
               "Omit the git-revision line (for reproducible test diffs)");
  CLI11_PARSE(app, argc, argv);

  try {
    std::ifstream din(driver);
    if (!din) {
      std::cerr << "cannot open driver: " << driver << "\n";
      return 1;
    }
    json d = json::parse(din, nullptr, true, /*skip_comments=*/true);
    std::filesystem::current_path(
        std::filesystem::absolute(driver).parent_path());

    Config cfg = load_config(d);
    cfg.out.omit_revision = omit_revision;
    setup_context(cfg);  // registry baked before any parse

    std::vector<std::pair<std::string, CellResult>> results;
    for (const auto& r : cfg.equations) {
      const ResultExpr res = parse_equation(read_file(r.equation_file));
      CellResult cell = process(cfg, res.expression(), res.result_as_tensor());

      if (cfg.out.dump_tree) {
        const std::string dump_path = r.name + ".tree.txt";
        std::ofstream to(dump_path);
        if (!to)
          throw std::runtime_error("cannot open dump file: " + dump_path);
        for (const auto& tree : cell.trees) to << full_expr(tree) << "\n";
        to.close();
        if (!to)
          throw std::runtime_error("failed writing dump file: " + dump_path);
      }
      results.emplace_back(r.name, std::move(cell));
    }

    const SimResult sim =
        cfg.cache.enabled ? simulate_cache(cfg, results) : SimResult{};

    std::ofstream out(cfg.out.path);
    if (!out)
      throw std::runtime_error("cannot open output file: " + cfg.out.path);
    write_report(cfg, results, sim, out);
    out.close();
    if (!out)
      throw std::runtime_error("failed writing output file: " + cfg.out.path);
    // Print the absolute path: the working directory was switched to the
    // driver's parent above, so a relative out.path lands there, not in the
    // user's invocation directory.
    std::cout << "Report written to " << std::filesystem::absolute(cfg.out.path)
              << "\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
