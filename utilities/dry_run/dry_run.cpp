//
// dry_run: replay the planning of an evaluation -- optimize, binarize,
// schedule, meter -- from a JSON driver, with no numerical data. See README.md.
//

#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/meter.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/kramers_blind.hpp>
#include <SeQuant/core/eval/loop_space_key.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using nlohmann::json;
using namespace sequant;

namespace {

namespace dr = sequant::eval::dryrun;

// ---------------------------------------------------------------------------
// the driver
// ---------------------------------------------------------------------------

struct SpaceSpec {
  std::wstring key;
  std::optional<long> type, qns;  // the producer's bits, verified when given
  std::size_t size = 0;
};

struct BatchSpaceSpec {
  std::size_t target_size = 0;
  bool contracted = false, external = false;
};

struct EquationSpec {
  std::string name, head, rhs;
  // the producer's fingerprints of the trees it built for this block: the
  // ones it evaluates and, when it ran its own dry run, the ones it metered
  std::optional<std::string> tree_fingerprint, dry_tree_fingerprint;
};

struct MeterNumbers {
  double peak_bytes = 0, flops_persistent = 0, flops_volatile = 0,
         cost_persistent = 0, cost_volatile = 0;
  std::size_t builds = 0;
};

struct Driver {
  // context
  std::string convention = "mpqc-sr";
  bool vbs = false;
  bool complex_field = true;
  bool strict_braket_symmetry = true;
  std::size_t first_dummy_index_ordinal =
      Context::Defaults::first_dummy_index_ordinal;
  // sizes
  std::vector<SpaceSpec> spaces;
  std::map<std::wstring, container::svector<std::size_t>> slices;
  std::map<std::size_t, std::array<double, 5>> csv_moments;  // raw M_0..M_4
  bool csv_union_composites = false;
  // optimize: the options as the optimizer receives them
  bool optimize = true;
  std::string objective = "dense_flops";
  double volatile_weight = 1.0;
  bool cse_subnet = false;
  double peak_flops_tolerance = 0.10, footprint_weight = 0.0;
  bool prune_outer_products = true;
  bool term_batch_axes = false;
  RooflineParams roofline;
  // batch: the policy, its index predicates tabulated per space
  bool batch = false;
  std::string scheduler = "forest_descent";
  double peak_threshold = std::numeric_limits<double>::infinity();
  double accumulation_factor = 0.0;
  bool persistent_only = false, spectator_indices = false;
  std::map<std::wstring, BatchSpaceSpec> batch_spaces;  // by space base key
  // cache
  bool cache = true;
  std::size_t cache_min_repeats = 2;
  double cache_max_footprint = 0.0;  // elements; 0 = no cap
  // binarize
  std::wstring volatile_label = L"t";
  bool fold_kramers_eval_leaves = false;
  std::wstring kramers_blind_leaf_label;
  bool kramers_blind_loop_fusion = false;
  bool kramers_fold_intermediates = false;
  // equations
  std::vector<EquationSpec> equations;
  // what the producer got, when it recorded it
  std::optional<std::string> expected_fingerprint;
  std::optional<MeterNumbers> expected_meter;
};

double number_or_inf(const json& j) {
  // JSON has no infinity: null (or the string "inf") spells an unlimited value
  if (j.is_null()) return std::numeric_limits<double>::infinity();
  if (j.is_string()) {
    const auto s = j.get<std::string>();
    if (s == "inf" || s == "infinity")
      return std::numeric_limits<double>::infinity();
    return std::stod(s);
  }
  return j.get<double>();
}

void read_sizes(const json& d, Driver& c) {
  if (d.contains("spaces")) {
    // sizes from a second file override the sizes of the spaces it names and
    // leave the others (and every recorded type/qns) alone
    for (const auto& s : d.at("spaces")) {
      SpaceSpec sp;
      sp.key = toUtf16(s.at("key").get<std::string>());
      if (s.contains("type")) sp.type = s.at("type").get<long>();
      if (s.contains("qns")) sp.qns = s.at("qns").get<long>();
      const auto sz = s.at("size").get<long long>();
      if (sz < 0)
        throw Exception("index-space size must be >= 0: " + toUtf8(sp.key));
      sp.size = static_cast<std::size_t>(sz);
      auto it = std::ranges::find(c.spaces, sp.key, &SpaceSpec::key);
      if (it == c.spaces.end()) {
        c.spaces.push_back(std::move(sp));
      } else {
        it->size = sp.size;
        if (sp.type) it->type = sp.type;
        if (sp.qns) it->qns = sp.qns;
      }
    }
  }
  if (d.contains("slices")) {
    c.slices.clear();
    for (const auto& [key, arr] : d.at("slices").items()) {
      container::svector<std::size_t> v;
      for (const auto& x : arr) v.push_back(x.get<std::size_t>());
      c.slices[toUtf16(key)] = std::move(v);
    }
  }
  if (d.contains("csv")) {
    const auto& csv = d.at("csv");
    if (csv.contains("moments")) {
      c.csv_moments.clear();
      for (const auto& [rank, arr] : csv.at("moments").items()) {
        if (arr.size() != 5)
          throw Exception("csv.moments." + rank + " must hold M_0..M_4");
        std::array<double, 5> M{};
        for (std::size_t k = 0; k != 5; ++k) M[k] = arr[k].get<double>();
        c.csv_moments[std::stoul(rank)] = M;
      }
    }
    if (csv.contains("union_composites"))
      c.csv_union_composites = csv.at("union_composites").get<bool>();
  }
}

Driver load_driver(const json& d) {
  Driver c;
  if (d.value("format", std::string{}) != "sequant-dry-run")
    throw Exception("not a dry-run driver (\"format\" != \"sequant-dry-run\")");
  if (d.value("version", 0) != 1)
    throw Exception("unsupported driver version (want 1)");

  if (d.contains("context")) {
    const auto& ctx = d.at("context");
    c.convention = ctx.value("convention", c.convention);
    c.vbs = ctx.value("vbs", false);
    const auto field = ctx.value("field", std::string("complex"));
    if (field != "real" && field != "complex")
      throw Exception("unknown field (want real|complex): " + field);
    c.complex_field = field == "complex";
    c.strict_braket_symmetry = ctx.value("strict_braket_symmetry", true);
    c.first_dummy_index_ordinal =
        ctx.value("first_dummy_index_ordinal", c.first_dummy_index_ordinal);
  }

  read_sizes(d, c);

  if (d.contains("optimize")) {
    const auto& o = d.at("optimize");
    c.optimize = o.value("enabled", true);
    c.objective = o.value("objective_function", c.objective);
    c.volatile_weight = o.value("volatile_weight", c.volatile_weight);
    c.cse_subnet = o.value("cse_subnet", false);
    c.peak_flops_tolerance =
        o.value("peak_flops_tolerance", c.peak_flops_tolerance);
    c.footprint_weight = o.value("footprint_weight", c.footprint_weight);
    c.prune_outer_products = o.value("prune_outer_products", true);
    c.term_batch_axes = o.value("term_batch_axes", false);
    if (o.contains("roofline")) {
      const auto& r = o.at("roofline");
      c.roofline.machine_balance = r.value("machine_balance", 0.0);
      c.roofline.fast_mem_elems = r.value("fast_mem_elems", 0.0);
      c.roofline.block_tiles = r.value("block_tiles", c.roofline.block_tiles);
      c.roofline.block_prefactor =
          r.value("block_prefactor", c.roofline.block_prefactor);
    }
  }
  if (d.contains("batch")) {
    const auto& b = d.at("batch");
    c.batch = b.value("enabled", false);
    c.scheduler = b.value("scheduler", c.scheduler);
    if (c.scheduler != "ordered" && c.scheduler != "forest_descent")
      throw Exception(
          "unknown batch.scheduler (want ordered|forest_descent): " +
          c.scheduler);
    if (b.contains("peak_threshold"))
      c.peak_threshold = number_or_inf(b.at("peak_threshold"));
    c.accumulation_factor = b.value("accumulation_factor", 0.0);
    c.persistent_only = b.value("persistent_only", false);
    c.spectator_indices = b.value("spectator_indices", false);
    if (b.contains("spaces"))
      for (const auto& [key, sp] : b.at("spaces").items()) {
        BatchSpaceSpec bs;
        bs.target_size = sp.value("target_size", std::size_t{0});
        bs.contracted = sp.value("contracted", false);
        bs.external = sp.value("external", false);
        c.batch_spaces[toUtf16(key)] = bs;
      }
  }
  if (d.contains("cache")) {
    const auto& ca = d.at("cache");
    c.cache = ca.value("enabled", true);
    c.cache_min_repeats = ca.value("min_repeats", std::size_t{2});
    c.cache_max_footprint = ca.value("max_footprint", 0.0);
  }
  if (d.contains("binarize")) {
    const auto& b = d.at("binarize");
    c.volatile_label = toUtf16(b.value("volatile_label", std::string("t")));
    c.fold_kramers_eval_leaves = b.value("fold_kramers_eval_leaves", false);
    c.kramers_blind_leaf_label =
        toUtf16(b.value("kramers_blind_leaf_label", std::string{}));
    c.kramers_blind_loop_fusion = b.value("kramers_blind_loop_fusion", false);
    c.kramers_fold_intermediates = b.value("kramers_fold_intermediates", false);
  }
  for (const auto& e : d.at("equations")) {
    EquationSpec es;
    es.name = e.at("name").get<std::string>();
    es.head = e.at("head").get<std::string>();
    es.rhs = e.at("rhs").get<std::string>();
    if (e.contains("tree_fingerprint"))
      es.tree_fingerprint = e.at("tree_fingerprint").get<std::string>();
    if (e.contains("dry_tree_fingerprint"))
      es.dry_tree_fingerprint = e.at("dry_tree_fingerprint").get<std::string>();
    c.equations.push_back(std::move(es));
  }
  if (d.contains("expected")) {
    const auto& x = d.at("expected");
    if (x.contains("eval_tree_fingerprint"))
      c.expected_fingerprint = x.at("eval_tree_fingerprint").get<std::string>();
    if (x.contains("meter")) {
      const auto& m = x.at("meter");
      MeterNumbers n;
      n.peak_bytes = m.at("peak_bytes").get<double>();
      n.flops_persistent = m.at("flops_persistent").get<double>();
      n.flops_volatile = m.at("flops_volatile").get<double>();
      n.cost_persistent = m.at("cost_persistent").get<double>();
      n.cost_volatile = m.at("cost_volatile").get<double>();
      n.builds = m.at("builds").get<std::size_t>();
      c.expected_meter = n;
    }
  }
  return c;
}

// ---------------------------------------------------------------------------
// the context
// ---------------------------------------------------------------------------

// MPQC's single-reference LCAO convention
// (mpqc/math/external/sequant/sequant.cpp, make_sr_spaces). The driver lists
// the producer's spaces with their type and quantum-number bits, and
// setup_context() refuses a driver whose listing this registry does not
// reproduce, so a change on either side cannot pass unnoticed.
std::shared_ptr<IndexSpaceRegistry> make_mpqc_sr_spaces(bool vbs) {
  using namespace sequant::mbpt;
  auto isr = std::make_shared<IndexSpaceRegistry>();
  const auto spin_any = IndexSpace::QuantumNumbers{Spin::any};
  isr->add(L"i", 0b01, spin_any, is_vacuum_occupied, is_reference_occupied,
           is_hole)
      .add(L"a", 0b10, spin_any, is_particle)
      .add_union(L"p", {L"i", L"a"}, is_complete);
  add_fermi_spin(*isr);
  add_ao_spaces(isr, spin_any, vbs, /* abs = */ false);
  add_pao_spaces(isr, spin_any);
  {
    const auto pao = isr->retrieve(L"μ̃");
    const IndexSpace pao_up{
        spinannotation_add(pao.base_key(), Spin::alpha), pao.type(),
        IndexSpace::QuantumNumbers{Spin::alpha} | LCAOQNS::pao,
        pao.approximate_size()};
    const IndexSpace pao_dn{
        spinannotation_add(pao.base_key(), Spin::beta), pao.type(),
        IndexSpace::QuantumNumbers{Spin::beta} | LCAOQNS::pao,
        pao.approximate_size()};
    isr->add(pao_up);
    isr->add(pao_dn);
    isr->add_kramers_partners(pao_up, pao_dn);
  }
  add_df_spaces(isr);
  isr->physical_particle_attribute_mask(mask_v<Spin>);
  return isr;
}

void setup_context(const Driver& cfg) {
  std::shared_ptr<IndexSpaceRegistry> isr;
  if (cfg.convention == "mpqc-sr")
    isr = make_mpqc_sr_spaces(cfg.vbs);
  else
    throw Exception("unknown context.convention (want mpqc-sr): " +
                    cfg.convention);

  const Field field = cfg.complex_field ? Field::Complex : Field::Real;
  std::vector<std::wstring> keys;
  for (const auto& sp : *isr) keys.push_back(std::wstring(sp.base_key()));
  for (const auto& key : keys)
    if (IndexSpace* sp = isr->retrieve_ptr(key)) sp->field(field);

  for (const auto& s : cfg.spaces) {
    IndexSpace* sp = isr->retrieve_ptr(s.key);
    if (!sp)
      throw Exception("the driver lists index space \"" + toUtf8(s.key) +
                      "\", which the \"" + cfg.convention +
                      "\" convention does not have");
    if ((s.type && *s.type != static_cast<long>(sp->type().to_int32())) ||
        (s.qns && *s.qns != static_cast<long>(sp->qns().to_int32())))
      throw Exception(
          "index space \"" + toUtf8(s.key) +
          "\" has other type/quantum-number bits here than in the program "
          "that wrote the driver; the two conventions have diverged");
    sp->approximate_size(s.size);
  }

  auto ctx = Context()
                 .set(Vacuum::SingleProduct)
                 .set(IndexSpaceMetric::Unit)
                 .set(SPBasis::Spinor)
                 .set(CanonicalizeOptions::default_options().copy_and_set(
                     CanonicalizationMethod::Complete))
                 .set(isr);
  ctx.set_first_dummy_index_ordinal(cfg.first_dummy_index_ordinal);
  if (!cfg.strict_braket_symmetry) ctx.set(AssertStrictBraKetSymmetry::No);
  set_default_context(ctx);
}

// ---------------------------------------------------------------------------
// the planner inputs, assembled as the producer assembles them
// ---------------------------------------------------------------------------

/// the size of a composite's domain: the k-th power mean of the per-cluster
/// domain sizes, never below 2, and for a pair (or larger cluster) never below
/// the RMS (k = 2)
double csv_extent(const Driver& cfg, std::size_t rank, std::size_t k) {
  static const std::array<double, 5> ones{1, 1, 1, 1, 1};
  auto it = cfg.csv_moments.find(std::max<std::size_t>(rank, 1));
  if (it == cfg.csv_moments.end() && rank > 2) it = cfg.csv_moments.find(2);
  const auto& M = it != cfg.csv_moments.end() ? it->second : ones;
  std::size_t k_eff = rank >= 2 ? std::max<std::size_t>(k, 2) : k;
  k_eff = std::clamp<std::size_t>(k_eff, 1, 4);
  return std::max(2.0, M[k_eff]);
}

std::array<double, 5> csv_floored(const Driver& cfg, std::size_t rank) {
  std::array<double, 5> a{1, 1, 1, 1, 1};
  for (std::size_t k = 1; k <= 4; ++k) a[k] = csv_extent(cfg, rank, k);
  return a;
}

/// a spin-free composite spans the members of both Kramers flavours
double union_factor(const Driver& cfg, const Index& idx) {
  using namespace sequant::mbpt;
  return cfg.csv_union_composites && to_spin(idx.space().qns()) == Spin::any
             ? 2.0
             : 1.0;
}

std::function<std::size_t(const Index&)> make_idx_to_extent(const Driver& cfg) {
  return [&cfg](const Index& idx) -> std::size_t {
    if (idx.has_proto_indices())
      return static_cast<std::size_t>(
          std::ceil(union_factor(cfg, idx) *
                    csv_extent(cfg, idx.proto_indices().size(), 1)));
    return idx.space().approximate_size();
  };
}

std::function<double(const Index&, std::size_t)> make_inner_pow(
    const Driver& cfg) {
  return [&cfg](const Index& idx, std::size_t k) -> double {
    return union_factor(cfg, idx) *
           csv_extent(cfg, idx.proto_indices().size(), k);
  };
}

std::function<bool(const Tensor&)> label_is(const std::wstring& label) {
  if (label.empty()) return {};
  return [label](const Tensor& t) { return t.label() == label; };
}

/// the batch policy the evaluation runs under: the default (nothing batches)
/// unless the driver's is in force
BatchPolicy make_batch_policy(const Driver& cfg) {
  BatchPolicy p;
  if (!cfg.batch) return p;
  auto table = std::make_shared<const std::map<std::wstring, BatchSpaceSpec>>(
      cfg.batch_spaces);
  auto resolve = [table](const Index& ix) -> const BatchSpaceSpec* {
    if (!ix.space()) return nullptr;
    auto it = table->find(std::wstring(ix.space().base_key()));
    return it == table->end() ? nullptr : &it->second;
  };
  p.is_batchable_contracted_index = [resolve](const Index& ix) {
    const auto* s = resolve(ix);
    return s && s->contracted;
  };
  p.is_batchable_external_index = [resolve](const Index& ix) {
    const auto* s = resolve(ix);
    return s && s->external;
  };
  p.batch_target_size = [resolve](const Index& ix) -> std::size_t {
    const auto* s = resolve(ix);
    return s ? s->target_size : std::size_t{0};
  };
  p.batch_spectator_indices = cfg.spectator_indices;
  p.is_volatile_leaf = label_is(cfg.volatile_label);
  p.persistent_only = cfg.persistent_only;
  p.accumulation_factor = cfg.accumulation_factor;
  p.peak_threshold = cfg.peak_threshold;
  p.scheduler = cfg.scheduler == "ordered" ? BatchScheduler::ordered
                                           : BatchScheduler::forest_descent;
  return p;
}

ObjectiveFunction objective_of(const std::string& s) {
  if (s == "dense_flops") return ObjectiveFunction::DenseFLOPs;
  if (s == "dense_size") return ObjectiveFunction::DenseSize;
  if (s == "dense_space_time") return ObjectiveFunction::DenseSpaceTime;
  if (s == "dense_space_time_batched")
    return ObjectiveFunction::DenseSpaceTimeBatched;
  if (s == "dense_time_space") return ObjectiveFunction::DenseTimeSpace;
  if (s == "dense_time_space_batched")
    return ObjectiveFunction::DenseTimeSpaceBatched;
  throw Exception(
      "unknown optimize.objective_function (want dense_flops|dense_size|"
      "dense_space_time[_batched]|dense_time_space[_batched]): " +
      s);
}

OptimizeOptions make_optimize_options(const Driver& cfg,
                                      const BatchPolicy& policy) {
  OptimizeOptions opts;
  opts.objective_function = objective_of(cfg.objective);
  opts.reorder = ReorderSum::Reorder;
  opts.CSE.subnet = cfg.cse_subnet;
  opts.idx_to_extent = make_idx_to_extent(cfg);
  opts.inner_pow = make_inner_pow(cfg);
  opts.batch_policy = policy;
  // an unbatched evaluation carries no policy, yet the optimizer still
  // weights the contractions that depend on the volatile leaf
  if (!opts.batch_policy.is_volatile_leaf)
    opts.batch_policy.is_volatile_leaf = label_is(cfg.volatile_label);
  opts.volatile_weight = cfg.volatile_weight;
  opts.peak_flops_tolerance = cfg.peak_flops_tolerance;
  opts.footprint_weight = cfg.footprint_weight;
  opts.prune_outer_products = cfg.prune_outer_products;
  opts.roofline = cfg.roofline;
  if (cfg.term_batch_axes) {
    using BatchAxesVec = container::vector<NodeBatchAnnotation>;
    using BatchAxesMap = std::unordered_map<Expr const*, BatchAxesVec>;
    opts.term_batch_axes = std::make_shared<BatchAxesMap>();
  }
  return opts;
}

/// node identity that is blind to the Kramers flavour of the outer
/// pure-occupied slots of the leaf labelled \p label (the CSV projector, which
/// one array serves for every pair flavour)
eval::KramersBlindness make_kramers_blindness(const std::wstring& label) {
  const auto isr = get_default_context().index_space_registry();
  return eval::KramersBlindness{
      .blind_slot =
          [label, isr](const Tensor& t, std::size_t slot) {
            if (t.label() != label) return false;
            const auto& ix = *(t.const_slots().begin() + slot);
            bool blind = false;
            if (ix.has_proto_indices()) {
              blind = !ix.proto_indices().empty();
              for (const auto& p : ix.proto_indices())
                blind = blind && isr->is_pure_occupied(p.space());
            } else {
              blind = isr->is_pure_occupied(ix.space());
            }
            return blind;
          },
      .erase_space =
          [](const IndexSpace& s) {
            using namespace sequant::mbpt;
            if ((s.qns().to_int32() & mask_v<Spin>) == 0) return s;
            if (to_spin(s.qns()) == Spin::any) return s;
            return make_spinalpha(Index(s, 1)).space();
          }};
}

dr::SizeRegime make_size_regime(const Driver& cfg,
                                const std::vector<dr::EvalNodeDryRun>& forest) {
  dr::SizeRegime regime;
  regime.csv_osv_moment = csv_floored(cfg, 1);
  regime.csv_pno_moment = csv_floored(cfg, 2);
  for (const auto& [rank, M] : cfg.csv_moments)
    if (rank >= 3) regime.csv_moment_by_rank[rank] = csv_floored(cfg, rank);
  regime.space_slice_extents = cfg.slices;
  const auto idx_to_extent = make_idx_to_extent(cfg);
  for (const auto& node : forest)
    node.visit([&](const auto& n) {
      if (!n.leaf() || !(*n).is_tensor()) return;
      for (const auto& ix : (*n).canon_indices())
        if (!ix.has_proto_indices())
          regime.space_extent[std::wstring(ix.space().base_key())] =
              idx_to_extent(ix);
    });
  return regime;
}

double bytes_per_element(const Driver& cfg) {
  return 8.0 * sequant::opt::detail::field_cost_factors(
                   cfg.complex_field ? Field::Complex : Field::Real)
                   .elem_mult;
}

// ---------------------------------------------------------------------------
// one equation: optimize and binarize it
// ---------------------------------------------------------------------------

std::size_t eval_node_fingerprint(const dr::EvalNodeDryRun& root) {
  std::size_t fp = 0;
  root.visit(
      [&fp](const dr::EvalNodeDryRun& n) {
        hash::combine(fp, n->hash_value());
        hash::combine(fp, n->layout_fingerprint());
        hash::combine(fp, static_cast<int>(n->canon_phase()));
        hash::combine(fp, n.leaf());
      },
      TreeTraversal::PreOrder);
  return fp;
}

struct Planned {
  std::vector<dr::EvalNodeDryRun> roots;
};

Planned plan_equation(const Driver& cfg, const BatchPolicy& policy,
                      const Tensor& head, const ExprPtr& expr) {
  auto input = expr->clone();
  auto popped = pop_tensor(input, reserved::symm_label());
  if (!popped.has_value()) pop_tensor(input, reserved::antisymm_label());

  auto to_node = [&](const ExprPtr& e) -> dr::EvalNodeDryRun {
    std::string msg;
    if (!is_valid(e, &msg)) throw Exception("invalid expression (" + msg + ")");
    ExprPtr rhs = e;
    BinarizationOptions bopts;
    bopts.kramers_fold_intermediates = cfg.kramers_fold_intermediates;
    if (!cfg.kramers_blind_leaf_label.empty()) {
      bopts.kramers_blindness =
          make_kramers_blindness(cfg.kramers_blind_leaf_label);
      eval::set_loop_space_erasure(cfg.kramers_blind_loop_fusion
                                       ? bopts.kramers_blindness.erase_space
                                       : nullptr);
    }
    if (cfg.optimize) {
      auto opts = make_optimize_options(cfg, policy);
      rhs = optimize_result(e, opts).expr;
      if (opts.term_batch_axes) {
        auto it = opts.term_batch_axes->find(rhs.get());
        if (it != opts.term_batch_axes->end())
          bopts.node_batch_axes = it->second;
      }
    }
    return binarize<dr::EvalExprDryRun>(ResultExpr{head, std::move(rhs)},
                                        bopts);
  };

  auto build = [&]() {
    Planned p;
    const bool single_root = policy.scheduler == BatchScheduler::ordered;
    if (!single_root && input->is<Sum>())
      for (const auto& s : input->as<Sum>()) p.roots.push_back(to_node(s));
    else
      p.roots.push_back(to_node(input));
    return p;
  };
  if (cfg.fold_kramers_eval_leaves) {
    auto fold_ctx = Context(get_default_context());
    fold_ctx.set(CanonicalizeOptions::default_options().copy_and_set(
        CanonicalizeOptions::FoldKramersEvalLeaves::Yes));
    auto scope = set_scoped_default_context(std::move(fold_ctx));
    return build();
  }
  return build();
}

// ---------------------------------------------------------------------------
// the report
// ---------------------------------------------------------------------------

/// an op's class: its trace label with index ordinals and Kramers flavours
/// dropped, so the instances of one contraction over blocks and batches add up
std::string op_class(std::string s, bool keep_flavours) {
  std::string out;
  out.reserve(s.size());
  for (std::size_t i = 0; i < s.size();) {
    if (s[i] == '_' && i + 1 < s.size() &&
        std::isdigit(static_cast<unsigned char>(s[i + 1]))) {
      ++i;
      while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i])))
        ++i;
      continue;
    }
    // U+2191 / U+2193 (the up / down flavour marks), UTF-8 encoded
    if (!keep_flavours && i + 2 < s.size() &&
        static_cast<unsigned char>(s[i]) == 0xE2 &&
        static_cast<unsigned char>(s[i + 1]) == 0x86 &&
        (static_cast<unsigned char>(s[i + 2]) == 0x91 ||
         static_cast<unsigned char>(s[i + 2]) == 0x93)) {
      i += 3;
      continue;
    }
    out.push_back(s[i++]);
  }
  return out;
}

struct ClassRow {
  std::size_t count = 0;
  double flops = 0, exec = 0, result_bytes = 0;
  bool volatile_ = false;
  std::size_t first_seen = 0;
};

/// sums the per-op records of the meter's trace (an `OpCost | flops | exec`
/// line, then the op's `Eval | Product | ...` line) by op class
std::vector<std::pair<std::string, ClassRow>> classify_ops(
    std::istream& trace,
    const std::unordered_map<std::string, bool>& volatile_of_label,
    bool keep_flavours) {
  std::map<std::string, ClassRow> rows;
  std::string line;
  double flops = 0, exec = 0;
  bool pending = false;
  std::size_t n = 0;
  auto split = [](const std::string& s) {
    std::vector<std::string> f;
    std::size_t b = 0;
    for (std::size_t e; (e = s.find(" | ", b)) != std::string::npos; b = e + 3)
      f.push_back(s.substr(b, e - b));
    f.push_back(s.substr(b));
    return f;
  };
  while (std::getline(trace, line)) {
    if (line.starts_with("OpCost | ")) {
      const auto f = split(line);
      if (f.size() >= 3) {
        flops = std::stod(f[1]);
        exec = std::stod(f[2]);
        pending = true;
      }
      continue;
    }
    if (!line.starts_with("Eval | Product | ")) continue;
    const auto f = split(line);
    std::string label, result;
    for (const auto& x : f) {
      if (const auto p = x.find(" canon="); p != std::string::npos)
        label = x.substr(0, p);
      if (x.starts_with("result=")) result = x.substr(7);
    }
    if (label.empty()) continue;
    auto& row = rows[op_class(label, keep_flavours)];
    if (row.count == 0) row.first_seen = n++;
    ++row.count;
    if (pending) {
      row.flops += flops;
      row.exec += exec;
    }
    if (!result.empty())
      row.result_bytes = std::max(row.result_bytes, std::stod(result));
    if (auto it = volatile_of_label.find(label); it != volatile_of_label.end())
      row.volatile_ = row.volatile_ || it->second;
    pending = false;
  }
  std::vector<std::pair<std::string, ClassRow>> sorted(rows.begin(),
                                                       rows.end());
  std::ranges::sort(sorted, [](const auto& a, const auto& b) {
    return a.second.flops > b.second.flops;
  });
  return sorted;
}

std::string human_bytes(double b) {
  if (b >= 1e12) return std::format("{:.2f} TB", b / 1e12);
  if (b >= 1e9) return std::format("{:.2f} GB", b / 1e9);
  if (b >= 1e6) return std::format("{:.2f} MB", b / 1e6);
  return std::format("{:.0f} B", b);
}

bool close(double a, double b) {
  // the meter line a producer prints carries ~12 significant digits
  return a == b || std::abs(a - b) <= 1e-9 * std::max(std::abs(a), std::abs(b));
}

}  // namespace

int main(int argc, char** argv) {
  set_locale();
  sequant::detail::OpIdRegistrar op_id_registrar;

  CLI::App app(
      "SeQuant dry run: plan an evaluation from a JSON driver and meter it");
  argv = app.ensure_utf8(argv);
  std::filesystem::path driver_path, sizes_path, trace_path;
  std::string peak_arg;
  std::size_t top = 12;
  app.add_option("--driver", driver_path, "the JSON driver")->required();
  app.add_option("--sizes", sizes_path,
                 "take the sizes (spaces, slices, csv) from this driver or "
                 "sizes-only JSON file instead");
  app.add_option("--peak-threshold", peak_arg,
                 "override batch.peak_threshold (bytes; \"inf\" = none)");
  app.add_option("--trace", trace_path,
                 "write the per-op trace of the metered replay here");
  app.add_option("--top", top, "rows of the op-class table");
  bool keep_flavours = false, print_equation_fingerprints = false;
  app.add_flag("--equation-fingerprints", print_equation_fingerprints,
               "print the tree fingerprint of every equation");
  std::string dump_tree_of;
  app.add_option("--dump-tree", dump_tree_of,
                 "print the planned tree of the equation with this name, one "
                 "node per line, children indented");
  app.add_flag("--flavours", keep_flavours,
               "keep the Kramers flavours apart in the op-class table");
  CLI11_PARSE(app, argc, argv);

  try {
    auto read_json = [](const std::filesystem::path& p) {
      std::ifstream in(p);
      if (!in) throw Exception("cannot open " + p.string());
      return json::parse(in, nullptr, true, /*skip_comments=*/true);
    };
    Driver cfg = load_driver(read_json(driver_path));
    bool overridden = false;
    if (!sizes_path.empty()) {
      read_sizes(read_json(sizes_path), cfg);
      overridden = true;
    }
    if (!peak_arg.empty()) {
      cfg.peak_threshold = number_or_inf(json(peak_arg));
      overridden = true;
    }
    setup_context(cfg);

    const BatchPolicy policy = make_batch_policy(cfg);

    // The equations carry dummy indices their producer generated, whose
    // ordinals a reader rejects by default: read them with that bound lifted,
    // then start this process's own generated ordinals past every one read.
    namespace ser = sequant::io::serialization;
    struct Equation {
      std::string name;
      ExprPtr head, rhs;
    };
    std::vector<Equation> equations;
    {
      auto lifted = Context(get_default_context());
      lifted.set_first_dummy_index_ordinal(
          std::numeric_limits<std::size_t>::max());
      const auto scope = set_scoped_default_context(std::move(lifted));
      for (const auto& eq : cfg.equations) {
        auto head = ser::from_string<ExprPtr>(eq.head);
        if (!head->is<Tensor>())
          throw Exception("equation \"" + eq.name +
                          "\": the head is not a tensor");
        equations.push_back(
            {eq.name, std::move(head), ser::from_string<ExprPtr>(eq.rhs)});
      }
    }
    std::size_t max_ordinal = 0;
    {
      auto note = [&max_ordinal](const Index& ix) {
        if (ix.ordinal())
          max_ordinal =
              std::max(max_ordinal, static_cast<std::size_t>(*ix.ordinal()));
        for (const auto& p : ix.proto_indices())
          if (p.ordinal())
            max_ordinal =
                std::max(max_ordinal, static_cast<std::size_t>(*p.ordinal()));
      };
      auto walk = [&note](const ExprPtr& e, auto& self) -> void {
        if (e->is<Tensor>())
          for (const auto& ix : e->as<Tensor>().const_braketaux_indices())
            note(ix);
        for (const auto& sub : *e) self(sub, self);
      };
      for (const auto& eq : equations) {
        walk(eq.head, walk);
        walk(eq.rhs, walk);
      }
    }
    Index::reset_tmp_index();
    if (max_ordinal >= Index::min_tmp_index())
      while (Index::next_tmp_index() < max_ordinal) {
      }

    std::vector<dr::EvalNodeDryRun> forest;
    std::size_t fingerprint = 0, n_roots = 0;
    std::vector<std::size_t> equation_fingerprints;
    for (const auto& eq : equations) {
      auto planned = plan_equation(cfg, policy, eq.head->as<Tensor>(), eq.rhs);
      // the producer orders each equation's roots by their fingerprint
      std::vector<std::size_t> fps;
      for (const auto& r : planned.roots)
        fps.push_back(eval_node_fingerprint(r));
      std::ranges::sort(fps);
      std::size_t fp_eq = 0;
      for (auto f : fps) hash::combine(fp_eq, f);
      hash::combine(fingerprint, fp_eq);
      equation_fingerprints.push_back(fp_eq);
      if (eq.name == dump_tree_of)
        for (const auto& root : planned.roots) {
          auto print = [](const dr::EvalNodeDryRun& n, std::size_t depth,
                          auto& self) -> void {
            std::string sliced;
            for (const auto& [ix, mode] : n->node_slice_mask())
              sliced +=
                  std::format(" {}{}", toUtf8(ix.full_label()),
                              mode == BatchModeType::External ? "(ext)" : "");
            std::cout << std::string(2 * depth, ' ') << n->label()
                      << (sliced.empty() ? "" : "   sliced" + sliced) << "\n";
            if (n.leaf()) return;
            self(n.left(), depth + 1, self);
            self(n.right(), depth + 1, self);
          };
          print(root, 0, print);
        }
      n_roots += planned.roots.size();
      // a scalar root (an observable) is not part of the summed replay
      for (auto& r : planned.roots)
        if (!r->canon_indices().empty()) forest.push_back(std::move(r));
    }

    std::cout << std::format("equations: {}  roots: {}  metered trees: {}\n",
                             cfg.equations.size(), n_roots, forest.size());
    std::cout << std::format("eval trees: fingerprint {:#018x}\n", fingerprint);
    if (print_equation_fingerprints)
      for (std::size_t e = 0; e != cfg.equations.size(); ++e)
        std::cout << std::format("  {:#018x}  {}\n", equation_fingerprints[e],
                                 cfg.equations[e].name);

    const dr::SizeRegime regime = make_size_regime(cfg, forest);
    dr::CacheConfig cache_cfg;
    cache_cfg.min_repeats = cfg.cache_min_repeats;
    cache_cfg.max_footprint = cfg.cache_max_footprint * bytes_per_element(cfg);
    const auto vpred = policy.is_volatile_leaf ? policy.is_volatile_leaf
                                               : label_is(cfg.volatile_label);
    if (vpred)
      cache_cfg.is_volatile = [vpred](const dr::EvalNodeDryRun& n) {
        if (!n.leaf() || !n->is_tensor()) return false;
        return vpred(n->as_tensor());
      };

    std::stringstream trace;
    const auto rep =
        dr::meter(forest, policy, regime, cache_cfg, &trace, cfg.roofline);
    std::cout << "dry-run meter ["
              << (rep.scheduler == BatchScheduler::ordered ? "ordered"
                                                           : "forest")
              << ", combined " << forest.size()
              << " trees]: peak=" << rep.peak_bytes / 1e9
              << " GB (dense model) | flops persistent=" << rep.flops_persistent
              << " volatile=" << rep.flops_volatile
              << " | time persistent=" << rep.cost_persistent
              << " volatile=" << rep.cost_volatile
              << " | builds=" << rep.builds_total << "\n";

    if (!trace_path.empty()) {
      std::ofstream out(trace_path);
      if (!out) throw Exception("cannot open " + trace_path.string());
      out << trace.str();
    }

    // which ops depend on the volatile leaf, by trace label
    std::unordered_map<std::string, bool> volatile_of_label;
    if (cache_cfg.is_volatile) {
      const auto volatile_of =
          dr::compute_volatility(forest, cache_cfg.is_volatile);
      for (const auto& root : forest)
        root.visit([&](const dr::EvalNodeDryRun& n) {
          if (n.leaf()) return;
          auto it = volatile_of.find(n->hash_value());
          const bool v = it != volatile_of.end() && it->second;
          auto [pos, fresh] =
              volatile_of_label.emplace(sequant::log::label(n), v);
          if (!fresh) pos->second = pos->second || v;
        });
    }
    const auto rows = classify_ops(trace, volatile_of_label, keep_flavours);
    double total_flops = 0;
    for (const auto& [name, row] : rows) total_flops += row.flops;
    std::cout << "\nproducts by class (flops summed over every execution of "
                 "the replay; * = depends on "
              << toUtf8(cfg.volatile_label) << "):\n";
    std::cout << std::format("  {:>10}  {:>6}  {:>7}  {:>11}  {}\n", "flops",
                             "share", "runs", "max result", "class");
    std::size_t shown = 0;
    for (const auto& [name, row] : rows) {
      if (shown++ == top) break;
      std::cout << std::format(
          "  {:>10.3e}  {:>5.1f}%  {:>7}  {:>11}  {}{}\n", row.flops,
          total_flops > 0 ? 100 * row.flops / total_flops : 0.0, row.count,
          human_bytes(row.result_bytes), row.volatile_ ? "* " : "  ", name);
    }
    if (rows.size() > top)
      std::cout << std::format("  ... {} more classes\n", rows.size() - top);

    // compare with what the producer recorded
    int rc = 0;
    if (overridden) {
      if (cfg.expected_fingerprint || cfg.expected_meter)
        std::cout << "\n(sizes or budget overridden: the producer's recorded "
                     "results do not apply)\n";
    } else {
      // equation by equation, against the trees the producer metered when
      // it recorded them, else against the trees it evaluates
      std::size_t n_compared = 0, n_same = 0;
      std::vector<std::string> differing;
      bool against_dry = true;
      for (std::size_t e = 0; e != cfg.equations.size(); ++e) {
        const auto& spec = cfg.equations[e];
        const auto& theirs = spec.dry_tree_fingerprint
                                 ? spec.dry_tree_fingerprint
                                 : spec.tree_fingerprint;
        if (!theirs) continue;
        if (!spec.dry_tree_fingerprint) against_dry = false;
        ++n_compared;
        if (std::format("{:#018x}", equation_fingerprints[e]) == *theirs)
          ++n_same;
        else
          differing.push_back(spec.name);
      }
      if (n_compared) {
        std::cout << "\nproducer's " << (against_dry ? "dry-run" : "eval")
                  << " trees, by equation: " << n_same << " of " << n_compared
                  << " identical"
                  << (n_same == n_compared ? "  MATCH" : "  MISMATCH") << "\n";
        for (const auto& name : differing)
          std::cout << "  differs: " << name << "\n";
        if (n_same != n_compared) rc = 2;
      } else if (cfg.expected_fingerprint) {
        const auto mine = std::format("{:#018x}", fingerprint);
        const bool same = mine == *cfg.expected_fingerprint;
        std::cout << "\nproducer's eval trees: " << *cfg.expected_fingerprint
                  << (same ? "  MATCH" : "  MISMATCH") << "\n";
        if (!same) rc = 2;
      }
      if (cfg.expected_meter) {
        const auto& x = *cfg.expected_meter;
        const bool same = close(x.peak_bytes, rep.peak_bytes) &&
                          close(x.flops_persistent, rep.flops_persistent) &&
                          close(x.flops_volatile, rep.flops_volatile) &&
                          close(x.cost_persistent, rep.cost_persistent) &&
                          close(x.cost_volatile, rep.cost_volatile) &&
                          x.builds == rep.builds_total;
        std::cout << "producer's meter: peak=" << x.peak_bytes / 1e9
                  << " GB | flops persistent=" << x.flops_persistent
                  << " volatile=" << x.flops_volatile
                  << " | builds=" << x.builds
                  << (same ? "  MATCH" : "  MISMATCH") << "\n";
        if (!same) rc = 2;
      }
    }
    return rc;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
