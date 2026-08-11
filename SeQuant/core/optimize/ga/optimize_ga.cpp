#include <SeQuant/core/optimize/ga/optimize_ga.hpp>

#include <SeQuant/core/space.hpp>

namespace sequant::opt::ga {

GAResult optimize_ga(container::svector<TargetInput> const& targets,
                     OptimizeOptions const& opts, GAOptions const& ga_opts,
                     CostModel const& cost, ProducerResolution resolution) {
  auto ixex = opts.idx_to_extent
                  ? opts.idx_to_extent
                  : index_to_extent_t{[](Index const& ix) -> std::size_t {
                      return ix.nonnull() ? ix.space().approximate_size() : 1;
                    }};
  // Replay weighting follows single-term optimization's rule (options.hpp):
  // consulted only when a volatile-leaf predicate is supplied, so callers that
  // pass none keep the historical volatility-blind cost bit for bit.
  CostModel cm = cost;
  if (opts.batch_policy.is_volatile_leaf) {
    cm.volatile_weight = opts.volatile_weight;
    // The matched pair (cost.hpp): widening the amortization class is only
    // meaningful once volatility is modelled, and it widens BOTH the
    // objective's replay rule and emission's naming rule at once.
    cm.amortize_persistent = opts.ga_amortize_persistent;
    cm.naming_cap_elems = opts.ga_naming_cap_elems;
  }
  GAResult out;
  auto kt = build_key_table(targets, ixex, opts.batch_policy.is_volatile_leaf);
  Fitness fitness(kt, cm, resolution, ga_opts.memo_capacity);
  SearchTrace trace;
  // The seed. Weighted only on request, and only where volatility is modelled
  // at all -- so a caller with no volatile-leaf predicate keeps the blind
  // per-term DP, and every number recorded against it, whatever it sets. The
  // two surfaces are a UNION, not a precedence: both default off, either turns
  // it on, and neither can turn the other off.
  const bool aware_seed =
      static_cast<bool>(opts.batch_policy.is_volatile_leaf) &&
      (ga_opts.volatility_aware_seed || opts.ga_volatility_aware_seed);
  auto [flops, genome] =
      run_ga(fitness, aware_seed ? seed_genome(kt, cm) : seed_genome(kt),
             ga_opts, &trace);
  Schedule schedule = fitness.explain(genome);
  out.flops = flops;
  out.seed_flops = trace.seed_cost;
  out.hill_climbed_flops = trace.hill_climbed_cost;
  out.restart_flops = std::move(trace.restart_costs);
  auto emission = emit_named(fitness, schedule);
  out.exprs = std::move(emission.targets);
  out.definitions = std::move(emission.definitions);
  return out;
}

}  // namespace sequant::opt::ga

namespace sequant {

GAOptimized optimize_ga(container::vector<ResultExpr> exprs,
                        OptimizeOptions const& opts,
                        opt::ga::GAOptions const& ga_opts) {
  using namespace opt::ga;
  container::svector<TargetInput> targets;
  for (std::size_t i = 0; i < exprs.size(); ++i) {
    auto const& r = exprs[i];
    TargetInput tgt;
    tgt.label = r.has_label() ? r.label() : L"Z" + std::to_wstring(i);
    FaceSet ext;
    for (auto const& ix : r.indices()) ext.emplace(ix);
    ExprPtr const& e = r.expression();
    if (e->is<Sum>())
      for (auto const& s : e->as<Sum>().summands()) tgt.summands.push_back(s);
    else
      tgt.summands.push_back(e);
    tgt.ext.assign(tgt.summands.size(), ext);
    targets.push_back(std::move(tgt));
  }
  auto result = optimize_ga(targets, opts, ga_opts);
  for (std::size_t i = 0; i < exprs.size(); ++i)
    exprs[i].expression() = result.exprs[i];
  GAOptimized out;
  out.targets = std::move(exprs);
  out.definitions.assign(result.definitions.begin(),
                         result.definitions.end());
  out.flops = result.flops;
  out.seed_flops = result.seed_flops;
  out.hill_climbed_flops = result.hill_climbed_flops;
  out.restart_flops.assign(result.restart_flops.begin(),
                           result.restart_flops.end());
  return out;
}

}  // namespace sequant
