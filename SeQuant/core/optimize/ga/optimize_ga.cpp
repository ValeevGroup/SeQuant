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
  if (opts.batch_policy.is_volatile_leaf)
    cm.volatile_weight = opts.volatile_weight;
  GAResult out;
  auto kt = build_key_table(targets, ixex, opts.batch_policy.is_volatile_leaf);
  Fitness fitness(kt, cm, resolution, ga_opts.memo_capacity);
  auto [flops, genome] = run_ga(fitness, seed_genome(kt), ga_opts);
  Schedule schedule = fitness.explain(genome);
  out.flops = flops;
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
  return out;
}

}  // namespace sequant
