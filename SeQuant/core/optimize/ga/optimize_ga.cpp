#include <SeQuant/core/optimize/ga/optimize_ga.hpp>

#include <SeQuant/core/space.hpp>

namespace sequant::opt::ga {

namespace {

// cost of one cluster's emitted subtree, every occurrence counted: emission
// renders each keyed array as its picked producer's route, so the walk
// follows pick at every internal node
double cl_tree_cost(KeyTable const& kt, Schedule const& sch,
                    CostModel const& cm, int d, NodeMask S) {
  if (std::popcount(S) == 1) return 0;
  if (auto it = sch.pick.find(kt.terms[d].key[S]); it != sch.pick.end()) {
    d = it->second.d;
    S = it->second.S;
  }
  auto const* cc = sch.forest.terms[d].children_of(S);
  return cm.merge(kt.terms[d], cc[0], cc[1]) +
         cl_tree_cost(kt, sch, cm, d, cc[0]) +
         cl_tree_cost(kt, sch, cm, d, cc[1]);
}

double no_sharing_cost(KeyTable const& kt, CostModel const& cm,
                       Schedule const& sch, ValPtr const& val) {
  Val const& v = *val;
  switch (v.kind) {
    case Val::Cl:
      return cl_tree_cost(kt, sch, cm, v.d, v.S);
    case Val::Fx:
      return no_sharing_cost(kt, cm, sch, v.inner) +
             cl_tree_cost(kt, sch, cm, v.d, v.V) +
             cm.merge(kt.terms[v.d], v.S, v.V);
    case Val::Sum:
      return no_sharing_cost(kt, cm, sch, v.s1.val) +
             no_sharing_cost(kt, cm, sch, v.s2.val) +
             cm.addition(kt.terms[v.d], v.S);
  }
  SEQUANT_UNREACHABLE_TOKEN;
}

}  // namespace

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
  Fitness fitness(kt, cm, resolution);
  auto [flops, genome] = run_ga(fitness, seed_genome(kt), ga_opts);
  out.schedule = fitness.explain(genome);
  out.genome = std::move(genome);
  out.flops = flops;
  for (auto const& root : out.schedule.roots)
    out.flops_no_sharing += no_sharing_cost(kt, cm, out.schedule, root.val);
  auto emission = emit_named(fitness, out.schedule);
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
