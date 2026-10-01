#include <SeQuant/core/wick_extended.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/core/wick.hpp>

#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <string>

namespace sequant::detail {

namespace {

using Blocks = container::svector<container::svector<std::size_t>>;

/// the core (R minus U), active (R ∩ U) and virtual (U minus R) parts of
/// the space, where R is the reference-occupied and U the vacuum-unoccupied
/// space
struct SpaceParts {
  IndexSpace::Type core, active, virt;
};

SpaceParts space_parts(const IndexSpaceRegistry &isr,
                       IndexSpace::QuantumNumbers qns) {
  const auto r = isr.reference_occupied_space(qns).type();
  const auto u = isr.vacuum_unoccupied_space(qns).type();
  const auto active = r.intersection(u);
  return {r.xOr(active), active, u.xOr(active)};
}

/// enumerates every way to group the ops of @p survivors into disjoint
/// cumulant blocks and (unless @p full) a remainder, and reports each via
/// @p sink(sign, blocks, remainder); blocks and remainder hold op positions
/// in storage order
template <Statistics S, typename Sink>
void for_each_block_assignment(const NormalOperator<S> &survivors,
                               const OpProvenance &provenance,
                               std::size_t max_rank, bool full, Sink &&sink) {
  const std::size_t n = survivors.size();
  constexpr std::size_t npos = static_cast<std::size_t>(-1);
  // assignment[i] = id of the block holding op i, or npos
  container::svector<std::size_t> assignment(n, npos);
  Blocks blocks;

  auto is_cre = [&](std::size_t i) {
    return survivors[i].action() == Action::Create;
  };
  const auto &isr = *get_default_context(S).index_space_registry();
  // only active ops can be cumulant legs
  auto is_active = [&](std::size_t i) {
    const auto &sp = survivors[i].index().space();
    return space_parts(isr, sp.qns()).active.includes(sp.type());
  };
  auto prov = [&](std::size_t i) {
    auto it = provenance.find(survivors[i].index());
    if (it == provenance.end())
      throw Exception("cumulant_expand: surviving operator index " +
                      toUtf8(survivors[i].index().full_label()) +
                      " has no provenance");
    return it->second;
  };

  auto emit = [&]() {
    container::svector<std::size_t> remainder;
    for (std::size_t i = 0; i != n; ++i)
      if (assignment[i] == npos) remainder.push_back(i);
    if (full && !remainder.empty()) return;
    // sign = parity of [block0 legs][block1 legs]...[remainder] relative to
    // storage order
    int sign = 1;
    if constexpr (S == Statistics::FermiDirac) {
      container::svector<std::size_t> perm;
      perm.reserve(n);
      for (const auto &b : blocks) perm.insert(perm.end(), b.begin(), b.end());
      perm.insert(perm.end(), remainder.begin(), remainder.end());
      for (std::size_t i = 0; i != n; ++i)
        for (std::size_t j = i + 1; j != n; ++j)
          if (perm[i] > perm[j]) sign = -sign;
    }
    sink(sign, blocks, remainder);
  };

  // decides the fate of the first undecided op at or after `first`; each
  // block is generated once, from its lowest leg
  std::function<void(std::size_t)> recurse = [&](std::size_t first) {
    while (first != n && assignment[first] != npos) ++first;
    if (first == n) {
      emit();
      return;
    }
    // `first` stays in the remainder
    if (!full) recurse(first + 1);
    // `first` is the lowest leg of a new block
    if (max_rank < 2 || !is_active(first)) return;
    const std::size_t block_id = blocks.size();
    container::svector<std::size_t> cre_cands, ann_cands;
    for (std::size_t i = first + 1; i != n; ++i)
      if (assignment[i] == npos && is_active(i))
        (is_cre(i) ? cre_cands : ann_cands).push_back(i);
    const auto &same = is_cre(first) ? cre_cands : ann_cands;
    const auto &other = is_cre(first) ? ann_cands : cre_cands;
    // pick k-1 more legs of first's kind and k of the other kind
    for (std::size_t k = 2; k <= max_rank; ++k) {
      if (same.size() + 1 < k || other.size() < k) break;
      container::svector<bool> sel_same(same.size(), false);
      std::fill(sel_same.begin(), sel_same.begin() + (k - 1), true);
      do {
        container::svector<bool> sel_other(other.size(), false);
        std::fill(sel_other.begin(), sel_other.begin() + k, true);
        do {
          container::svector<std::size_t> block{first};
          for (std::size_t i = 0; i != same.size(); ++i)
            if (sel_same[i]) block.push_back(same[i]);
          for (std::size_t i = 0; i != other.size(); ++i)
            if (sel_other[i]) block.push_back(other[i]);
          std::sort(block.begin(), block.end());
          // a block within a single GNO string vanishes
          const auto p0 = prov(block[0]);
          if (std::any_of(block.begin(), block.end(),
                          [&](std::size_t i) { return prov(i) != p0; })) {
            for (auto i : block) assignment[i] = block_id;
            blocks.push_back(block);
            recurse(first + 1);
            blocks.pop_back();
            for (auto i : block) assignment[i] = npos;
          }
        } while (std::prev_permutation(sel_other.begin(), sel_other.end()));
      } while (std::prev_permutation(sel_same.begin(), sel_same.end()));
    }
  };
  recurse(0);
}

/// @return @p ops, given in storage order, as a MultiProduct-vacuum
/// NormalOperator
template <Statistics S>
NormalOperator<S> make_nop(const container::svector<Op<S>> &ops) {
  container::svector<Op<S>> cre_ops, ann_ops;
  for (const auto &op : ops)
    (op.action() == Action::Create ? cre_ops : ann_ops).push_back(op);
  // the ctor takes annihilators in particle order, the reverse of storage
  std::reverse(ann_ops.begin(), ann_ops.end());
  return NormalOperator<S>(cre(std::move(cre_ops)), ann(std::move(ann_ops)),
                           Vacuum::MultiProduct);
}

/// @return the ops of @p nop at storage positions @p positions as a
/// MultiProduct-vacuum NormalOperator
template <Statistics S>
NormalOperator<S> subset(const NormalOperator<S> &nop,
                         const container::svector<std::size_t> &positions) {
  container::svector<Op<S>> ops;
  for (auto i : positions) ops.push_back(nop[i]);
  return make_nop<S>(ops);
}

/// @return the index of every Op of @p expr (a Product, NormalOperator<S> or
/// NormalOperatorSequence<S>) mapped to the ordinal of its NormalOperator
template <Statistics S>
OpProvenance make_provenance(const Expr &expr) {
  OpProvenance prov;
  std::size_t ord = 0;
  auto record = [&](const NormalOperator<S> &nop) {
    for (const auto &op : nop) prov.emplace(op.index(), ord);
    ++ord;
  };
  if (expr.is<NormalOperatorSequence<S>>()) {
    for (const auto &nop : expr.as<NormalOperatorSequence<S>>()) record(nop);
  } else if (expr.is<Product>()) {
    for (const auto &f : expr.as<Product>())
      if (f->template is<NormalOperator<S>>())
        record(f->template as<NormalOperator<S>>());
  } else if (expr.is<NormalOperator<S>>()) {
    record(expr.as<NormalOperator<S>>());
  }
  return prov;
}

/// renames, in every NormalOperator<S> of @p term after the first one
/// carrying it, an op index that an earlier one carries, so that each op index
/// belongs to one operator
/// @return the δs binding each new index to the one it replaces
template <Statistics S>
container::svector<ExprPtr> separate_shared_indices(Expr &term) {
  container::svector<ExprPtr> deltas;
  container::set<Index> seen;
  auto separate = [&](NormalOperator<S> &nop) {
    container::svector<Op<S>> ops;
    bool renamed = false;
    for (const auto &op : nop) {
      const Index &idx = op.index();
      if (!seen.contains(idx)) {
        ops.push_back(op);
        continue;
      }
      const auto j = Index::make_tmp_index(idx.space(), idx.proto_indices());
      ops.emplace_back(j, op.action());
      deltas.push_back(op.action() == Action::Create ? make_kronecker(j, idx)
                                                     : make_kronecker(idx, j));
      renamed = true;
    }
    for (const auto &op : nop) seen.insert(op.index());
    if (renamed) nop = make_nop<S>(ops);
  };
  if (term.is<NormalOperatorSequence<S>>()) {
    for (auto &nop : term.as<NormalOperatorSequence<S>>()) separate(nop);
  } else if (term.is<Product>()) {
    for (auto &f : term.as<Product>())
      if (f->template is<NormalOperator<S>>())
        separate(f->template as<NormalOperator<S>>());
  }
  return deltas;
}

/// @return registered spaces that partition @p type: the space of that type
/// if registered, else its base spaces
container::svector<IndexSpace> registered_pieces(
    const IndexSpaceRegistry &isr, IndexSpace::Type type,
    IndexSpace::QuantumNumbers qns) {
  container::svector<IndexSpace> result;
  if (!type) return result;
  if (const auto *sp = isr.retrieve_ptr(type, qns)) {
    result.push_back(*sp);
    return result;
  }
  for (const auto &t : isr.base_space_types())
    if (type.includes(t)) result.push_back(isr.retrieve(t, qns));
  return result;
}

/// a sum of products, each held as its list of factors
using Alternatives = container::svector<container::svector<ExprPtr>>;

/// @return the split of a 1-body γ (@p is_gamma) or η {@p bra; @p ket} into
/// a δ over its core (γ) or virtual (η) part and a γ/η over its active
/// part, or nullopt if both indices are already active
std::optional<Alternatives> split_density(const IndexSpaceRegistry &isr,
                                          const Index &bra, const Index &ket,
                                          bool is_gamma) {
  const auto parts = space_parts(isr, bra.space().qns());
  if (parts.active.includes(bra.space().type()) &&
      parts.active.includes(ket.space().type()))
    return std::nullopt;
  const auto &common = isr.intersection(bra.space(), ket.space());
  const auto inactive = is_gamma ? parts.core : parts.virt;
  SEQUANT_ASSERT(inactive.unIon(parts.active).includes(common.type()));
  Alternatives result;
  for (const auto &sp : registered_pieces(
           isr, common.type().intersection(inactive), common.qns())) {
    const auto d = Index::make_tmp_index(sp);
    result.push_back({make_kronecker(bra, d), make_kronecker(d, ket)});
  }
  if (const auto active = common.type().intersection(parts.active)) {
    const auto &sp = isr.retrieve(active, common.qns());
    const auto b = Index::make_tmp_index(sp);
    const auto k = Index::make_tmp_index(sp);
    result.push_back(
        {make_kronecker(bra, b),
         is_gamma ? density::make_rdm(b, k) : density::make_hole_rdm(b, k),
         make_kronecker(k, ket)});
  }
  return result;
}

/// @return the projections of @p nop in which every op is active or, unless
/// @p full, pure core or pure virtual; each is the projected NormalOperator
/// preceded by the δs binding projected indices to the original ones
template <Statistics S>
Alternatives split_survivors(const IndexSpaceRegistry &isr,
                             const NormalOperator<S> &nop, bool full) {
  // the projections so far: their ops and the δs they need
  container::svector<
      std::pair<container::svector<Op<S>>, container::svector<ExprPtr>>>
      partials(1);
  for (const auto &op : nop) {
    const Index &idx = op.index();
    const auto type = idx.space().type();
    const auto qns = idx.space().qns();
    const auto parts = space_parts(isr, qns);
    container::svector<IndexSpace::Type> allowed{parts.active};
    if (!full) allowed.insert(allowed.end(), {parts.core, parts.virt});
    const bool keep = std::any_of(allowed.begin(), allowed.end(),
                                  [&](auto t) { return t.includes(type); });
    container::svector<IndexSpace> targets;
    if (!keep)
      for (const auto &t : allowed)
        for (const auto &sp : registered_pieces(isr, type.intersection(t), qns))
          targets.push_back(sp);
    decltype(partials) next;
    for (const auto &[ops, deltas] : partials) {
      if (keep) {
        next.emplace_back(ops, deltas).first.push_back(op);
        continue;
      }
      for (const auto &sp : targets) {
        const auto j = Index::make_tmp_index(sp, idx.proto_indices());
        auto &[ops2, deltas2] = next.emplace_back(ops, deltas);
        ops2.emplace_back(j, op.action());
        deltas2.push_back(op.action() == Action::Create
                              ? make_kronecker(j, idx)
                              : make_kronecker(idx, j));
      }
    }
    partials = std::move(next);
  }
  Alternatives result;
  for (auto &[ops, deltas] : partials) {
    deltas.push_back(ex<NormalOperator<S>>(make_nop<S>(ops)));
    result.push_back(std::move(deltas));
  }
  return result;
}

/// rewrites @p term so that every γ and η index is active and every
/// surviving op index is active or, unless @p full, pure core or pure
/// virtual
/// @return the rewritten term as a list of Products
template <Statistics S>
container::svector<std::shared_ptr<Product>> split_mixed_spaces(
    const ExprPtr &term, const IndexSpaceRegistry &isr, bool full) {
  const auto product =
      term->is<Product>()
          ? std::static_pointer_cast<Product>(term->clone().as_shared_ptr())
          : std::make_shared<Product>(ExprPtrList{term->clone()});
  container::svector<std::shared_ptr<Product>> partials{
      std::make_shared<Product>(product->scalar(), ExprPtrList{})};
  for (const auto &f : product->factors()) {
    std::optional<Alternatives> alternatives;
    if (f->is<Tensor>()) {
      const auto &t = f->as<Tensor>();
      const bool is_gamma = t.label() == density::rdm_label();
      if ((is_gamma || t.label() == density::hole_rdm_label()) &&
          t.bra_rank() == 1 && t.ket_rank() == 1)
        alternatives = split_density(isr, t.bra()[0], t.ket()[0], is_gamma);
    } else if (f->is<NormalOperator<S>>()) {
      alternatives = split_survivors<S>(isr, f->as<NormalOperator<S>>(), full);
    }
    if (!alternatives) {
      for (auto &p : partials) p->append(1, f);
      continue;
    }
    decltype(partials) next;
    for (const auto &p : partials)
      for (const auto &alt : *alternatives) {
        auto q = std::static_pointer_cast<Product>(p->clone().as_shared_ptr());
        for (const auto &x : alt) q->append(1, x);
        next.push_back(std::move(q));
      }
    partials = std::move(next);
  }
  return partials;
}

/// adds to @p prov every index bound by a chain of Kronecker deltas of
/// @p product to an index already in @p prov
void extend_provenance(const Product &product, OpProvenance &prov) {
  bool extended;
  do {
    extended = false;
    for (const auto &f : product.factors()) {
      if (!f->is<Tensor>()) continue;
      const auto &t = f->as<Tensor>();
      if (t.label() != reserved::kronecker_label()) continue;
      const Index &b = t.bra()[0], &k = t.ket()[0];
      const auto b_it = prov.find(b), k_it = prov.find(k);
      if (b_it != prov.end() && k_it == prov.end()) {
        prov.emplace(k, b_it->second);
        extended = true;
      } else if (k_it != prov.end() && b_it == prov.end()) {
        prov.emplace(b, k_it->second);
        extended = true;
      }
    }
  } while (extended);
}

/// @return whether @p term realizes every `opts.nop_connections` pair and no
/// `opts.nop_avoided_connections` pair; a γ, η, κ, δ or overlap, i.e. a factor
/// the theorem produces, connects the input operators of all its indices;
/// other tensors, e.g. coefficients, connect nothing
bool satisfies_connectivity(const Product &term, const OpProvenance &prov,
                            const ExtendedWickOptions &opts) {
  if (opts.nop_connections.empty() && opts.nop_avoided_connections.empty())
    return true;
  using Edge = std::pair<std::size_t, std::size_t>;
  const auto edge = [](std::size_t a, std::size_t b) {
    return Edge{std::min(a, b), std::max(a, b)};
  };
  const container::set<std::wstring> connecting{
      density::rdm_label(), density::hole_rdm_label(),
      density::cumulant_label(), reserved::kronecker_label(),
      reserved::overlap_label()};
  container::set<Edge> edges;
  for (const auto &f : term.factors()) {
    if (!f->is<Tensor>() ||
        !connecting.contains(std::wstring(f->as<Tensor>().label())))
      continue;
    container::svector<std::size_t> ords;
    for (const auto &idx : f->as<Tensor>().const_braket())
      if (const auto it = prov.find(idx); it != prov.end())
        ords.push_back(it->second);
    for (std::size_t i = 0; i != ords.size(); ++i)
      for (std::size_t j = i + 1; j != ords.size(); ++j)
        if (ords[i] != ords[j]) edges.insert(edge(ords[i], ords[j]));
  }
  for (const auto &[a, b] : opts.nop_connections)
    if (!edges.contains(edge(a, b))) return false;
  for (const auto &[a, b] : opts.nop_avoided_connections)
    if (edges.contains(edge(a, b))) return false;
  return true;
}

/// @return @p expr with every δ over an index summed within its term applied
template <Statistics S>
ExprPtr apply_dummy_deltas(const ExprPtr &expr) {
  auto result = std::make_shared<Sum>();
  const auto terms = expr->is<Sum>() ? expr : ex<Sum>(ExprPtrList{expr});
  for (const auto &term : *terms) {
    if (!term->is<Product>()) {
      result->append(term);
      continue;
    }
    // the term's own index counts tell its dummies from its externals
    ExprPtr reduced = term->clone();
    WickTheorem<S> reducer{reduced};
    reducer.reduce(reduced);
    if (!reduced->is<Product>()) continue;  // vanished
    // an applied δ is left as a factor of 1; append() folds it into the scalar
    const auto &product = reduced->as<Product>();
    ExprPtr folded = std::make_shared<Product>(product.scalar(), ExprPtrList{});
    for (const auto &f : product) folded->as<Product>().append(1, f);
    // a lone tensor left by simplify would keep the reducer's dummy names
    result->append(canonicalize(folded));
  }
  return result;
}

/// rewrites every 1-body η of @p expr, a Sum, as δ - γ
void rewrite_eta(ExprPtr &expr) {
  expr->visit(
      [](ExprPtr &e) {
        if (e->is<Tensor>() &&
            e->as<Tensor>().label() == density::hole_rdm_label()) {
          const auto &t = e->as<Tensor>();
          const Index &b = t.bra()[0], &k = t.ket()[0];
          e = make_kronecker(b, k) - density::make_rdm(b, k);
        }
      },
      /*atoms_only=*/true);
  expand(expr);
}

}  // namespace

template <Statistics S>
ExprPtr cumulant_expand(const ExprPtr &wick_output,
                        const OpProvenance &provenance,
                        const ExtendedWickOptions &opts) {
  if (!wick_output || wick_output->is<Constant>())
    return wick_output ? wick_output : ex<Constant>(0);

  const std::size_t max_rank =
      opts.max_cumulant_rank.value_or(static_cast<std::size_t>(-1));
  auto result = std::make_shared<Sum>();

  auto expand_term = [&](const ExprPtr &term) {
    const auto product =
        term->is<Product>()
            ? std::static_pointer_cast<Product>(term->clone().as_shared_ptr())
            : std::make_shared<Product>(ExprPtrList{term->clone()});
    const auto &factors = product->factors();
    const auto nop_it = std::find_if(
        factors.begin(), factors.end(),
        [](const ExprPtr &f) { return f->is<NormalOperator<S>>(); });
    if (nop_it == factors.end()) {
      if (satisfies_connectivity(*product, provenance, opts))
        result->append(product);
      return;
    }
    const auto &survivors = (*nop_it)->template as<NormalOperator<S>>();

    for_each_block_assignment<S>(
        survivors, provenance, max_rank, opts.full_contractions,
        [&](int sign, const Blocks &blocks,
            const container::svector<std::size_t> &remainder) {
          auto summand = std::make_shared<Product>(
              product->scalar() * sign *
                  detail::term_weight<S>(survivors, blocks),
              ExprPtrList{});
          for (auto it = factors.begin(); it != factors.end(); ++it)
            if (it != nop_it) summand->append(1, *it);
          for (const auto &b : blocks)
            summand->append(1, detail::block_value<S>(subset(survivors, b)));
          if (!remainder.empty())
            summand->append(
                1, ex<NormalOperator<S>>(subset(survivors, remainder)));
          if (!satisfies_connectivity(*summand, provenance, opts)) return;
          result->append(summand);
        });
  };

  if (wick_output->is<Sum>()) {
    for (const auto &term : *wick_output) expand_term(term);
  } else {
    expand_term(wick_output);
  }
  ExprPtr out = result;
  simplify(out);
  if (out->is<Sum>() && out->as<Sum>().empty()) return ex<Constant>(0);
  return out;
}

template ExprPtr cumulant_expand<Statistics::FermiDirac>(
    const ExprPtr &, const OpProvenance &, const ExtendedWickOptions &);

template <Statistics S>
ExprPtr extended_wick(ExprPtr input, const ExtendedWickOptions &opts,
                      WickTheorem<S> &stats_sink) {
  const auto &ctx = get_default_context(S);
  SEQUANT_ASSERT(ctx.vacuum() == Vacuum::MultiProduct);
  const auto &isr = *ctx.index_space_registry();

  // provenance is per input term
  auto per_term = [&](ExprPtr term) -> ExprPtr {
    if (term->is<NormalOperator<S>>()) term = ex<Product>(ExprPtrList{term});
    // canonicalize first, so that the provenance sees the final indices
    if (term->is<Product>()) {
      [[maybe_unused]] auto bp = term->rapid_canonicalize();
      SEQUANT_ASSERT(bp == nullptr);
    }
    // the operators' equivalences are those of the term itself, in which an
    // index shared by two operators is a dummy; renaming it below leaves
    // every op at its ordinal
    std::optional<typename WickTheorem<S>::TopologicalPartitions> partitions;
    if (opts.use_topology && term->is<Product>())
      partitions = WickTheorem<S>::analyze_topology(term->as<Product>());
    // a summed index shared by two operators is two indices bound by a δ,
    // which multiplies the result so that it does not count as a connection
    const auto shared_deltas = separate_shared_indices<S>(*term);
    const auto provenance = make_provenance<S>(*term);

    // WickTheorem sees only the operators, so every index is external to it
    // and none is renamed; the c-number factors multiply its result. reduce
    // keeps every input index too, so every projected index stays δ-bound to
    // an input op index; apply_dummy_deltas applies the δs over true dummies
    auto nopseq = std::make_shared<NormalOperatorSequence<S>>();
    ExprPtr prefactor = ex<Constant>(1);
    container::set<Index> fixed_indices;
    for (const auto &[idx, ord] : provenance) fixed_indices.insert(idx);
    if (term->is<NormalOperatorSequence<S>>()) {
      *nopseq = term->as<NormalOperatorSequence<S>>();
    } else if (term->is<Product>()) {
      prefactor = ex<Constant>(term->as<Product>().scalar());
      for (const auto &f : term->as<Product>()) {
        if (f->is<NormalOperator<S>>()) {
          nopseq->push_back(f->as<NormalOperator<S>>());
          continue;
        }
        prefactor = prefactor * f;
        if (f->is<Tensor>())
          for (const auto &idx : f->as<Tensor>().const_braket())
            fixed_indices.insert(idx);
      }
    }
    if (nopseq->empty()) return term;
    for (const auto &pairs :
         {opts.nop_connections, opts.nop_avoided_connections})
      for (const auto &[a, b] : pairs)
        if (std::max(a, b) >= nopseq->size())
          throw Exception("WickTheorem::compute: connection ordinal " +
                          std::to_string(std::max(a, b)) + " exceeds the " +
                          std::to_string(nopseq->size()) +
                          " input operators of a term");

    WickTheorem<S> wick{nopseq};
    wick.full_contractions(false).use_topology(opts.use_topology);
    if (partitions && !partitions->nop_partitions.empty())
      wick.set_nop_partitions(partitions->nop_partitions);
    if (partitions && !partitions->op_partitions.empty())
      wick.set_op_partitions(partitions->op_partitions);
    const ExprPtr raw = wick.compute_contractions(
        /*count_only=*/false, /*skip_input_canonicalization=*/false);
    stats_sink.stats() += wick.stats();

    auto result = std::make_shared<Sum>();
    auto process = [&](const ExprPtr &t) {
      for (const auto &p :
           split_mixed_spaces<S>(ex<Product>(ExprPtrList{prefactor, t}), isr,
                                 opts.full_contractions)) {
        ExprPtr reduced = p;
        // the operator and coefficient indices are external to this reduction
        const auto scope =
            set_scoped_modified_default_context([&fixed_indices](Context &ctx) {
              ctx.set(ctx.canonicalization_options()
                          .value_or(CanonicalizeOptions{})
                          .copy_and_set(fixed_indices));
            });
        WickTheorem<S> reducer{reduced};
        reducer.reduce(reduced);
        if (!reduced->is<Product>()) continue;  // vanished
        OpProvenance prov = provenance;
        extend_provenance(reduced->as<Product>(), prov);
        ExprPtr expanded = cumulant_expand<S>(reduced, prov, opts);
        for (const auto &d : shared_deltas) expanded = expanded * d;
        expand(expanded);
        result->append(expanded);
      }
    };
    if (raw->is<Sum>()) {
      for (const auto &t : *raw) process(t);
    } else if (raw->is<Constant>()) {
      ExprPtr value = prefactor * raw;
      for (const auto &d : shared_deltas) value = value * d;
      result->append(value);
    } else {
      process(raw);
    }
    return result;
  };

  input = input->clone();
  expand(input);
  auto result = std::make_shared<Sum>();
  if (input->is<Sum>()) {
    for (const auto &t : *input) result->append(per_term(t));
  } else {
    result->append(per_term(input));
  }
  ExprPtr out = result;
  if (opts.eta_as_delta_minus_gamma) rewrite_eta(out);
  out = apply_dummy_deltas<S>(out);
  simplify(out);
  if (out->is<Sum>() && out->as<Sum>().empty()) return ex<Constant>(0);
  return out;
}

template ExprPtr extended_wick<Statistics::FermiDirac>(
    ExprPtr, const ExtendedWickOptions &,
    WickTheorem<Statistics::FermiDirac> &);

}  // namespace sequant::detail
