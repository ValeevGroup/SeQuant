#include <SeQuant/core/wick_extended.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <functional>
#include <memory>

namespace sequant {

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

/// @return the ops of @p nop at storage positions @p positions as a
/// MultiProduct-vacuum NormalOperator
template <Statistics S>
NormalOperator<S> subset(const NormalOperator<S> &nop,
                         const container::svector<std::size_t> &positions) {
  container::svector<Op<S>> cre_ops, ann_ops;
  for (auto i : positions)
    (nop[i].action() == Action::Create ? cre_ops : ann_ops).push_back(nop[i]);
  // the ctor takes annihilators in particle order, the reverse of storage
  std::reverse(ann_ops.begin(), ann_ops.end());
  return NormalOperator<S>(cre(std::move(cre_ops)), ann(std::move(ann_ops)),
                           Vacuum::MultiProduct);
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

}  // namespace sequant
