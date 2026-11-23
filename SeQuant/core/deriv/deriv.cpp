#include <SeQuant/core/deriv/deriv.hpp>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/tensor_power.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <ranges>
#include <string>
#include <tuple>
#include <utility>

namespace sequant::deriv {

namespace {

/// @return the derivative of the factor @p f by the scalar variable @p y, or
/// nullptr if @p f is a constant with respect to @p y
ExprPtr scalar_rule(const ExprPtr& f, const ExprPtr& y) {
  if (f->is<Variable>()) return *f == *y ? ex<Constant>(1) : nullptr;
  if (!f->is<Power>()) return nullptr;
  const auto& power = f->as<Power>();
  auto base = power.base()->clone();
  if (power.conjugated()) base->adjoint();
  if (*base != *y) return nullptr;
  ExprPtr lowered = ex<Power>(power.base(), power.exponent() - 1);
  if (power.conjugated()) lowered->as<Power>().conjugate();
  Power::flatten(lowered);
  return ex<Product>(power.exponent(), ExprPtrList{lowered});
}

/// @return whether @p f is a symmetrizer (Â or Ŝ)
bool is_symmetrizer(const ExprPtr& f) {
  if (!f->is<Tensor>()) return false;
  const auto label = f->as<Tensor>().label();
  return label == reserved::antisymm_label() || label == reserved::symm_label();
}

/// @return the scalar and the factors of @p term , a term of an expanded
/// expression
std::pair<Constant::scalar_type, container::svector<ExprPtr>> split(
    const ExprPtr& term) {
  Constant::scalar_type scalar = 1;
  container::svector<ExprPtr> factors;
  const auto walk = [&](const auto& self, const ExprPtr& e) -> void {
    if (e->is<Product>()) {
      scalar *= e->as<Product>().scalar();
      for (const auto& f : e->as<Product>().factors()) self(self, f);
    } else if (e->is<Constant>()) {
      scalar *= e->as<Constant>().value();
    } else {
      SEQUANT_ASSERT(!e->is<Sum>());
      factors.push_back(e);
    }
  };
  walk(walk, term);
  return {scalar, std::move(factors)};
}

/// @return the terms of the expanded @p E
container::svector<ExprPtr> terms_of(const ExprPtr& E) {
  if (!E->is<Sum>()) return {E};
  const auto& summands = E->as<Sum>().summands();
  return {summands.begin(), summands.end()};
}

/// a spelling of a tensor variable, and its phase
using Spelling = std::pair<ExprPtr, Constant::scalar_type>;

/// core's spelling of the tensor @p x and its phase
Spelling spelled(const ExprPtr& x) {
  auto result = x->clone();
  const auto phase =
      TensorBlockCanonicalizer{}.apply(result->as<AbstractTensor>());
  return {result,
          phase ? phase->as<Constant>().value() : Constant::scalar_type{1}};
}

/// the external indices of a term of the function of diff(): those in one
/// slot of the factors other than Â and Ŝ, those of @p known that the term
/// holds, and the protoindices of those. An index in more slots that is not
/// known is summed: in two, a contraction; in more, a Hadamard product
container::set<Index> externals_of(const container::svector<ExprPtr>& factors,
                                   const container::set<Index>& known) {
  container::map<Index, std::int64_t> slots;
  // the indices of a symmetrizer are those of a residual: external, whatever
  // the slots of the other factors that hold them
  container::set<Index> symmetrized;
  for (const auto& f : factors) {
    if (is_symmetrizer(f)) {
      for (const auto& i : get_used_indices(f)) symmetrized.insert(i);
      continue;
    }
    for (const auto& [i, count] : get_used_indices_with_counts(*f))
      slots[i] += count.nonproto();
  }
  container::set<Index> result;
  for (const auto& [i, n] : slots) {
    if (n >= 2 && !known.contains(i) && !symmetrized.contains(i)) continue;
    result.insert(i);
    for (const auto& p : i.proto_indices()) result.insert(p);
  }
  return result;
}

/// a term of an expanded expression: its scalar, its factors and its external
/// indices
using Term = std::tuple<Constant::scalar_type, container::svector<ExprPtr>,
                        container::set<Index>>;

/// the terms of the expanded @p E ; their external indices (externals_of(),
/// with the indices @p known to be external) are found only if
/// @p with_externals : the scalar rule needs none
container::svector<Term> split_terms(const ExprPtr& E, bool with_externals,
                                     const container::set<Index>& known = {}) {
  container::svector<Term> result;
  for (const auto& term : terms_of(E)) {
    auto [scalar, factors] = split(term);
    auto externals =
        with_externals ? externals_of(factors, known) : container::set<Index>{};
    result.emplace_back(scalar, std::move(factors), std::move(externals));
  }
  return result;
}

/// @return a temporary index of the basis (and protoindices) of @p i whose
/// label @p unused accepts
template <typename Validator>
Index temporary(const Index& i, const Validator& unused) {
  Index result;
  do
    result = i.has_proto_indices()
                 ? Index::make_tmp_index(i.basis(), i.proto_indices(),
                                         i.symmetric_proto_indices())
                 : Index::make_tmp_index(i.basis());
  while (!unused(result));
  return result;
}

/// renames the summed indices of each of @p terms that are named like an
/// index of @p x to temporary ones named unlike any index of x and of the term
void rename_apart(container::svector<Term>& terms, const ExprPtr& x) {
  container::set<std::wstring> in_x;
  for (const auto& i : get_used_indices(x)) in_x.emplace(i.label());
  for (auto& [scalar, factors, externals] : terms) {
    container::set<Index> used;
    for (const auto& f : factors) {
      const auto in_f = get_used_indices(f);
      used.insert(in_f.begin(), in_f.end());
    }
    auto taken = in_x;
    for (const auto& i : used) taken.emplace(i.label());
    const auto unused = unused_label_validator(std::move(taken));
    container::map<Index, Index> fresh;
    for (const auto& i : used)
      if (!externals.contains(i) && in_x.contains(std::wstring(i.label())))
        fresh.emplace(i, temporary(i, unused));
    if (fresh.empty()) continue;
    const auto renames = replace_proto_indices_of_images(fresh);
    for (auto& f : factors) f = transform_expr(f, renames);
  }
}

/// the position of @p p among the slots of @p t , or -1 if it is not a slot
std::ptrdiff_t slot_of(const Tensor& t, const Index& p) {
  const auto slots = t.slots();
  const auto it = std::ranges::find(slots, p);
  return it == slots.end() ? -1 : it - slots.begin();
}

/// the position of the first slot of @p t that carries the protoindices
/// @p list
std::size_t first_with(const Tensor& t, const Index::index_vector& list) {
  const auto slots = t.slots();
  return std::ranges::find(slots, list, &Index::proto_indices) - slots.begin();
}

/// ∂ leaf / ∂ y for a leaf of a term, as what it does to the rest of the term
struct Elementary {
  /// summed labels of the leaf -> labels of y
  container::map<Index, Index> renames;
  /// make_kronecker factors that stay
  container::svector<ExprPtr> deltas;
  /// whether a protoindex of the leaf that is not a slot is not renamed
  bool unrenamed = false;
};

/// whether the spaces of @p l and @p r overlap
bool overlap(const Index& l, const Index& r) {
  const auto registry = get_default_index_basis_registry();
  SEQUANT_ASSERT(registry);
  return registry->valid_intersection(l.space(), r.space());
}

/// whether the slot @p l of a leaf can hold the slot @p y_l of the variable: a
/// basis that is no other instance of y_l's, a space that overlaps y_l's and
/// as many protoindices
bool meets(const Index& l, const Index& y_l) {
  return !different_instances(l.basis(), y_l.basis()) && overlap(l, y_l) &&
         l.proto_indices().size() == y_l.proto_indices().size();
}

/// whether @p f , as it is written, is a leaf of the block of the spelling
/// @p y : the label and ranks of y, and in each slot a slot that meets() y's
bool leaf_of(const ExprPtr& f, const Tensor& y) {
  if (!f->is<Tensor>()) return false;
  const auto& t = f->as<Tensor>();
  return t.label() == y.label() && t.bra_rank() == y.bra_rank() &&
         t.ket_rank() == y.ket_rank() && t.aux_rank() == y.aux_rank() &&
         std::ranges::equal(t.slots(), y.slots(), meets);
}

/// the derivative of @p leaf by @p y : the Kronecker delta of each label of
/// @p leaf (slot index, protoindex) with that of @p y , contracted away where
/// the label is summed (not in @p externals ) and runs over all of y's
/// @return one Elementary per choice of an order of each distinct protoindex
/// list of y that is not its slots (one, unless such a list has two or more
/// labels), but none that has a Kronecker delta between disjoint spaces, which
/// is zero. @p leaf is a leaf of y's block (leaf_of())
container::svector<Elementary> elementary(
    const Tensor& leaf, const Tensor& y,
    const container::set<Index>& externals) {
  // (leaf label, y label, is it a bra slot)
  container::svector<std::tuple<Index, Index, bool>> pairs;
  const auto add_slots = [&](const auto& ls, const auto& ys, bool bra) {
    for (std::size_t k = 0; k != ls.size(); ++k)
      pairs.emplace_back(ls[k], ys[k], bra);
  };
  add_slots(leaf.bra(), y.bra(), true);
  add_slots(leaf.ket(), y.ket(), false);
  add_slots(leaf.aux(), y.aux(), false);

  // per distinct protoindex list of y that is not y's slots: the leaf's list
  // at its first slot, y's list, and the range of the list in both
  container::svector<Index> free, y_free;
  container::svector<std::pair<std::size_t, std::size_t>> lists;
  for (std::size_t k = 0; k != pairs.size(); ++k) {
    const auto& [l, y_l, bra] = pairs[k];
    const auto& Q = y_l.proto_indices();
    if (Q.empty() || first_with(y, Q) != k || slot_of(y, Q.front()) >= 0)
      continue;
    lists.emplace_back(free.size(), free.size() + Q.size());
    free.insert(free.end(), l.proto_indices().begin(), l.proto_indices().end());
    y_free.insert(y_free.end(), Q.begin(), Q.end());
  }

  container::svector<Elementary> result;
  container::svector<std::size_t> order(y_free.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  // the next order of the first list that has one; the lists before it start
  // over
  const auto next = [&order](const auto& list) {
    return std::ranges::next_permutation(order.begin() + list.first,
                                         order.begin() + list.second)
        .found;
  };
  do {
    Elementary r;
    bool zero = false;
    // every rename first, so that each delta is spelled under all of them
    const auto rename = [&](const Index& l, const Index& y_l) {
      if (l != y_l && !externals.contains(l) &&
          includes(l.basis(), y_l.basis()) && !r.renames.contains(l))
        r.renames.emplace(l, y_l);
    };
    const auto delta = [&](const Index& l, const Index& y_l, bool bra) {
      auto renamed = l;
      renamed.transform(r.renames);
      if (renamed == y_l) return;
      zero = zero || !overlap(renamed, y_l);
      r.deltas.push_back(bra ? make_kronecker(renamed, y_l)
                             : make_kronecker(y_l, renamed));
    };
    for (const auto& [l, y_l, bra] : pairs) rename(l, y_l);
    for (std::size_t k = 0; k != free.size(); ++k)
      rename(free[k], y_free[order[k]]);
    for (const auto& l : free) r.unrenamed |= !r.renames.contains(l);
    for (const auto& [l, y_l, bra] : pairs) delta(l, y_l, bra);
    for (std::size_t k = 0; k != free.size(); ++k)
      delta(free[k], y_free[order[k]], true);
    if (!zero) result.push_back(std::move(r));
  } while (std::ranges::any_of(lists, next));
  return result;
}

/// the definitions of tensor leaves that diff() differentiates through
using Definitions = container::svector<ResultExpr>;

/// @return the definition in @p defs of a leaf @p t : the one whose head has
/// t's label and ranks, or nullptr if there is none
const ResultExpr* definition_of(const Tensor& t, const Definitions& defs) {
  for (const auto& d : defs)
    if (d.label() == t.label() && d.bra().size() == t.bra_rank() &&
        d.ket().size() == t.ket_rank() && d.aux().size() == t.aux_rank())
      return &d;
  return nullptr;
}

/// @return the expression of the definition @p def of the leaf @p leaf ,
/// with the slots of the head renamed to those of the leaf, slot by slot, and
/// the summed indices to temporary ones whose labels are not in @p taken
ExprPtr instantiate(const ResultExpr& def, const Tensor& leaf,
                    container::set<std::wstring> taken) {
  container::map<Index, Index> renames;
  const auto map_slots = [&renames](const auto& from, const auto& to) {
    for (std::size_t k = 0; k != from.size(); ++k)
      renames.emplace(from[k], to[k]);
  };
  map_slots(def.bra(), leaf.bra());
  map_slots(def.ket(), leaf.ket());
  map_slots(def.aux(), leaf.aux());
  const auto used = get_used_indices(def.expression());
  for (const auto& i : used) taken.emplace(i.label());
  for (const auto& i : leaf.slots()) taken.emplace(i.label());
  const auto unused = unused_label_validator(std::move(taken));
  for (const auto& i : used)
    if (!renames.contains(i)) renames.emplace(i, temporary(i, unused));
  auto result = transform_expr(def.expression(),
                               replace_proto_indices_of_images(renames));
  expand(result);
  return result;
}

/// the labels of the indices of the @p factors of a term and of @p y
container::set<std::wstring> labels_of(
    const container::svector<ExprPtr>& factors, const ExprPtr& y) {
  container::set<std::wstring> result;
  for (const auto& i : get_used_indices(y)) result.emplace(i.label());
  for (const auto& f : factors)
    for (const auto& i : get_used_indices(f)) result.emplace(i.label());
  return result;
}

/// D_y(E): the product rule over the factors of each of the @p terms of E,
/// with elementary() for a leaf of y's label, the scalar rule for a scalar
/// @p y , the power rule for a TensorPower of a leaf of y's label, and the
/// chain rule through the definitions @p defs for a leaf (or a TensorPower of
/// one) that has one
ExprPtr one_way(const container::svector<Term>& terms, const ExprPtr& y,
                const Definitions& defs) {
  container::vector<ExprPtr> result;
  for (const auto& [scalar, factors, externals] : terms) {
    const auto product = [&](std::size_t skip,
                             const container::map<Index, Index>& renames) {
      auto p = ex<Product>(scalar, ExprPtrList{});
      for (std::size_t j = 0; j != factors.size(); ++j)
        if (j != skip)
          p->as<Product>().append(1, renames.empty()
                                         ? factors[j]
                                         : transform_expr(factors[j], renames));
      return p;
    };
    for (std::size_t k = 0; k != factors.size(); ++k) {
      const auto& f = factors[k];
      if (y->is<Variable>()) {
        if (auto d = scalar_rule(f, y)) {
          auto p = product(k, {});
          p->as<Product>().append(1, d);
          result.push_back(p);
          continue;
        }
      }
      const bool power = f->is<TensorPower>();
      const auto& base = power ? f->as<TensorPower>().base() : f;
      if (!base->is<Tensor>()) continue;
      const auto& leaf = base->as<Tensor>();
      // the power rule: e base^(e-1), its indices renamed by @p renames
      const auto lowered =
          [&](const container::map<Index, Index>& renames) -> ExprPtr {
        const auto e = f->as<TensorPower>().exponent();
        const auto b = renames.empty() ? base : transform_expr(base, renames);
        return e == 1 ? ex<Constant>(1)
                      : ex<Constant>(e) * pow(b->clone(), e - 1);
      };
      if (y->is<Tensor>() && leaf_of(base, y->as<Tensor>())) {
        for (const auto& r : elementary(leaf, y->as<Tensor>(), externals)) {
          SEQUANT_ASSERT(!r.unrenamed);
          auto p = product(k, r.renames);
          if (power) p->as<Product>().append(1, lowered(r.renames));
          for (const auto& delta : r.deltas) p->as<Product>().append(1, delta);
          result.push_back(p);
        }
        continue;
      }
      const auto* def = definition_of(leaf, defs);
      if (!def) continue;
      if (std::ranges::any_of(leaf.slots(), &Index::has_proto_indices))
        throw Exception(
            "deriv::diff: a leaf with a definition has protoindices");
      container::set<Index> leaf_slots;
      for (const auto& i : leaf.slots()) leaf_slots.insert(i);
      const auto inner_terms =
          split_terms(instantiate(*def, leaf, labels_of(factors, y)),
                      /*with_externals=*/y->is<Tensor>(), leaf_slots);
      auto inner = one_way(inner_terms, y, defs);
      expand(inner);
      // per term of the derivative of the definition, a Kronecker delta of a
      // summed index of this term with an index of y is summed away, as
      // elementary() does for a leaf of y's block
      const auto in_y = get_used_indices(y);
      for (const auto& term : terms_of(inner)) {
        if (term->is_zero()) continue;
        auto [inner_scalar, inner_factors] = split(term);
        container::map<Index, Index> renames;
        container::svector<ExprPtr> rest;
        for (const auto& g : inner_factors) {
          if (g->is<Tensor>() &&
              g->as<Tensor>().label() == reserved::kronecker_label()) {
            const auto& delta = g->as<Tensor>();
            const auto& u = delta.bra()[0];
            const auto& v = delta.ket()[0];
            const auto summed = [&](const Index& l) {
              return !externals.contains(l) && !l.has_proto_indices() &&
                     !renames.contains(l) && slot_of(leaf, l) >= 0;
            };
            if (summed(u) && in_y.contains(v) &&
                includes(u.basis(), v.basis())) {
              renames.emplace(u, v);
              continue;
            }
            if (summed(v) && in_y.contains(u) &&
                includes(v.basis(), u.basis())) {
              renames.emplace(v, u);
              continue;
            }
          }
          rest.push_back(g);
        }
        auto p = product(k, renames);
        p->as<Product>().scale(inner_scalar);
        if (power) p->as<Product>().append(1, lowered(renames));
        for (const auto& g : rest)
          p->as<Product>().append(
              1, renames.empty() ? g : transform_expr(g, renames));
        result.push_back(p);
      }
    }
  }
  return result.empty() ? ex<Constant>(0)
                        : ex<Sum>(result.begin(), result.end());
}

/// @return whether the adjoint spelling of @p y holds the elements of @p y
bool adjoint_holds_elements(const Tensor& y) {
  return y.braket_symmetry() == BraKetSymmetry::Symm ||
         base_field(y) == Field::Real;
}

/// whether the slots @p l and @p r of the variable @p y are like: one basis,
/// the same protoindices, and each protoindex list of y holds both or neither
bool alike(const Tensor& y, const Index& l, const Index& r) {
  return l.basis() == r.basis() && l.proto_indices() == r.proto_indices() &&
         std::ranges::all_of(y.slots(), [&](const Index& i) {
           return std::ranges::count(i.proto_indices(), l) ==
                  std::ranges::count(i.proto_indices(), r);
         });
}

/// whether @p y is column-symmetric without slot symmetry
bool column_symmetric(const Tensor& y) {
  return y.symmetry() == Symmetry::Nonsymm &&
         y.column_symmetry() == ColumnSymmetry::Symm;
}

/// the pairs of positions among the slots of @p x of two slots of one side
/// that a symmetry of x relates: any two (Antisymm); those of two columns
/// (column symmetry)
container::svector<std::pair<std::size_t, std::size_t>> related_slots(
    const Tensor& x) {
  const bool antisymm = x.symmetry() == Symmetry::Antisymm;
  const bool columns = column_symmetric(x);
  const auto nb = x.bra_rank();
  const auto nk = x.ket_rank();
  // whether slot l of a side, paired with an earlier one, is related to it
  const auto related = [&](std::size_t l) {
    return antisymm || (columns && l < std::min(nb, nk));
  };
  container::svector<std::pair<std::size_t, std::size_t>> result;
  for (std::size_t k = 0; k != nb; ++k)
    for (std::size_t l = k + 1; l != nb; ++l)
      if (related(l)) result.emplace_back(k, l);
  for (std::size_t k = 0; k != nk; ++k)
    for (std::size_t l = k + 1; l != nk; ++l)
      if (related(l)) result.emplace_back(nb + k, nb + l);
  return result;
}

/// the slots of @p y in its bra and in its ket that its symmetry permutes
/// among like ones: those of a related pair that are alike, and for column
/// symmetry whose column partners are alike too
std::pair<container::svector<Index>, container::svector<Index>> like_slots(
    const Tensor& y) {
  const auto slots = y.slots();
  const auto nb = y.bra_rank();
  const auto partner = [nb](std::size_t k) { return k < nb ? nb + k : k - nb; };
  container::svector<bool> like(y.num_slots(), false);
  for (const auto& [k, l] : related_slots(y))
    if (alike(y, slots[k], slots[l]) &&
        (!column_symmetric(y) ||
         alike(y, slots[partner(k)], slots[partner(l)])))
      like[k] = like[l] = true;
  container::svector<Index> B, K;
  for (std::size_t k = 0; k != like.size(); ++k)
    if (like[k]) (k < nb ? B : K).push_back(slots[k]);
  return {B, K};
}

/// P_y: @p G summed over the permutations of the like slots (or columns) of
/// @p y , i.e. each term of @p G times SeQuant's normalized Â (or Ŝ) over them
/// and the number of permutations it covers
ExprPtr block_sum(const ExprPtr& G, const Tensor& y) {
  const auto [B, K] = like_slots(y);
  if (B.empty() && K.empty()) return G;
  const bool columns = column_symmetric(y);
  const auto count =
      columns ? factorial(B.size()) : factorial(B.size()) * factorial(K.size());
  const auto P = columns ? make_symmetrizer(bra(B), ket(K))
                         : make_antisymmetrizer(bra(B), ket(K));
  container::vector<ExprPtr> result;
  for (const auto& term : terms_of(G))
    if (!term->is_zero())
      result.push_back(ex<Product>(count, ExprPtrList{P->clone(), term}));
  return result.empty() ? ex<Constant>(0)
                        : ex<Sum>(result.begin(), result.end());
}

/// per slot of @p t : the positions among the slots of t of its protoindices,
/// sorted, -1 for one that is not a slot; and the position of the first slot
/// that carries the same protoindex list
auto proto_pattern(const Tensor& t) {
  container::svector<container::svector<std::ptrdiff_t>> positions;
  container::svector<std::size_t> first;
  for (const auto& i : t.slots()) {
    auto& at = positions.emplace_back();
    for (const auto& p : i.proto_indices()) at.push_back(slot_of(t, p));
    std::ranges::sort(at);
    first.push_back(first_with(t, i.proto_indices()));
  }
  return std::pair{std::move(positions), std::move(first)};
}

/// throws, with the message starting with @p where , for a variable @p x that
/// is rejected whatever the function
void check_variable(const ExprPtr& x, const std::string& where) {
  const auto reject = [&where](const char* text) {
    throw Exception(where + text);
  };
  if (!x->is<Variable>() &&
      (!x->is<Tensor>() || is_symmetrizer(x) ||
       x->as<Tensor>().label() == reserved::kronecker_label()))
    reject("the variable is not a Variable or a tensor leaf");
  if (!x->is<Tensor>()) return;
  const auto& t = x->as<Tensor>();
  if (std::ranges::any_of(get_used_indices_with_counts(x), [](const auto& c) {
        return c.second.nonproto() > 1;
      }))
    reject("the variable repeats an index");
  for (const auto& positions : proto_pattern(t).first) {
    const auto others = std::ranges::count(positions, -1);
    if (others != 0 && others != std::ssize(positions))
      reject(
          "an index of the variable mixes protoindices that are its slots "
          "with others");
  }
  if (t.symmetry() == Symmetry::Symm)
    reject("a variable with Symmetry::Symm is not supported");
}

/// throws, with the message starting with @p where , for a @p leaf of the
/// label and slot counts of @p x (a variable or its adjoint) that declares
/// other symmetries: one whose slots each overlap a slot of x, and vice versa,
/// may hold elements of x
void check_symmetries(const Tensor& leaf, const Tensor& x,
                      const std::string& where) {
  // whether each slot of l overlaps some slot of r
  const auto covers = [](const Tensor& l, const Tensor& r) {
    return std::ranges::all_of(l.slots(), [&r](const Index& i) {
      return std::ranges::any_of(
          r.slots(), [&i](const Index& j) { return overlap(i, j); });
    });
  };
  if (leaf.label() == x.label() && leaf.aux_rank() == x.aux_rank() &&
      std::minmax(leaf.bra_rank(), leaf.ket_rank()) ==
          std::minmax(x.bra_rank(), x.ket_rank()) &&
      covers(leaf, x) && covers(x, leaf) &&
      (leaf.symmetry() != x.symmetry() ||
       leaf.braket_symmetry() != x.braket_symmetry() ||
       leaf.column_symmetry() != x.column_symmetry()))
    throw Exception(where +
                    "a leaf declares other symmetries than the variable");
}

/// throws, with the message starting with @p where , for a @p leaf of the
/// block of the spelling @p y (leaf_of()) whose protoindices do not follow its
/// slots as y's do: the same slots as protoindices, and one list where y has
/// one
void check_layout(const Tensor& leaf, const Tensor& y,
                  const std::string& where) {
  const auto [positions, first] = proto_pattern(leaf);
  const auto [y_positions, y_first] = proto_pattern(y);
  if (positions != y_positions ||
      std::ranges::any_of(
          std::views::iota(std::size_t{0}, first.size()),
          [&](std::size_t k) { return first[k] != first[y_first[k]]; }))
    throw Exception(where +
                    "a leaf's protoindices do not follow its slots as the "
                    "variable's do");
}

/// whether the bra-ket exchange is an operation of @p y : y has bra-ket
/// symmetry, its adjoint spelling holds its elements, and its bra and its ket
/// have one nonzero rank
bool exchanges(const Tensor& y) {
  return y.braket_symmetry() != BraKetSymmetry::Nonsymm &&
         adjoint_holds_elements(y) && y.bra_rank() == y.ket_rank() &&
         y.bra_rank() != 0;
}

/// the spellings of the elements of @p y that the operations of y make, each
/// with its phase, one per set of operations that block_sum() sums as one
/// factor, @p y first. The operations are the permutations of the slots of
/// each side (Symmetry::Antisymm, with their signs) or of the columns (column
/// symmetry, equal ranks), each also followed by the exchange of bra and ket
/// if that is an operation (exchanges())
container::svector<Spelling> spellings(const Spelling& y) {
  const auto& t = y.first->as<Tensor>();
  const container::svector<Index> B(t.bra().begin(), t.bra().end()),
      K(t.ket().begin(), t.ket().end());
  const bool antisymm = t.symmetry() == Symmetry::Antisymm;
  const bool columns = column_symmetric(t) && B.size() == K.size();
  const bool exchange = exchanges(t);
  // a spelling up to permutations of like slots: a like slot stands for its
  // group
  const auto [like_B, like_K] = like_slots(t);
  const auto key = [&](const Tensor& s) {
    container::svector<Index> result;
    for (const auto& i : s.slots())
      result.push_back(std::ranges::count(like_B, i)   ? like_B.front()
                       : std::ranges::count(like_K, i) ? like_K.front()
                                                       : i);
    return result;
  };
  const auto permuted = [](const auto& slots, const auto& order) {
    container::svector<Index> result;
    for (const auto k : order) result.push_back(slots[k]);
    return result;
  };

  container::svector<Spelling> result;
  container::svector<container::svector<Index>> keys;
  container::svector<std::size_t> order_B(B.size()), order_K(K.size());
  std::iota(order_B.begin(), order_B.end(), std::size_t{0});
  std::iota(order_K.begin(), order_K.end(), std::size_t{0});
  int parity_B = 0, parity_K = 0;
  do {
    do {
      const auto b = permuted(B, order_B);
      const auto k = permuted(K, columns ? order_B : order_K);
      const auto phase =
          antisymm && parity_B != parity_K ? -y.second : y.second;
      for (const bool exchanged : {false, true}) {
        if (exchanged && !exchange) break;
        auto s = y.first->clone();
        auto& s_t = s->as<Tensor>();
        s_t.set_bra(
            {exchanged ? k.begin() : b.begin(), exchanged ? k.end() : b.end()});
        s_t.set_ket(
            {exchanged ? b.begin() : k.begin(), exchanged ? b.end() : k.end()});
        auto s_key = key(s_t);
        if (std::ranges::find(keys, s_key) != keys.end()) continue;
        keys.push_back(std::move(s_key));
        result.emplace_back(std::move(s), phase);
      }
    } while (antisymm &&
             next_permutation_parity(parity_K, order_K.begin(), order_K.end()));
  } while ((antisymm || columns) &&
           next_permutation_parity(parity_B, order_B.begin(), order_B.end()));
  return result;
}

/// the spellings() of core's spelling y of the tensor variable @p x , then,
/// if @p with_adjoint , those of the adjoint of y, with y's phase, where that
/// is another spelling. The adjoint is not respelled: core may spell it as y
/// where the bra-ket exchange is no operation of x
container::svector<Spelling> all_spellings(const ExprPtr& x,
                                           bool with_adjoint) {
  const auto y = spelled(x);
  auto result = spellings(y);
  if (with_adjoint) {
    const Spelling y_adjoint{adjoint(y.first), y.second};
    if (*y_adjoint.first != *y.first)
      for (auto& s : spellings(y_adjoint)) result.push_back(std::move(s));
  }
  return result;
}

/// throws, with the message starting with @p where , for an input whose
/// derivative is not defined: the tensor @p x , accepted by check_variable(),
/// its spellings @p ys (by_one()), and the @p terms of the function, renamed
/// apart from x
void check(const container::svector<Term>& terms, const ExprPtr& x,
           const container::svector<Spelling>& ys, const std::string& where) {
  const auto reject = [&where](const char* text) {
    throw Exception(where + text);
  };
  const auto counts = get_used_indices_with_counts(x);
  const auto& t = x->as<Tensor>();
  for (const auto& i : t.slots())
    if (i.has_proto_indices() && !i.symmetric_proto_indices())
      reject("ordered protoindices are not supported");
  // a protoindex that is not a slot stands in just the slots with its list
  const auto slots = t.slots();
  const auto first_slot = proto_pattern(t).second;
  for (std::size_t k = 0; k != first_slot.size(); ++k)
    for (const auto& p : slots[k].proto_indices())
      if (slot_of(t, p) < 0 &&
          counts.at(p).proto != std::ranges::count(first_slot, first_slot[k]))
        reject(
            "the protoindex lists of the variable's slots share some but not "
            "all labels");
  const auto [B, K] = like_slots(t);
  const auto one_group = [&t](const auto& like) {
    return std::ranges::all_of(
        like, [&](const Index& i) { return alike(t, i, like.front()); });
  };
  if (!one_group(B) || !one_group(K))
    reject("the variable has two groups of like slots on one side");

  const auto x_adjoint = adjoint(x);
  container::set<std::wstring> names;
  for (const auto& [i, count] : counts) names.emplace(i.label());
  for (const auto& [scalar, factors, externals] : terms) {
    for (const auto& i : externals)
      if (names.contains(std::wstring(i.label())))
        reject("the variable names an external index of the function");
    for (const auto& f : factors) {
      const bool power = f->is<TensorPower>();
      const auto& base = power ? f->as<TensorPower>().base() : f;
      if (!base->is<Tensor>()) continue;
      const auto& b = base->as<Tensor>();
      for (const auto& r : {x, x_adjoint})
        check_symmetries(b, r->as<Tensor>(), where);
      for (const auto& y : ys) {
        const auto& y_t = y.first->as<Tensor>();
        if (!leaf_of(base, y_t)) continue;
        check_layout(b, y_t, where);
        if (std::ranges::any_of(elementary(b, y_t, externals),
                                &Elementary::unrenamed))
          reject(
              "a protoindex label of a leaf of the variable must be summed "
              "over the space of the variable's label");
      }
    }
  }

  // with at most one slot per side the column symmetry permutes nothing
  if (column_symmetric(t) && t.bra_rank() != t.ket_rank() &&
      std::max(t.bra_rank(), t.ket_rank()) > 1)
    reject("a column-symmetric variable with unequal ranks is not supported");
}

/// throws, with the message starting with @p where , for definitions @p defs
/// that diff() cannot differentiate through by the variable @p x : a head
/// that is not a labeled tensor, has protoindices, has the label of @p x or
/// of another head with its ranks; an expression that holds operators; or a
/// definition that depends on itself, directly or through others
void check_definitions(const Definitions& defs, const ExprPtr& x,
                       const std::string& where) {
  const auto reject = [&where](const char* text) {
    throw Exception(where + text);
  };
  // the leaves of each definition's expression that have a definition
  container::svector<container::svector<std::size_t>> uses(defs.size());
  for (std::size_t k = 0; k != defs.size(); ++k) {
    const auto& d = defs[k];
    if (!d.has_label() || !d.produces_tensor())
      reject("a definition's head is not a labeled tensor");
    for (const auto& side : {d.bra(), d.ket(), d.aux()})
      if (std::ranges::any_of(side, &Index::has_proto_indices))
        reject("a definition's head has protoindices");
    if (!d.expression() || !d.expression()->is_cnumber())
      reject("a definition holds operators");
    const auto head = d.result_as_tensor();
    container::set<Index> head_slots;
    for (const auto& i : head.slots())
      if (!head_slots.insert(i).second)
        reject("a definition's head repeats an index");
    if (x->is<Tensor>() && x->as<Tensor>().label() == head.label())
      reject("a definition's head has the label of the variable");
    if (definition_of(head, Definitions(defs.begin(), defs.begin() + k)))
      reject("two definitions have one head");
    auto e = d.expression()->clone();
    expand(e);
    for (const auto& term : terms_of(e)) {
      for (const auto& i : externals_of(split(term).second, head_slots))
        if (!head_slots.contains(i))
          reject("a definition's expression has a free index its head lacks");
    }
    for (const auto& term : terms_of(e))
      for (const auto& f : split(term).second) {
        const auto& base =
            f->is<TensorPower>() ? f->as<TensorPower>().base() : f;
        if (!base->is<Tensor>()) continue;
        for (std::size_t l = 0; l != defs.size(); ++l)
          if (definition_of(base->as<Tensor>(), Definitions{defs[l]}))
            uses[k].push_back(l);
      }
  }
  // depth-first search for a cycle
  container::svector<int> state(defs.size(), 0);  // 0 new, 1 open, 2 done
  const auto visit = [&](const auto& self, std::size_t k) -> void {
    if (state[k] == 2) return;
    if (state[k] == 1) reject("a definition depends on itself");
    state[k] = 1;
    for (const auto l : uses[k]) self(self, l);
    state[k] = 2;
  };
  for (std::size_t k = 0; k != defs.size(); ++k) visit(visit, k);
}

/// the derivative of @p E by the one variable @p x : the sum over the
/// spellings of x of their phase times the derivative by each, through the
/// definitions @p defs ; the indices @p known are external wherever they
/// stand
ExprPtr by_one(const ExprPtr& E, const ExprPtr& x, const Definitions& defs,
               const container::set<Index>& known) {
  const std::string where = "deriv::diff: ";
  auto e = E->clone();
  expand(e);
  check_variable(x, where);
  check_definitions(defs, x, where);
  if (!e->is_cnumber()) throw Exception(where + "the function holds operators");
  auto terms = split_terms(e, /*with_externals=*/x->is<Tensor>(), known);
  ExprPtr result;
  if (x->is<Variable>()) {
    result = one_way(terms, x, defs);
  } else {
    // the spellings of x, and those of its adjoint where that holds the same
    // elements and the bra-ket exchange is no operation of x
    const auto& t = x->as<Tensor>();
    const auto ys =
        all_spellings(x, adjoint_holds_elements(t) && !exchanges(t));
    rename_apart(terms, x);
    check(terms, x, ys, where);
    // the leaves the definitions differentiate through are checked as those
    // of the function are
    for (const auto& d : defs) {
      // the head's names are dummies: renamed apart from x like the rest
      auto head = d.result_as_tensor();
      container::set<std::wstring> taken;
      for (const auto& i : get_used_indices(x)) taken.emplace(i.label());
      for (const auto& i : get_used_indices(d.expression()))
        taken.emplace(i.label());
      for (const auto& i : head.slots()) taken.emplace(i.label());
      const auto unused = unused_label_validator(taken);
      container::map<Index, Index> fresh;
      for (const auto& i : head.slots()) fresh.emplace(i, temporary(i, unused));
      head.transform_indices(fresh);
      auto d_terms = split_terms(instantiate(d, head, std::move(taken)), true);
      rename_apart(d_terms, x);
      check(d_terms, x, ys, where);
    }
    for (const auto& [y, phase] : ys) {
      auto by_y = ex<Constant>(phase) *
                  block_sum(one_way(terms, y, defs), y->as<Tensor>());
      result = result ? result + by_y : by_y;
    }
  }
  expand(result);
  rapid_simplify(result);
  return result;
}

}  // namespace

ExprPtr diff(const ExprPtr& E, const container::svector<ExprPtr>& variables) {
  return diff(E, variables, Definitions{});
}

ExprPtr diff(const ExprPtr& E, const container::svector<ExprPtr>& variables,
             const container::svector<ResultExpr>& definitions) {
  if (!E || std::ranges::any_of(variables, [](const auto& x) { return !x; }))
    throw Exception("deriv::diff: null argument");
  if (variables.empty()) return E->clone();
  // by_one() works on a clone of its input; the indices of each variable are
  // free in the derivatives by it, even where they stand in two slots or more
  auto result = E;
  container::set<Index> known;
  for (const auto& x : variables) {
    result = by_one(result, x, definitions, known);
    for (const auto& i : get_used_indices(x)) known.insert(i);
  }
  return result;
}

namespace {

/// the bases of the slots of a block: bra, ket, auxiliary
struct Block {
  container::svector<IndexBasis> bra, ket, aux;
  bool operator==(const Block&) const = default;
};

Block block_of(const Tensor& t) {
  Block result;
  for (const auto& i : t.bra()) result.bra.push_back(i.basis());
  for (const auto& i : t.ket()) result.ket.push_back(i.basis());
  for (const auto& i : t.aux()) result.aux.push_back(i.basis());
  return result;
}

/// the images of the block @p b of the tensor @p t under the operations of
/// t, one per operation: the permutations of the slots of each side
/// (Symmetry::Antisymm) or of the columns (column symmetry, equal ranks),
/// each also followed by the exchange of bra and ket where diff() reads the
/// leaves at the exchanged position as elements of t: as an operation
/// (exchanges()), or as the adjoint spelling of a tensor with bra-ket
/// symmetry
container::svector<Block> images(const Tensor& t, const Block& b) {
  const bool antisymm = t.symmetry() == Symmetry::Antisymm;
  const bool columns = column_symmetric(t) && b.bra.size() == b.ket.size();
  const bool exchange = t.braket_symmetry() != BraKetSymmetry::Nonsymm &&
                        adjoint_holds_elements(t);
  container::svector<std::size_t> order_B(b.bra.size()), order_K(b.ket.size());
  std::iota(order_B.begin(), order_B.end(), std::size_t{0});
  std::iota(order_K.begin(), order_K.end(), std::size_t{0});
  const auto permuted = [](const auto& slots, const auto& order) {
    container::svector<IndexBasis> result;
    for (const auto k : order) result.push_back(slots[k]);
    return result;
  };
  container::svector<Block> result;
  do {
    do {
      Block image{permuted(b.bra, order_B),
                  permuted(b.ket, columns ? order_B : order_K), b.aux};
      result.push_back(image);
      if (exchange) {
        std::swap(image.bra, image.ket);
        result.push_back(std::move(image));
      }
    } while (antisymm && std::next_permutation(order_K.begin(), order_K.end()));
  } while ((antisymm || columns) &&
           std::next_permutation(order_B.begin(), order_B.end()));
  return result;
}

/// whether the blocks @p l and @p r share elements: the same ranks, and
/// overlapping spaces in every slot
bool overlap(const Block& l, const Block& r) {
  const auto registry = get_default_index_basis_registry();
  SEQUANT_ASSERT(registry);
  const auto side = [&registry](const auto& a, const auto& b) {
    return std::ranges::equal(
        a, b, [&registry](const IndexBasis& x, const IndexBasis& y) {
          return !different_instances(x, y) &&
                 registry->valid_intersection(x.space(), y.space());
        });
  };
  return side(l.bra, r.bra) && side(l.ket, r.ket) && side(l.aux, r.aux);
}

}  // namespace

container::svector<BlockGradient> block_gradients(
    const ExprPtr& E, std::wstring_view label,
    const container::svector<ResultExpr>& definitions) {
  const std::string where = "deriv::block_gradients: ";
  if (!E) throw Exception(where + "null argument");
  // the leaves with the label, one per class of blocks, and the labels in
  // use; a leaf spelled as the adjoint, where that holds the elements of the
  // tensor with the label, counts as its adjoint
  container::svector<ExprPtr> leaves;
  container::svector<container::svector<Block>> classes;
  container::set<std::wstring> taken;
  const auto visit = [&](const ExprPtr& expr) {
    auto e = expr->clone();
    expand(e);
    for (const auto& term : terms_of(e))
      for (const auto& f : split(term).second) {
        for (const auto& i : get_used_indices(f)) taken.emplace(i.label());
        auto base = f->is<TensorPower>() ? f->as<TensorPower>().base() : f;
        if (!base->is<Tensor>()) continue;
        if (base->as<Tensor>().label() != label) {
          auto flipped = adjoint(base);
          if (flipped->as<Tensor>().label() != label ||
              !adjoint_holds_elements(flipped->as<Tensor>()))
            continue;
          base = flipped;
        }
        const auto& t = base->as<Tensor>();
        if (std::ranges::any_of(t.slots(), &Index::has_proto_indices))
          throw Exception(where + "a leaf with the label has protoindices");
        const auto b = block_of(t);
        bool known = false;
        for (const auto& c : classes) {
          if (std::ranges::find(c, b) != c.end()) {
            known = true;
            break;
          }
          if (std::ranges::any_of(
                  c, [&b](const Block& other) { return overlap(b, other); }))
            throw Exception(where +
                            "two blocks overlap without being one block");
        }
        if (known) continue;
        leaves.push_back(base);
        classes.push_back(images(t, b));
      }
  };
  visit(E);
  // the leaves that E holds through its definitions
  for (const auto& d : definitions)
    if (d.expression()) visit(d.expression());

  container::svector<BlockGradient> result;
  const auto unused = unused_label_validator(std::move(taken));
  for (std::size_t k = 0; k != leaves.size(); ++k) {
    // the leaf under index names of its own
    auto x = leaves[k]->clone();
    container::svector<Index> bra, ket, aux;
    for (const auto& i : x->as<Tensor>().bra())
      bra.push_back(temporary(i, unused));
    for (const auto& i : x->as<Tensor>().ket())
      ket.push_back(temporary(i, unused));
    for (const auto& i : x->as<Tensor>().aux())
      aux.push_back(temporary(i, unused));
    auto& x_t = x->as<Tensor>();
    x_t.set_bra({bra.begin(), bra.end()});
    x_t.set_ket({ket.begin(), ket.end()});
    x_t.set_aux({aux.begin(), aux.end()});
    const auto& block = classes[k].front();
    const auto n = std::ranges::count(classes[k], block);
    auto gradient = diff(E, {x}, definitions);
    if (gradient->is_zero()) continue;
    gradient = ex<Constant>(rational(1, n)) * gradient;
    expand(gradient);
    rapid_simplify(gradient);
    result.push_back({std::move(x), std::move(gradient)});
  }
  return result;
}

container::svector<BlockGradient> block_gradients(const ExprPtr& E,
                                                  std::wstring_view label) {
  return block_gradients(E, label, {});
}

ExprPtr at_zero(const ExprPtr& E,
                const container::svector<ExprPtr>& variables) {
  const std::string where = "deriv::at_zero: ";
  if (!E || std::ranges::any_of(variables, [](const auto& x) { return !x; }))
    throw Exception(where + "null argument");
  // the scalar variables, and the spellings of each tensor and of its
  // adjoint; and each tensor and its adjoint
  container::svector<ExprPtr> ys, tensors;
  for (const auto& x : variables) {
    check_variable(x, where);
    if (x->is<Variable>()) {
      ys.push_back(x);
    } else {
      for (const auto& [y, phase] : all_spellings(x, true)) ys.push_back(y);
      tensors.push_back(x);
      tensors.push_back(adjoint(x));
    }
  }
  // whether the factor f vanishes: it lies inside some spelling of a
  // variable; throws for a variable under a non-positive power, for a leaf
  // that meets a spelling and lies inside none, and for one whose
  // protoindices do not follow the slots of a spelling it meets
  const auto zero = [&](const ExprPtr& f) {
    const auto [base, exponent] =
        f->is<Power>()
            ? std::pair{f->as<Power>().base(), f->as<Power>().exponent()}
        : f->is<TensorPower>() ? std::pair{f->as<TensorPower>().base(),
                                           f->as<TensorPower>().exponent()}
                               : std::pair{f, rational(1)};
    if (base->is<Tensor>())
      for (const auto& x : tensors)
        check_symmetries(base->as<Tensor>(), x->as<Tensor>(), where);
    bool met = false, held = false;
    for (const auto& y : ys) {
      if (y->is<Variable>()) {
        held = held || (base->is<Variable>() && base->as<Variable>().label() ==
                                                    y->as<Variable>().label());
        continue;
      }
      if (!leaf_of(base, y->as<Tensor>())) continue;
      check_layout(base->as<Tensor>(), y->as<Tensor>(), where);
      met = true;
      held = held || std::ranges::equal(
                         y->as<Tensor>().slots(), base->as<Tensor>().slots(),
                         [](const Index& y_l, const Index& l) {
                           return includes(y_l.basis(), l.basis());
                         });
    }
    if (met && !held)
      throw Exception(where +
                      "a leaf holds elements of a variable and lies inside "
                      "none of the variables");
    if (held && exponent <= 0)
      throw Exception(where + "a variable stands under a non-positive power");
    return held;
  };
  auto e = E->clone();
  expand(e);
  container::vector<ExprPtr> result;
  for (const auto& term : terms_of(e)) {
    // every factor, so that each is checked
    bool vanishes = false;
    for (const auto& f : split(term).second) vanishes = zero(f) || vanishes;
    if (!vanishes) result.push_back(term);
  }
  auto r =
      result.empty() ? ex<Constant>(0) : ex<Sum>(result.begin(), result.end());
  rapid_simplify(r);
  return r;
}

}  // namespace sequant::deriv
