#include <SeQuant/domain/mbpt/spbasis_rotation_gradient.hpp>

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/tensor_power.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace sequant::mbpt {

namespace {

/// appends the factors of the term @p e of an expanded expression to
/// @p factors and multiplies @p scalar by its scalars
void split(const ExprPtr& e, Constant::scalar_type& scalar,
           container::svector<ExprPtr>& factors) {
  if (e->is<Product>()) {
    scalar *= e->as<Product>().scalar();
    for (const auto& f : e->as<Product>().factors()) split(f, scalar, factors);
  } else if (e->is<Constant>()) {
    scalar *= e->as<Constant>().value();
  } else {
    SEQUANT_ASSERT(!e->is<Sum>());
    factors.push_back(e);
  }
}

[[noreturn]] void reject(const char* text) {
  throw Exception(std::string("mbpt::spbasis_rotation_gradient: ") + text);
}

/// what spbasis_rotation_gradient() is given
struct Rotation {
  const container::svector<std::wstring>& labels;
  const Index& p;
  const Index& q;
  const container::svector<ResultExpr>& definitions;

  /// whether @p f is a tensor labelled as one to rotate, or the adjoint
  /// spelling of one: every slot transforms alike
  bool labelled(const ExprPtr& f) const {
    const auto named = [this](const ExprPtr& leaf) {
      return std::ranges::find(labels, leaf->as<Tensor>().label()) !=
             labels.end();
    };
    return f->is<Tensor>() && (named(f) || named(adjoint(f)));
  }

  /// the definition of @p f , or nullptr if it has none
  const ResultExpr* definition_of(const ExprPtr& f) const {
    return f->is<Tensor>()
               ? deriv::detail::definition_of(f->as<Tensor>(), definitions)
               : nullptr;
  }

  /// appends to @p result the terms of the derivative of the expanded @p e .
  /// The slots that hold an index of @p fixed do not transform. @p head holds
  /// the indices of the head of the definition that @p e is the expression of,
  /// @p depth definitions below the function; @p outer holds the labels of
  /// the indices of the term that expression stands in
  void differentiate(const ExprPtr& e, const container::set<Index>& fixed,
                     const container::set<Index>& head, std::size_t depth,
                     const container::set<std::wstring>& outer,
                     container::vector<ExprPtr>& result) const {
    const auto terms = e->is<Sum>() ? container::svector<ExprPtr>(
                                          e->as<Sum>().summands().begin(),
                                          e->as<Sum>().summands().end())
                                    : container::svector<ExprPtr>{e};
    for (const auto& term : terms) {
      Constant::scalar_type scalar = 1;
      container::svector<ExprPtr> factors;
      split(term, scalar, factors);
      // in how many bra, ket and auxiliary slots of the term an index stands,
      // the indices that are protoindices in it, and the labels it uses
      container::map<Index, std::int64_t> slots;
      container::set<Index> protos;
      auto taken = outer;
      taken.emplace(p.label());
      taken.emplace(q.label());
      for (const auto& f : factors) {
        if (f->is<Tensor>() &&
            (f->as<Tensor>().label() == reserved::antisymm_label() ||
             f->as<Tensor>().label() == reserved::symm_label()))
          reject("the function holds an Â or Ŝ");
        for (const auto& [i, count] : get_used_indices_with_counts(*f)) {
          if (i.label() == p.label() || i.label() == q.label())
            reject("the function uses the label of p or q");
          slots[i] += count.nonproto();
          if (count.proto > 0) protos.insert(i);
          taken.emplace(i.label());
        }
      }
      for (std::size_t k = 0; k != factors.size(); ++k) {
        const auto& f = factors[k];
        const bool power = f->is<TensorPower>();
        const auto& leaf = power ? f->as<TensorPower>().base() : f;
        if (power && (labelled(leaf) || definition_of(leaf)))
          reject("a leaf that rotates stands under a TensorPower");
        if (!leaf->is<Tensor>()) continue;
        if (!definition_of(leaf) && definition_of(adjoint(leaf)))
          reject("a leaf is spelled as the adjoint of one with a definition");
        const auto* definition = definition_of(f);
        if (!labelled(f) && !definition) {
          for (const auto& u : std::as_const(leaf->as<Tensor>()).braket())
            if (head.contains(u))
              reject(
                  "a tensor of the expression of a definition holds an index "
                  "of the head and does not rotate");
          continue;
        }
        const auto& t = f->as<Tensor>();
        const auto nb = t.bra_rank();
        for (std::size_t s = 0; s != nb + t.ket_rank(); ++s) {
          const auto& u = s < nb ? t.bra()[s] : t.ket()[s - nb];
          if (fixed.contains(u)) continue;
          if (u.has_proto_indices() || protos.contains(u))
            reject(
                "a slot of a leaf that rotates has protoindices or is a "
                "protoindex");
          if (slots.at(u) > 2)
            reject(
                "a slot of a leaf that rotates holds an index that stands in "
                "more than two slots");
          if (!includes(u.basis(), p.basis())) {
            if (get_default_index_basis_registry()->valid_intersection(
                    u.space(), p.space()))
              reject(
                  "a slot of a leaf that rotates overlaps the space of p "
                  "without including it");
            continue;
          }
          // a free index stays, with its Kronecker delta with p; a summed one
          // becomes p elsewhere
          const bool free = slots.at(u) == 1;
          const container::map<Index, Index> u_to_p =
              free ? container::map<Index, Index>{}
                   : container::map<Index, Index>{{u, p}};
          auto rotated = f->clone();
          auto& r = rotated->as<Tensor>();
          container::svector<Index> side(
              s < nb ? r.bra().begin() : r.ket().begin(),
              s < nb ? r.bra().end() : r.ket().end());
          side[s < nb ? s : s - nb] = q;
          if (s < nb)
            r.set_bra({side.begin(), side.end()});
          else
            r.set_ket({side.begin(), side.end()});
          rotated = transform_expr(rotated, u_to_p);
          auto product = ex<Product>(scalar, ExprPtrList{});
          for (std::size_t j = 0; j != factors.size(); ++j)
            product->as<Product>().append(
                1, j == k ? rotated : transform_expr(factors[j], u_to_p));
          if (free)
            product->as<Product>().append(
                1, s < nb ? make_kronecker(u, p) : make_kronecker(p, u));
          result.push_back(product);
        }
        if (!definition) continue;
        // through the summed indices of the leaf's expression: its own slots
        // are those above
        if (depth == definitions.size())
          reject("a definition depends on itself");
        container::set<Index> of_head(t.braket().begin(), t.braket().end());
        auto held = fixed;
        held.insert(of_head.begin(), of_head.end());
        container::vector<ExprPtr> inner;
        differentiate(deriv::detail::instantiate(*definition, t, taken), held,
                      of_head, depth + 1, taken, inner);
        for (const auto& g : inner) {
          auto product = ex<Product>(scalar, ExprPtrList{});
          for (std::size_t j = 0; j != factors.size(); ++j)
            product->as<Product>().append(1, j == k ? g : factors[j]);
          result.push_back(product);
        }
      }
    }
  }
};

}  // namespace

ExprPtr spbasis_rotation_gradient(
    const ExprPtr& L, const container::svector<std::wstring>& labels,
    const Index& p, const Index& q,
    const container::svector<ResultExpr>& definitions,
    const container::svector<Index>& fixed) {
  if (!L) reject("null argument");
  if (p.has_proto_indices() || q.has_proto_indices())
    reject("p and q must not have protoindices");
  if (p.label() == q.label()) reject("p and q are one index");
  auto e = L->clone();
  expand(e);
  if (!e->is_cnumber()) reject("the function holds operators");
  for (const auto& d : definitions)
    if (!d.expression() || !d.expression()->is_cnumber())
      reject("a definition holds operators");
  const Rotation rotation{labels, p, q, definitions};
  container::vector<ExprPtr> result;
  rotation.differentiate(e, container::set<Index>(fixed.begin(), fixed.end()),
                         {}, 0, {}, result);
  auto r =
      result.empty() ? ex<Constant>(0) : ex<Sum>(result.begin(), result.end());
  rapid_simplify(r);
  return r;
}

ExprPtr spbasis_rotation_gradient(
    const ExprPtr& L, const container::svector<std::wstring>& labels,
    const Index& p, const Index& q) {
  return spbasis_rotation_gradient(L, labels, p, q,
                                   container::svector<ResultExpr>{});
}

}  // namespace sequant::mbpt
