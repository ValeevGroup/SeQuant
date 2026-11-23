#include <SeQuant/domain/mbpt/rotation_gradient.hpp>

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/tensor_power.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>

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

}  // namespace

ExprPtr rotation_gradient(const ExprPtr& L,
                          const container::svector<std::wstring>& labels,
                          const Index& p, const Index& q) {
  const auto reject = [](const char* text) {
    throw Exception(std::string("mbpt::rotation_gradient: ") + text);
  };
  if (!L) reject("null argument");
  if (p.has_proto_indices() || q.has_proto_indices())
    reject("p and q must not have protoindices");
  if (p.label() == q.label()) reject("p and q are one index");
  auto e = L->clone();
  expand(e);
  if (!e->is_cnumber()) reject("the function holds operators");
  const auto registry = get_default_index_basis_registry();
  SEQUANT_ASSERT(registry);
  const auto rotates = [&labels](const Tensor& t) {
    return std::ranges::find(labels, t.label()) != labels.end();
  };

  container::vector<ExprPtr> result;
  const auto terms = e->is<Sum>() ? container::svector<ExprPtr>(
                                        e->as<Sum>().summands().begin(),
                                        e->as<Sum>().summands().end())
                                  : container::svector<ExprPtr>{e};
  for (const auto& term : terms) {
    Constant::scalar_type scalar = 1;
    container::svector<ExprPtr> factors;
    split(term, scalar, factors);
    // in how many bra, ket and auxiliary slots of the term an index stands
    container::map<Index, std::int64_t> slots;
    for (const auto& f : factors)
      for (const auto& [i, count] : get_used_indices_with_counts(*f)) {
        if (i.label() == p.label() || i.label() == q.label())
          reject("the function uses the label of p or q");
        slots[i] += count.nonproto();
      }
    for (std::size_t k = 0; k != factors.size(); ++k) {
      const auto& f = factors[k];
      if (f->is<TensorPower>() && f->as<TensorPower>().base()->is<Tensor>() &&
          rotates(f->as<TensorPower>().base()->as<Tensor>()))
        reject("a leaf to rotate stands under a TensorPower");
      if (!f->is<Tensor>() || !rotates(f->as<Tensor>())) continue;
      const auto& t = f->as<Tensor>();
      const auto nb = t.bra_rank();
      for (std::size_t s = 0; s != nb + t.ket_rank(); ++s) {
        const auto& u = s < nb ? t.bra()[s] : t.ket()[s - nb];
        if (u.has_proto_indices())
          reject("a slot of a leaf to rotate has protoindices");
        if (!includes(u.basis(), p.basis())) {
          if (registry->valid_intersection(u.space(), p.space()))
            reject(
                "a slot of a leaf to rotate overlaps the space of p without "
                "including it");
          continue;
        }
        if (slots.at(u) < 2) reject("a slot of a leaf to rotate is free");
        // the term with q in this slot and p for u elsewhere
        const container::map<Index, Index> u_to_p{{u, p}};
        auto rotated = transform_expr(f, u_to_p);
        auto& r = rotated->as<Tensor>();
        container::svector<Index> side(
            s < nb ? r.bra().begin() : r.ket().begin(),
            s < nb ? r.bra().end() : r.ket().end());
        side[s < nb ? s : s - nb] = q;
        if (s < nb)
          r.set_bra({side.begin(), side.end()});
        else
          r.set_ket({side.begin(), side.end()});
        auto product = ex<Product>(scalar, ExprPtrList{});
        for (std::size_t j = 0; j != factors.size(); ++j)
          product->as<Product>().append(
              1, j == k ? rotated : transform_expr(factors[j], u_to_p));
        result.push_back(product);
      }
    }
  }
  auto r =
      result.empty() ? ex<Constant>(0) : ex<Sum>(result.begin(), result.end());
  rapid_simplify(r);
  return r;
}

}  // namespace sequant::mbpt
