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
#include <SeQuant/core/reserved.hpp>
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
  // a leaf spelled as the adjoint of a tensor to rotate rotates too: every
  // slot transforms alike
  const auto rotates = [&labels](const ExprPtr& f) {
    const auto named = [&labels](const ExprPtr& leaf) {
      return std::ranges::find(labels, leaf->as<Tensor>().label()) !=
             labels.end();
    };
    return f->is<Tensor>() && (named(f) || named(adjoint(f)));
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
    // in how many bra, ket and auxiliary slots of the term an index stands,
    // and the indices that are protoindices in it
    container::map<Index, std::int64_t> slots;
    container::set<Index> protos;
    for (const auto& f : factors) {
      if (f->is<Tensor>() &&
          (f->as<Tensor>().label() == reserved::antisymm_label() ||
           f->as<Tensor>().label() == reserved::symm_label()))
        reject("the function is not a scalar");
      for (const auto& [i, count] : get_used_indices_with_counts(*f)) {
        if (i.label() == p.label() || i.label() == q.label())
          reject("the function uses the label of p or q");
        slots[i] += count.nonproto();
        if (count.proto > 0) protos.insert(i);
      }
    }
    for (const auto& [i, n] : slots)
      if (n == 1 && !protos.contains(i)) reject("the function is not a scalar");
    for (std::size_t k = 0; k != factors.size(); ++k) {
      const auto& f = factors[k];
      if (f->is<TensorPower>() && rotates(f->as<TensorPower>().base()))
        reject("a leaf to rotate stands under a TensorPower");
      if (!rotates(f)) continue;
      const auto& t = f->as<Tensor>();
      const auto nb = t.bra_rank();
      for (std::size_t s = 0; s != nb + t.ket_rank(); ++s) {
        const auto& u = s < nb ? t.bra()[s] : t.ket()[s - nb];
        if (u.has_proto_indices() || protos.contains(u))
          reject(
              "a slot of a leaf to rotate has protoindices or is a "
              "protoindex");
        if (!includes(u.basis(), p.basis())) {
          if (registry->valid_intersection(u.space(), p.space()))
            reject(
                "a slot of a leaf to rotate overlaps the space of p without "
                "including it");
          continue;
        }
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
