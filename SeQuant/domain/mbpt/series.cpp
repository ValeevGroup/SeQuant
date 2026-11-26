#include <SeQuant/domain/mbpt/series.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

#include <string>
#include <utility>

namespace sequant::mbpt {

namespace detail {

const ExprPtr& base_of(const ExprPtr& atom) {
  if (atom->is<TensorPower>()) return atom->as<TensorPower>().base();
  if (atom->is<Power>()) return atom->as<Power>().base();
  return atom;
}

bool is_operator(const Expr& e) {
  return dynamic_cast<const FOperatorBase*>(&e) ||
         dynamic_cast<const BOperatorBase*>(&e);
}

ExprPtr tensor_like(const Tensor& model, std::wstring_view label,
                    container::svector<Index> b, container::svector<Index> k,
                    container::svector<Index> a) {
  return ex<Tensor>(std::wstring(label), bra(std::move(b)), ket(std::move(k)),
                    aux(std::move(a)), model.symmetry(),
                    model.braket_symmetry(), model.column_symmetry());
}

ExprPtr substituted(const ExprPtr& expr,
                    const std::function<ExprPtr(const Tensor&)>& make) {
  auto result = expr->clone();
  for_each_atom(result, [&make](ExprPtr& leaf) {
    if (!leaf->is<Tensor>()) return;
    if (auto replacement = make(leaf->as<Tensor>()))
      leaf = std::move(replacement);
  });
  return result;
}

}  // namespace detail

namespace {

using Indices = container::svector<Index>;

/// the k-th coefficient of the leaf @p t , its base label being @p label ,
/// with the perturbation indices @p z after its own auxiliary indices; that
/// of p⁺ is the adjoint of that of p
ExprPtr coefficient(const Tensor& t, std::wstring_view label, std::size_t k,
                    const Indices& z) {
  auto c = detail::decorate_with_pert_order(label, static_cast<int>(k));
  if (is_adjoint_label(t.label())) toggle_adjoint_label(c);
  Indices aux(t.aux().begin(), t.aux().end());
  aux.insert(aux.end(), z.begin(), z.end());
  return detail::tensor_like(t, c, {t.bra().begin(), t.bra().end()},
                             {t.ket().begin(), t.ket().end()}, std::move(aux));
}

}  // namespace

namespace detail {

std::function<Index(const IndexBasis&)> fresh_indices(const ExprPtr& expr) {
  return [unused = unused_label_validator(*expr)](const IndexBasis& basis) {
    Index result;
    do result = Index::make_tmp_index(basis);
    while (!unused(result));
    return result;
  };
}

void check_series(const ExprPtr& expr, std::wstring_view label,
                  const ExprPtr& strength, std::size_t order, bool expanded) {
  const std::string where = "mbpt::series: ";
  const auto reject = [&where](const std::string& text) {
    throw Exception(where + text);
  };
  if (!expr || !strength) reject("null argument");
  if (label.empty() || is_adjoint_label(label))
    reject("the label must be nonempty and not an adjoint label");
  const bool scalar = strength->is<Variable>();
  const auto* s = strength->is<Tensor>() ? &strength->as<Tensor>() : nullptr;
  if (!scalar &&
      !(s && s->bra_rank() == 0 && s->ket_rank() == 0 && s->aux_rank() == 1))
    reject(
        "the strength must be a Variable or a tensor with one auxiliary index "
        "and no other slot");
  if (scalar && strength->as<Variable>().conjugated())
    reject("the strength must not be conjugated");
  if (s && s->aux().begin()->has_proto_indices())
    reject("the index of the strength must not carry protoindices");
  if (order > 9)
    reject("the order must be at most 9, not " + std::to_string(order));

  container::set<std::wstring> coefficients;
  for (std::size_t k = 1; k <= order; ++k)
    coefficients.emplace(
        detail::decorate_with_pert_order(label, static_cast<int>(k)));
  const auto strength_label = strip_adjoint_label(
      scalar ? strength->as<Variable>().label() : s->label());
  if (strength_label == label ||
      coefficients.contains(std::wstring(strength_label)))
    reject(
        "the strength must not be labelled like a tensor to expand or a "
        "coefficient");
  for_each_atom(expr, [&](const ExprPtr& leaf) {
    if (is_operator(*leaf))
      reject("the expression holds operators; expand its tensor form");
    const auto& t = base_of(leaf);
    if (expanded && t->is<Variable>() && t->as<Variable>().label() == label)
      reject("the expression holds a Variable labelled " + toUtf8(label));
    if (!t->is<Tensor>()) return;
    const auto base = strip_adjoint_label(t->as<Tensor>().label());
    if (expanded && leaf->is<TensorPower>() && base == label)
      reject("a tensor to expand stands under a power");
    if (coefficients.contains(std::wstring(base)))
      reject("the expression already holds a tensor labelled " +
             toUtf8(t->as<Tensor>().label()));
  });
}

ExprPtr series_term(const Tensor& t, std::wstring_view label,
                    const ExprPtr& strength, std::size_t k,
                    const std::function<Index(const IndexBasis&)>& fresh) {
  if (k == 0) return t.clone();
  const auto* s = strength->is<Tensor>() ? &strength->as<Tensor>() : nullptr;
  Indices z;
  if (s)
    for (std::size_t j = 0; j != k; ++j)
      z.push_back(fresh(s->aux().begin()->basis()));
  rational scale = 1;
  for (std::size_t j = 2; j <= k; ++j) scale /= j;
  container::svector<ExprPtr> factors{coefficient(t, label, k, z)};
  if (!s)
    factors.push_back(k == 1 ? strength->clone()
                             : ex<Power>(strength->clone(), rational(k)));
  for (const auto& zj : z)
    factors.push_back(tensor_like(*s, s->label(), {}, {}, {zj}));
  return ex<Product>(scale, std::move(factors));
}

}  // namespace detail

ExprPtr series(const ExprPtr& expr, std::wstring_view label,
               const ExprPtr& strength, std::size_t order) {
  detail::check_series(expr, label, strength, order);
  // perturbation indices are named apart from every index of expr
  const auto fresh = detail::fresh_indices(expr);
  return detail::substituted(expr, [&](const Tensor& t) -> ExprPtr {
    if (strip_adjoint_label(t.label()) != label || order == 0) return nullptr;
    container::svector<ExprPtr> terms;
    for (std::size_t k = 0; k <= order; ++k)
      terms.push_back(detail::series_term(t, label, strength, k, fresh));
    return ex<Sum>(std::move(terms));
  });
}

}  // namespace sequant::mbpt
