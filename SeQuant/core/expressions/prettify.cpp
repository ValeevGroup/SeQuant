#include <SeQuant/core/expressions/prettify.hpp>

#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/lexicographic_order.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <memory>

#include <range/v3/algorithm/none_of.hpp>

namespace sequant {

namespace {

void recanonicalize(ExprPtr& expr, const NamedIndexSpec& named_indices,
                    const CanonicalizeOptions& opts);

/// recanonicalizes the sums that are factors of @p product , e.g. the sum
/// that a symmetrizer multiplies, with the indices that a sum shares with the
/// other factors named, too; only done if the other factors have no dummy
/// indices of their own, which the dummy indices of the sum, being renamed
/// independently, could clash with
void recanonicalize_sum_factors(Product& product,
                                const NamedIndexSpec& named_indices,
                                const CanonicalizeOptions& opts) {
  if (named_indices.deduce()) return;
  const auto& named = named_indices.indices();
  // N.B. mutable iteration resets the memoized hash of the product
  container::svector<ExprPtr*> factors;
  for (auto& factor : product) factors.push_back(&factor);
  if (ranges::none_of(factors, [](const ExprPtr* f) { return f->is<Sum>(); }))
    return;
  // N.B. the indices of each factor are collected once
  container::svector<container::set<Index>> factor_indices;
  for (const auto* factor : factors)
    factor_indices.push_back(get_used_indices(*factor));
  for (std::size_t f = 0; f != factors.size(); ++f) {
    if (!factors[f]->is<Sum>()) continue;
    const auto& sum_indices = factor_indices[f];
    container::set<Index> sum_named = named;
    bool other_dummies = false;
    for (std::size_t g = 0; g != factors.size() && !other_dummies; ++g) {
      if (g == f) continue;
      for (auto&& idx : factor_indices[g]) {
        if (sum_indices.contains(idx))
          sum_named.insert(idx);
        else if (!named.contains(idx)) {
          other_dummies = true;
          break;
        }
      }
    }
    if (other_dummies) continue;
    recanonicalize(*factors[f],
                   detail::explicit_named_indices(std::move(sum_named)), opts);
  }
}

/// canonicalizes the tensor trees (see detail::is_tensor_tree(); their
/// canonicalization preserves the parenthesization) of @p expr , which are the
/// summands of the (nested) sums of @p expr or @p expr itself, and
/// the summands of the sums that are factors of its products (see
/// recanonicalize_sum_factors()), with @p opts
/// @note sums are neither flattened nor are their summands combined, i.e.
/// the structure of @p expr is preserved
void recanonicalize(ExprPtr& expr, const NamedIndexSpec& named_indices,
                    const CanonicalizeOptions& opts) {
  if (expr.is<Sum>()) {
    // N.B. mutable iteration resets the memoized hash of the sum
    for (auto& summand : *expr) recanonicalize(summand, named_indices, opts);
  } else if (detail::is_tensor_tree(expr)) {
    canonicalize(expr, named_indices, opts);
  } else if (expr.is<Product>()) {
    recanonicalize_sum_factors(expr.as<Product>(), named_indices, opts);
  }
}

CanonicalizeOptions canonicalize_options(const PrettifyOptions& opts) {
  // N.B. as when canonicalizing the summands of a sum, the labels of the named
  // indices are taken into account
  return CanonicalizeOptions{
      .method = CanonicalizationMethod::Complete,
      .ignore_named_index_labels =
          CanonicalizeOptions::IgnoreNamedIndexLabel::No,
      .presentation = std::make_shared<const CanonicalizationPresentation>(
          opts.presentation)};
}

LexicographicExprCompare ordering(const PrettifyOptions& opts) {
  return LexicographicExprCompare(LexicographicCompareOptions{
      .cardinal_tensor_labels =
          opts.presentation.cardinal_tensor_labels.value_or(
              default_cardinal_tensor_labels())});
}

/// prettifies @p expr in place
void prettify_in_place(ExprPtr& expr, const NamedIndexSpec& named_indices,
                       const PrettifyOptions& opts) {
  if (opts.recanonicalize)
    recanonicalize(expr, named_indices, canonicalize_options(opts));
  // N.B. the canonical marks need no reset: canonicalization with
  // presentation settings does not mark, and a mark copied by the clone is
  // not valid anymore once the expression has changed
  sort_expr(*expr, ordering(opts));
}

/// @return the prettified image of @p term under the permutation of @p group
/// (see detail::symmetrizer_group()) that is first in the order of @p opts ,
/// with the phase of the permutation
ExprPtr prettify_term_for_result(const ExprPtr& term,
                                 const detail::IndexPermutationGroup& group,
                                 const NamedIndexSpec& named_indices,
                                 const PrettifyOptions& opts) {
  const auto cmp = ordering(opts);
  ExprPtr best;
  for (const auto& [map, sign] : group) {
    ExprPtr image =
        map.empty() ? term->clone() : transform_expr(term, map, sign);
    prettify_in_place(image, named_indices, opts);
    if (image->is_zero()) return ex<Constant>(0);
    if (!best || cmp(image, best)) best = std::move(image);
  }
  return best;
}

}  // namespace

ExprPtr prettify(const ExprPtr& expr, const NamedIndexSpec& named_indices,
                 const PrettifyOptions& opts) {
  if (!expr) return expr;
  ExprPtr result = expr->clone();
  prettify_in_place(result, named_indices, opts);
  return result;
}

ResultExpr prettify(const ResultExpr& expr, const PrettifyOptions& opts) {
  const auto named = named_indices(expr);
  const auto& rhs = expr.expression();

  // shares the expression with expr until it is replaced by the prettified
  // copy
  ResultExpr result = expr;
  if (!rhs || !opts.recanonicalize) {
    result.expression() = prettify(rhs, named, opts);
    return result;
  }

  // the parts of the terms that a symmetrizer of the result multiplies are
  // presented by their images that come first (see
  // detail::for_each_symmetrized_term()); the other terms are not replaced by
  // images
  // N.B. the summands are replaced in place, so that nested sums stay nested
  auto present = [&](ExprPtr& term) {
    bool replaced = false;
    detail::for_each_symmetrized_term(
        term, expr,
        [&](ExprPtr& part, const container::svector<Tensor>& symmetrizers) {
          const auto group = detail::symmetrizers_group(symmetrizers);
          if (group.size() == 1) return;
          part = prettify_term_for_result(part, group, named, opts);
          replaced = true;
        });
    if (!replaced) term = prettify(term, named, opts);
  };
  ExprPtr prettified = rhs->clone();
  if (prettified.is<Sum>()) {
    for (auto& summand : *prettified) present(summand);
  } else {
    present(prettified);
  }
  sort_expr(*prettified, ordering(opts));
  result.expression() = std::move(prettified);
  return result;
}

}  // namespace sequant
