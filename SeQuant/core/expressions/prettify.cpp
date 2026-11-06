#include <SeQuant/core/expressions/prettify.hpp>

#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/lexicographic_order.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <memory>

#include <range/v3/algorithm/all_of.hpp>

namespace sequant {

namespace {

/// @return true if @p expr is a tensor, or a product of tensors, scalars and
/// (recursively) nested such products; only these are canonicalized by
/// prettify(), because their canonicalization preserves the parenthesization
bool is_tensor_tree(const ExprPtr& expr) {
  if (expr.is<AbstractTensor>()) return true;
  if (!expr.is<Product>()) return false;
  for (const auto& factor : expr.as<Product>().factors()) {
    if (factor->is_scalar()) continue;
    if (!is_tensor_tree(factor)) return false;
  }
  return true;
}

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
  for (std::size_t f = 0; f != factors.size(); ++f) {
    if (!factors[f]->is<Sum>()) continue;
    const auto sum_indices = get_used_indices(*factors[f]);
    container::set<Index> sum_named = named;
    bool other_dummies = false;
    for (std::size_t g = 0; g != factors.size() && !other_dummies; ++g) {
      if (g == f) continue;
      for (auto&& idx : get_used_indices(*factors[g])) {
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

/// canonicalizes the tensor trees (see is_tensor_tree()) of @p expr , which
/// are the summands of the (nested) sums of @p expr or @p expr itself, and
/// the summands of the sums that are factors of its products (see
/// recanonicalize_sum_factors()), with @p opts
/// @note sums are neither flattened nor are their summands combined, i.e.
/// the structure of @p expr is preserved
void recanonicalize(ExprPtr& expr, const NamedIndexSpec& named_indices,
                    const CanonicalizeOptions& opts) {
  if (expr.is<Sum>()) {
    // N.B. mutable iteration resets the memoized hash of the sum
    for (auto& summand : *expr) recanonicalize(summand, named_indices, opts);
  } else if (is_tensor_tree(expr)) {
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

/// @return the prettified image of @p term under the permutation of the
/// result indices of @p group that is first in the order of @p opts
ExprPtr prettify_term_for_result(
    const ExprPtr& term,
    const std::vector<std::pair<container::map<Index, Index>, int>>& group,
    const NamedIndexSpec& named_indices, const PrettifyOptions& opts) {
  const auto cmp = ordering(opts);
  ExprPtr best;
  for (const auto& [map, sign] : group) {
    ExprPtr image =
        map.empty() ? term->clone() : transform_expr(term, map, sign);
    prettify_in_place(image, named_indices, opts);
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
  const auto group = detail::result_symmetry_group(expr);
  const auto& rhs = expr.expression();

  // shares the expression with expr until it is replaced by the prettified
  // copy
  ResultExpr result = expr;
  if (!rhs || group.size() == 1 || !opts.recanonicalize) {
    result.expression() = prettify(rhs, named, opts);
    return result;
  }

  // N.B. the summands are replaced in place, so that nested sums stay nested
  auto present = [&](ExprPtr& term) {
    // a (anti)symmetrizer of the result times a sum, S * (X_1 + X_2 + ...):
    // S * X = S * P(X) (up to the phase of P) for every P of the group, so
    // each X_k is presented by its image that comes first
    if (term.is<Product>()) {
      const auto& product = term.as<Product>();
      ExprPtr const* sum = nullptr;
      bool symmetrizer = false, other = false;
      for (const auto& factor : product.factors()) {
        if (factor.is<Sum>() && !sum) {
          sum = &factor;
        } else if (factor.is<Tensor>() &&
                   (factor.as<Tensor>().label() == reserved::symm_label() ||
                    factor.as<Tensor>().label() ==
                        reserved::antisymm_label()) &&
                   ranges::all_of(factor.as<Tensor>().const_indices(),
                                  [&](const Index& idx) {
                                    return named.indices().contains(idx);
                                  })) {
          symmetrizer = true;
        } else {
          other = true;
        }
      }
      if (sum && symmetrizer && !other) {
        term = term->clone();
        // N.B. mutable iteration resets the memoized hashes
        for (auto& factor : *term) {
          if (!factor.is<Sum>()) continue;
          for (auto& summand : *factor)
            summand = prettify_term_for_result(summand, group, named, opts);
        }
        return;
      }
    }
    // N.B. a term that is not under a symmetrizer is not replaced by an
    // image: that the result has the symmetry does not make its terms have
    // it, and t + u = P(t) + u does not hold in general
    term = prettify(term, named, opts);
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
