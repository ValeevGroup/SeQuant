#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_operators.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/range/primitives.hpp>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace sequant {

std::wstring to_latex_align(const ExprPtr& exprptr, size_t max_lines_per_align,
                            size_t max_terms_per_line) {
  return to_latex_align(exprptr, get_default_context().print_options(),
                        max_lines_per_align, max_terms_per_line);
}

std::wstring to_latex_align(const ExprPtr& exprptr, const PrintOptions& opts,
                            size_t max_lines_per_align,
                            size_t max_terms_per_line) {
  std::wstring result = io::latex::to_string(exprptr, opts);
  if (exprptr->is<Sum>()) {
    result.erase(0, 7);  // remove leading  "{ \bigl"
    result.replace(result.size() - 8, 8,
                   L")");  // replace trailing "\bigr) }" with ")"
    result = std::wstring(L"\\begin{align}\n& ") + result;
    // assume no inner sums
    size_t line_counter = 0;
    size_t term_counter = 0;
    std::wstring::size_type pos = 0;
    std::wstring::size_type plus_pos = 0;
    std::wstring::size_type minus_pos = 0;
    bool last_pos_has_plus = false;
    bool have_next_term = true;
    auto insert_into_result_at = [&](std::wstring::size_type at,
                                     const auto& str) {
      SEQUANT_ASSERT(pos != std::wstring::npos);
      result.insert(at, str);
      const auto str_nchar = std::size(str) - 1;  // neglect end-of-string
      pos += str_nchar;
      if (plus_pos != std::wstring::npos) plus_pos += str_nchar;
      if (minus_pos != std::wstring::npos) minus_pos += str_nchar;
      if (pos != plus_pos)
        SEQUANT_ASSERT(plus_pos == result.find(L" + ", plus_pos));
      if (pos != minus_pos)
        SEQUANT_ASSERT(minus_pos == result.find(L" - ", minus_pos));
    };
    while (have_next_term) {
      if (max_lines_per_align > 0 &&
          line_counter == max_lines_per_align) {  // start new align block?
        insert_into_result_at(pos + 1, L"\n\\end{align}\n\\begin{align}\n& ");
        line_counter = 0;
      } else {
        // break the line if needed
        if (term_counter != 0 && term_counter % max_terms_per_line == 0) {
          insert_into_result_at(pos + 1, L"\\\\\n& ");
          ++line_counter;
        }
      }
      // next term, plz
      if (plus_pos == 0 || last_pos_has_plus)
        plus_pos = result.find(L" + ", plus_pos + 1);
      if (minus_pos == 0 || !last_pos_has_plus)
        minus_pos = result.find(L" - ", minus_pos + 1);
      pos = std::min(plus_pos, minus_pos);
      last_pos_has_plus = (pos == plus_pos);
      if (pos != std::wstring::npos)
        ++term_counter;
      else
        have_next_term = false;
    }
  } else {
    result = std::wstring(L"\\begin{align}\n& ") + result;
  }
  result += L"\n\\end{align}";
  return result;
}

std::size_t size(const Expr& expr) { return ranges::size(expr); }

std::size_t size(const ExprPtr& exprptr) {
  if (exprptr) {
    return size(*exprptr);
  }

  return 0;
}

namespace {

/// @return the term that @p expr stands for if it is a Sum of at most one
/// summand or a Product of one factor with a unit scalar, else nullptr
ExprPtr lone_term(const ExprPtr& expr) {
  if (expr->is<Sum>()) {
    const auto& sum = expr->as<Sum>();
    if (sum.empty()) return ex<Constant>(0);
    if (sum.size() == 1) return sum.summand(0);
  } else if (expr->is<Product>()) {
    const auto& product = expr->as<Product>();
    if (product.size() == 1 && product.scalar() == 1) return product.factor(0);
  }
  return nullptr;
}

}  // namespace

ExprPtr& canonicalize(ExprPtr& expr, const NamedIndexSpec& named_indices,
                      const CanonicalizeOptions& opts) {
  // the contexts are read at several points of this call and of the hooks
  const auto contexts_in_effect = pin_default_contexts();
  if (expr->is_canonical(named_indices, opts)) return expr;
  const auto contexts_version = current_contexts_version();
  const auto byproduct = expr->canonicalize(named_indices, opts);
  if (byproduct) {
    if (!byproduct->is<Constant>())
      throw Exception(
          "canonicalize: the byproduct of canonicalization must be a "
          "Constant");
    expr = byproduct * expr;
  }
  // a Sum or Product left with a single term is that term, whose canonical
  // form on its own may differ from its form as a summand or factor
  if (auto term = lone_term(expr)) {
    expr = std::move(term);
    return canonicalize(expr, named_indices, opts);
  }
  expr->mark_canonical(named_indices, opts, contexts_version);
  return expr;
}

ExprPtr canonicalize(ExprPtr&& expr_rv, const NamedIndexSpec& named_indices,
                     const CanonicalizeOptions& opts) {
  canonicalize(expr_rv, named_indices, opts);
  return std::move(expr_rv);
}

NamedIndexSpec named_indices(const ResultExpr& expr) {
  return detail::explicit_named_indices(expr.indices() |
                                        ranges::to<container::set<Index>>);
}

detail::IndexPermutationGroup detail::symmetrizer_group(
    const Tensor& symmetrizer) {
  const bool antisymm = symmetrizer.label() == reserved::antisymm_label();
  SEQUANT_ASSERT(antisymm || symmetrizer.label() == reserved::symm_label());
  const auto& bra = symmetrizer.bra();
  const auto& ket = symmetrizer.ket();
  using Perm = std::vector<std::size_t>;
  auto all_perms = [](std::size_t n) {
    std::vector<Perm> perms;
    Perm p(n);
    std::iota(p.begin(), p.end(), 0);
    do perms.push_back(p);
    while (std::next_permutation(p.begin(), p.end()));
    return perms;
  };
  auto parity = [](const Perm& p) {
    int sign = 1;
    for (std::size_t i = 0; i != p.size(); ++i)
      for (std::size_t j = i + 1; j != p.size(); ++j)
        if (p[i] > p[j]) sign = -sign;
    return sign;
  };

  // pairs of permutations of the bra and of the ket slots, with phase:
  // Â (see mbpt::A_maps()) permutes the bra and the ket independently, with
  // the parities of the permutations as phase; Ŝ (see
  // mbpt::S_replacement_maps()) permutes the columns, i.e. the bra and the ket
  // together, without phase
  std::vector<std::tuple<Perm, Perm, int>> slot_perms;
  if (antisymm) {
    for (const auto& pb : all_perms(bra.size()))
      for (const auto& pk : all_perms(ket.size()))
        slot_perms.emplace_back(pb, pk, parity(pb) * parity(pk));
  } else {
    if (bra.size() != ket.size()) return {{{}, 1}};
    for (const auto& p : all_perms(bra.size()))
      slot_perms.emplace_back(p, p, 1);
  }

  IndexPermutationGroup group;
  for (const auto& [pb, pk, sign] : slot_perms) {
    container::map<Index, Index> map;
    for (std::size_t k = 0; k != bra.size(); ++k)
      if (bra[k].nonnull() && bra[k] != bra[pb[k]])
        map.emplace(bra[k], bra[pb[k]]);
    for (std::size_t k = 0; k != ket.size(); ++k)
      if (ket[k].nonnull() && ket[k] != ket[pk[k]])
        map.emplace(ket[k], ket[pk[k]]);
    // only the permutations that keep the space of every index: they relate
    // elements of the same block
    if (ranges::all_of(map, [](const auto& source_target) {
          return source_target.first.space() == source_target.second.space();
        }))
      group.emplace_back(std::move(map), sign);
  }
  return group;
}

namespace {

/// @return the canonical form of @p term (a part of the expression of a
/// result that a symmetrizer multiplies, see
/// detail::for_each_symmetrized_term()) that is the same for all its images
/// under @p group , the group of the symmetrizer (see
/// detail::symmetrizer_group()), with the phase of the permutation absorbed:
/// of the canonical forms of the images the one that is smallest by hash
/// (then by its text) is chosen; the indices of the result keep their labels
ExprPtr canonicalize_term_for_result(const ExprPtr& term,
                                     const detail::IndexPermutationGroup& group,
                                     const NamedIndexSpec& named,
                                     const CanonicalizeOptions& opts) {
  // the key ignores the scalar, so that images that differ only by it
  // compare equal
  auto unscaled = [](const ExprPtr& image) {
    if (!image->is<Product>() || image->as<Product>().scalar() == 1)
      return image;
    ExprPtr result = image->clone();
    result->as<Product>().scale(1 / image->as<Product>().scalar());
    return result;
  };
  ExprPtr best;
  std::optional<std::pair<std::size_t, std::wstring>> best_key;
  for (const auto& [map, sign] : group) {
    ExprPtr image =
        map.empty() ? term->clone() : transform_expr(term, map, sign);
    canonicalize(image, named, opts);
    // an image that vanishes (e.g. by the symmetry of its tensors) is the
    // image of a term that vanishes
    if (image->is_zero()) return ex<Constant>(0);
    std::pair<std::size_t, std::wstring> key{unscaled(image)->hash_value(), {}};
    if (best_key && key.first == best_key->first) {
      key.second = unscaled(image)->to_latex();
      if (best_key->second.empty())
        best_key->second = unscaled(best)->to_latex();
    }
    if (!best_key || key < *best_key) {
      best = std::move(image);
      best_key = std::move(key);
    }
  }
  return best;
}

}  // namespace

bool detail::for_each_symmetrized_term(
    ExprPtr& term, const ResultExpr& result,
    const std::function<void(ExprPtr& part,
                             const IndexPermutationGroup& group)>& f) {
  if (!term.is<Product>()) return false;
  const auto result_indices =
      result.indices() | ranges::to<container::set<Index>>;
  auto is_symmetrizer = [&result_indices](const ExprPtr& factor) {
    if (!factor.is<Tensor>()) return false;
    const auto& tensor = factor.as<Tensor>();
    return (tensor.label() == reserved::symm_label() ||
            tensor.label() == reserved::antisymm_label()) &&
           ranges::all_of(tensor.const_indices(), [&](const Index& idx) {
             return result_indices.contains(idx);
           });
  };
  std::size_t nsums = 0;
  const auto& product = term.as<Product>();
  // the group of the symmetrizers: the product of their groups, which act on
  // disjoint sets of indices
  IndexPermutationGroup group{{{}, 1}};
  container::set<Index> symmetrized_indices;
  std::size_t nsymmetrizers = 0;
  for (const auto& factor : product.factors()) {
    if (is_symmetrizer(factor)) {
      ++nsymmetrizers;
      for (const auto& idx : factor.as<Tensor>().const_indices()) {
        // the symmetrizers do not commute if they share indices
        if (!symmetrized_indices.insert(idx).second) return false;
      }
      IndexPermutationGroup product_group;
      for (const auto& [map1, sign1] : group)
        for (const auto& [map2, sign2] :
             symmetrizer_group(factor.as<Tensor>())) {
          auto map = map1;
          map.insert(map2.begin(), map2.end());
          product_group.emplace_back(std::move(map), sign1 * sign2);
        }
      group = std::move(product_group);
    } else if (factor.is<Sum>())
      ++nsums;
  }
  if (nsymmetrizers == 0) return false;
  if (nsums == 1 && nsymmetrizers + 1 == product.factors().size()) {
    // N.B. mutable iteration resets the memoized hashes
    for (auto& factor : *term) {
      if (!factor.is<Sum>()) continue;
      for (auto& summand : *factor) f(summand, group);
    }
    return true;
  }
  // S * X : f is applied to X, the symmetrizers being kept as they are
  auto rest = std::make_shared<Product>();
  rest->scale(product.scalar());
  auto result_term = std::make_shared<Product>();
  for (const auto& factor : product.factors()) {
    if (is_symmetrizer(factor))
      result_term->append(1, factor, Product::Flatten::No);
    else
      rest->append(1, factor, Product::Flatten::No);
  }
  ExprPtr part = rest;
  f(part, group);
  if (part->is_zero()) {
    term = ex<Constant>(0);
    return true;
  }
  if (part.is<Product>()) {
    result_term->scale(part.as<Product>().scalar());
    for (const auto& factor : part.as<Product>().factors())
      result_term->append(1, factor, Product::Flatten::No);
  } else {
    result_term->append(1, part, Product::Flatten::No);
  }
  term = result_term;
  return true;
}

ResultExpr& canonicalize(ResultExpr& expr, const CanonicalizeOptions& opts) {
  const auto named = named_indices(expr);
  // terms that differ by a permutation of the indices of a symmetrizer that
  // multiplies them are the same term (see detail::for_each_symmetrized_term())
  auto present = [&](ExprPtr& term,
                     const detail::IndexPermutationGroup& group) {
    if (group.size() > 1)
      term = canonicalize_term_for_result(term, group, named, opts);
  };
  auto& rhs = expr.expression();
  if (!rhs) return expr;
  if (rhs->is<Sum>()) {
    for (auto& summand : *rhs)
      detail::for_each_symmetrized_term(summand, expr, present);
  } else {
    detail::for_each_symmetrized_term(rhs, expr, present);
  }
  canonicalize(rhs, named, opts);
  return expr;
}

ResultExpr& canonicalize(ResultExpr&& expr, const CanonicalizeOptions& opts) {
  return canonicalize(expr, opts);
}

struct ExpandVisitor {
  void operator()(ExprPtr& expr) {
    if (Logger::instance().expand)
      std::wcout << "expand_visitor received " << io::latex::to_string(expr)
                 << std::endl;
    // apply expand() iteratively until done
    while (expand(expr)) {
      if (Logger::instance().expand)
        std::wcout << "after 1 round of expansion have "
                   << io::latex::to_string(expr) << std::endl;
    }
    if (Logger::instance().expand)
      std::wcout << "expansion result = " << io::latex::to_string(expr)
                 << std::endl;
    // simplification and canonicalization are to be done by other visitors
  }

  /// expands the first Sum in a Product
  /// @param[in,out] expr (shared_ptr to ) a Product whose first Sum gets
  /// expanded; on return @c expr contains the result
  bool expand_product(ExprPtr& expr) {
    auto& expr_ref = *expr;
    std::shared_ptr<Sum> result;
    const auto nsubexpr = size(expr);
    for (std::size_t i = 0; i != nsubexpr; ++i) {
      if (expr_ref[i]->is<Sum>()) {
        // make template for expr cloning to avoid cloning the Sum we are about
        // to expand
        auto scalar = std::static_pointer_cast<Product>(expr)->scalar();
        auto exprseq_clone_template = container::svector<ExprPtr>(
            ranges::begin(*expr), ranges::end(*expr));
        exprseq_clone_template[i].reset();
        // allocate the result, if not done yet
        if (!result) result = std::make_shared<Sum>();
        ExprPtr subexpr_to_expand = expr_ref[i];
        for (auto& subsubexpr : *subexpr_to_expand) {
          auto exprseq_clone =
              clone(exprseq_clone_template);  // clone the product factors
                                              // without the expanded sum
          exprseq_clone[i] = subsubexpr;      // scavenging summands here
          using std::begin;
          using std::end;
          result->append(
              ex<Product>(scalar, begin(exprseq_clone), end(exprseq_clone)));
        }
        expr =
            std::static_pointer_cast<Expr>(result);  // expanded one Sum, return
        return true;
      }
    }
    return false;
  }

  /// expands a Sum
  bool expand_sum(ExprPtr& expr) {
    auto& expr_ref = *expr;
    std::shared_ptr<Sum>
        result;  // will keep the result if one or more summands is expanded
    const auto nsubexpr = size(expr);
    if (Logger::instance().expand)
      std::wcout << "in expand_sum: expr = " << io::latex::to_string(expr)
                 << std::endl;
    for (std::size_t i = 0; i != nsubexpr; ++i) {
      // if summand is a Product, expand it
      if (expr_ref[i]->is<Product>()) {
        const auto this_term_expanded = expand_product(expr_ref[i]);
        // if this is the first term that was expanded, create a result and copy
        // all preceeding subexpressions into it
        if (!result && this_term_expanded) {
          result = std::make_shared<Sum>();
          for (std::size_t j = 0; j != i; ++j) result->append(expr_ref[j]);
        }
        // if expr != expanded result append current subexpr
        if (result) result->append(expr_ref[i]);
        if (Logger::instance().expand)
          std::wcout << "in expand_sum: after expand_product("
                     << (this_term_expanded ? "true)" : "false)")
                     << " result = "
                     << io::latex::to_string(result ? result : expr)
                     << std::endl;
      }
      // if summand is a Sum, flatten it
      else if (expr_ref[i]->is<Sum>()) {
        // create a result, if not yet created, by copying all preceeding
        // subexpressions into it
        if (!result) {
          result = std::make_shared<Sum>();
          for (std::size_t j = 0; j != i; ++j) result->append(expr_ref[j]);
        }
        if (result) result->append(expr_ref[i]);
        if (Logger::instance().expand)
          std::wcout << "in expand_sum: after flattening Sum summand result = "
                     << io::latex::to_string(result ? result : expr)
                     << std::endl;
      } else {  // nothing to expand? if expanded previously (i.e. result is
                // nonnull) append to result
        if (result) result->append(expr_ref[i]);
      }
    }
    bool expr_changed = false;
    if (result) {  // if any summand was expanded or flattened, copy result into
                   // expr
      expr = std::static_pointer_cast<Expr>(result);
      expr_changed = true;
    }
    if (size(expr) == 1) {  // if sum contains 1 element, raise it
      expr = (*expr)[0];
      expr_changed = true;
    } else if (size(expr) == 0) {  // if sum contains 0 elements, turn to 0
      expr = ex<Constant>(0);
      expr_changed = true;
    }
    return expr_changed;
  }

  // @return true if expanded Product of Sum into Sum of Product
  bool expand(ExprPtr& expr) {
    if (expr->is<Product>()) {
      return expand_product(expr);
    } else if (expr->is<Sum>()) {
      return expand_sum(expr);
    } else
      return false;
  }
};

ExprPtr& expand(ExprPtr& expr) {
  ExpandVisitor expander{};
  expr->visit(expander);
  expander(expr);
  return expr;
}

ExprPtr expand(ExprPtr&& expr) { return expand(expr); }

ResultExpr& expand(ResultExpr& expr) {
  expr.expression() = expand(expr.expression());

  return expr;
}

ResultExpr& expand(ResultExpr&& expr) { return expand(expr); }

ExprPtr& flatten(ExprPtr& expr) {
  auto impl = []<typename E>(std::shared_ptr<E> expr) {
    static_assert(std::is_base_of_v<Expr, E> &&
                  (std::is_same_v<E, Product> || std::is_same_v<E, Sum>));

    bool mutated = false;
    std::shared_ptr<E> flattened_expr;
    for (auto it = expr->begin(); it != expr->end(); ++it) {
      auto& subexpr = *it;
      if (mutated) {
        flattened_expr->append(flatten(subexpr));
        continue;
      }
      auto flattened_subexpr = flatten(subexpr);
      bool rebuild = flattened_subexpr.template is<E>() ||
                     (flattened_subexpr.get() != subexpr.get());
      if (rebuild) {
        mutated = true;
        if constexpr (std::is_same_v<E, Product>) {
          flattened_expr =
              std::make_shared<E>(expr->scalar(), expr->begin(), it);
        } else {
          flattened_expr = std::make_shared<E>(expr->begin(), it);
        }
        flattened_expr->append(flattened_subexpr);
      }
    }
    return mutated ? flattened_expr : expr;
  };

  if (expr.is<Product>()) {
    expr = impl(expr.as_shared_ptr<Product>());
    return expr;
  } else if (expr.is<Sum>()) {
    expr = impl(expr.as_shared_ptr<Sum>());
    return expr;
  } else
    return expr;
}

ExprPtr flatten(ExprPtr&& expr) { return flatten(expr); }

ResultExpr& flatten(ResultExpr& expr) {
  expr.expression() = flatten(expr.expression());

  return expr;
}

ResultExpr& flatten(ResultExpr&& expr) { return flatten(expr); }

struct RapidSimplifyVisitor {
  void operator()(ExprPtr& expr) {
    if (Logger::instance().simplify)
      std::wcout << "rapid_simplify_visitor received "
                 << io::latex::to_string(expr) << std::endl;
    // apply simplify() iteratively until done
    while (simplify(expr)) {
      if (Logger::instance().simplify)
        std::wcout << "after 1 round of simplification have "
                   << io::latex::to_string(expr) << std::endl;
    }
    if (Logger::instance().simplify)
      std::wcout << "simplification result = " << io::latex::to_string(expr)
                 << std::endl;
  }

  /// fuses the Powers of @p product with positive real Constant bases into
  /// simplest radical form (see simplest_radical_form()): the rational part
  /// of their value is multiplied into the scalar, the rest becomes a single
  /// Power `R^(1/D)`. The scalar itself is never fused in, so the rational
  /// part of the value always sits in the scalar and the radicand depends
  /// only on the value; numerically equal Products thus hash and compare
  /// equal, and terms that differ by a rational factor merge when the Sum is
  /// next rebuilt, as by canonicalize() and hence simplify() (rapid_simplify()
  /// leaves summands it has just made identical side by side). Powers
  /// with negative, complex or Variable bases are left alone, and a Product
  /// whose form would exceed the size cap of simplest_radical_form() is left
  /// whole, none of its Powers fused.
  /// @param[in] product a Product
  /// @return the fused Product, or null if @p product is already in this form
  static ExprPtr fuse_rational_powers(const Product& product) {
    auto is_fusible = [](const ExprPtr& factor) {
      if (!factor->is<Power>()) return false;
      const auto& base = factor->as<Power>().base();
      if (!base->is<Constant>()) return false;
      const auto& value = base->as<Constant>().value();
      return value.imag() == 0 && value.real() > 0;
    };
    const auto& factors = product.factors();
    const auto nfusible =
        std::count_if(factors.begin(), factors.end(), is_fusible);
    if (nfusible == 0) return {};

    container::svector<std::pair<rational, rational>> powers;
    for (const auto& factor : factors) {
      if (!is_fusible(factor)) continue;
      const auto& pw = factor->as<Power>();
      powers.emplace_back(pw.base()->as<Constant>().value().real(),
                          pw.exponent());
    }
    const auto form = simplest_radical_form(powers);
    if (!form) return {};

    ExprPtr fused_power;
    if (form->root_order != 1)
      fused_power = ex<Power>(ex<Constant>(form->radicand),
                              rational{1, form->root_order});

    const auto first_fusible =
        std::find_if(factors.begin(), factors.end(), is_fusible);
    // already in this form iff the only fusible factor is R^(1/D) itself
    // and nothing goes into the scalar
    const bool unchanged = form->rational_part == 1 &&
                           nfusible == (fused_power ? 1 : 0) &&
                           (!fused_power || **first_fusible == *fused_power);
    if (unchanged) return {};

    container::svector<ExprPtr> fused_factors;
    fused_factors.reserve(factors.size());
    for (auto it = factors.begin(); it != factors.end(); ++it) {
      if (!is_fusible(*it))
        fused_factors.push_back(*it);
      else if (it == first_fusible && fused_power)
        fused_factors.push_back(fused_power);
    }
    return ex<Product>(product.scalar() * form->rational_part,
                       fused_factors.begin(), fused_factors.end());
  }

  /// fuses the factors of @p product that are a Variable or a Power of one
  /// into one factor per Variable: the exponents are summed, a bare Variable
  /// counting as 1, grouped by label and by the conjugation of the base
  /// combined with that of the Power (as in Power::operator*=). A sum of 0
  /// drops the factor, 1 gives the bare Variable, anything else a Power of
  /// the Variable, with the conjugation on the Variable in either case, so
  /// `(x^2)*` becomes `(x*)^2`. Fused factors keep the position of their
  /// first occurrence.
  /// @param[in] product a Product
  /// @return the fused Product, or null if @p product is already in this form
  static ExprPtr fuse_variable_powers(const Product& product) {
    struct Key {
      std::wstring label;
      bool conjugated;  // of the base combined with that of the Power
      bool operator==(const Key&) const = default;
    };
    auto key_of = [](const ExprPtr& factor) -> std::optional<Key> {
      if (factor->is<Variable>()) {
        const auto& v = factor->as<Variable>();
        return Key{std::wstring(v.label()), v.conjugated()};
      }
      if (factor->is<Power>() && factor->as<Power>().base()->is<Variable>()) {
        const auto& pw = factor->as<Power>();
        const auto& v = pw.base()->as<Variable>();
        return Key{std::wstring(v.label()), v.conjugated() != pw.conjugated()};
      }
      return std::nullopt;
    };

    const auto& factors = product.factors();
    // the exponent of each Key, in order of first occurrence
    container::svector<std::pair<Key, rational>> exponents;
    container::svector<ExprPtr> old_factors;
    std::size_t first_fusible = factors.size();
    for (std::size_t i = 0; i != factors.size(); ++i) {
      const auto& factor = factors[i];
      const auto key = key_of(factor);
      if (!key) continue;
      if (old_factors.empty()) first_fusible = i;
      old_factors.push_back(factor);
      auto it = std::find_if(exponents.begin(), exponents.end(),
                             [&key](const auto& e) { return e.first == *key; });
      if (it == exponents.end())
        it = exponents.insert(exponents.end(), {*key, rational{0}});
      it->second +=
          factor->is<Power>() ? factor->as<Power>().exponent() : rational{1};
    }
    if (old_factors.empty()) return {};

    container::svector<ExprPtr> new_factors;
    for (const auto& [key, exponent] : exponents) {
      if (exponent == 0) continue;
      auto variable = ex<Variable>(key.label);
      if (key.conjugated) variable->adjoint();
      if (exponent == 1)
        new_factors.push_back(std::move(variable));
      else
        new_factors.push_back(ex<Power>(std::move(variable), exponent));
    }
    if (std::equal(old_factors.begin(), old_factors.end(), new_factors.begin(),
                   new_factors.end(),
                   [](const auto& a, const auto& b) { return *a == *b; }))
      return {};

    container::svector<ExprPtr> fused_factors;
    fused_factors.reserve(factors.size());
    for (std::size_t i = 0; i != factors.size(); ++i) {
      if (!key_of(factors[i]))
        fused_factors.push_back(factors[i]);
      else if (i == first_fusible)
        fused_factors.insert(fused_factors.end(), new_factors.begin(),
                             new_factors.end());
    }
    return ex<Product>(product.scalar(), fused_factors.begin(),
                       fused_factors.end());
  }

  /// simplifies a Product by:
  /// - fusing its rational Powers into one, and its Variable factors into
  ///   one per Variable
  /// - flattening subproducts
  /// - factoring in constants
  /// @param[in,out] expr (shared_ptr to ) a Product
  bool simplify_product(ExprPtr& expr) {
    bool expr_changed = false;
    if (auto fused = fuse_rational_powers(expr->as<Product>())) {
      expr = std::move(fused);
      expr_changed = true;
    }
    if (auto fused = fuse_variable_powers(expr->as<Product>())) {
      expr = std::move(fused);
      expr_changed = true;
    }
    auto& expr_ref = *expr;

    // Power::flatten below matters only for the Powers the fusions leave
    // alone, such as those with Variable, negative or complex bases; a
    // fused Power is already in normal form
    // need to rebuild if any factor is a constant or product
    bool need_to_rebuild = false;
    const auto nsubexpr = size(expr);
    for (std::size_t i = 0; i != nsubexpr; ++i) {
      // try to flatten Power: mutates Power -> Constant if possible
      if (expr_ref[i]->is<Power>()) Power::flatten(expr_ref[i]);

      const auto expr_i_is_product = expr_ref[i]->is<Product>();
      const auto expr_i_is_constant = expr_ref[i]->is<Constant>();
      if (expr_i_is_product || expr_i_is_constant) {
        need_to_rebuild = true;
        break;
      }
    }
    if (need_to_rebuild) {
      expr = ex<Product>(expr->as<Product>().scalar(), begin(expr->expr()),
                         end(expr->expr()));
      expr_changed = true;
    }
    const auto expr_size = size(expr);
    auto expr_product = std::static_pointer_cast<Product>(expr);
    if (expr_product->scalar() ==
        0) {  // if scalar = 0, make it 0 (too aggressive?)
      expr = ex<Constant>(0);
      expr_changed = true;
    } else if (expr_size ==
               0) {  // if product reduced to a constant make it a constant
      expr = ex<Constant>(expr_product->scalar());
      expr_changed = true;
    } else if (expr_size == 1 &&
               expr_product->scalar() == 1) {  // if product has 1 term and the
                                               // scalar is 1, lift the factor
      expr = (*expr)[0];
      expr_changed = true;
    }
    return expr_changed;
  }

  /// simplifies a Sum ... generally Sum::{ap,pre}pend simplify automatically,
  /// but the user code may transform sums in a way that the same
  /// simplifications need to be applied here
  bool simplify_sum(ExprPtr& expr) {
    bool mutated = false;
    const Sum& expr_sum = expr->as<Sum>();

    // simplify sums with 0 and 1 arguments
    if (expr_sum.empty()) {
      expr = ex<Constant>(0);
      mutated = true;
    } else if (expr_sum.summands().size() == 1) {
      expr = expr_sum.summands()[0];
      mutated = true;
    } else {  // sums can be simplified if any of its summands are sums or two
      // or more summands are Constants (or have a zero Constant)
      size_t nconst = 0;
      bool need_to_rebuild = false;
      for (auto&& summand : expr_sum.summands()) {
        if (summand->is<Sum>()) {
          need_to_rebuild = true;
          break;
        } else if (summand->is<Constant>()) {
          if (summand->as<Constant>().is_zero()) {
            need_to_rebuild = true;
            break;
          }
          ++nconst;
          if (nconst == 2) {
            need_to_rebuild = true;
            break;
          }
        }
      }
      if (need_to_rebuild) {  // rebuilding will automatically simplify the sum
        auto summands = expr_sum.summands();
        expr = ex<Sum>(begin(summands), end(summands));
        mutated = true;
      }
    }
    return mutated;
  }

  /// simplifies a lone Power as simplify_product() would one inside a
  /// Product: folds it if possible, else brings it to the same form, so that
  /// e.g. `2^(3/2)` becomes `2 2^(1/2)` and merges with it in a Sum
  /// @param[in,out] expr (shared_ptr to ) a Power
  bool simplify_power(ExprPtr& expr) {
    Power::flatten(expr);
    if (!expr->is<Power>()) return true;
    const Product one_factor(ExprPtrList{expr});
    auto fused = fuse_rational_powers(one_factor);
    if (!fused) fused = fuse_variable_powers(one_factor);
    if (!fused) return false;
    expr = std::move(fused);
    return true;
  }

  /// @return true if any simplifications were performed
  bool simplify(ExprPtr& expr) {
    if (expr->is<Product>()) {
      return simplify_product(expr);
    } else if (expr->is<Sum>()) {
      return simplify_sum(expr);
    } else if (expr->is<Power>()) {
      return simplify_power(expr);
    } else
      return false;
  }
};

ExprPtr& rapid_simplify(ExprPtr& expr) {
  RapidSimplifyVisitor simplifier;
  expr->visit(simplifier);
  simplifier(expr);
  return expr;
}

ResultExpr& rapid_simplify(ResultExpr& expr) {
  expr.expression() = rapid_simplify(expr.expression());

  return expr;
}

ResultExpr& rapid_simplify(ResultExpr&& expr) { return rapid_simplify(expr); }

ExprPtr& simplify(ExprPtr& expr, const NamedIndexSpec& named_indices,
                  const CanonicalizeOptions& opts) {
  if (!opts.is_internal_form())
    throw Exception(
        "simplify: presentation settings are not supported; canonicalize "
        "explicitly to choose a presentation");
  // every stage reads the contexts
  const auto contexts_in_effect = pin_default_contexts();
  expand(expr);
  rapid_simplify(expr);
  canonicalize(expr, named_indices, opts);
  rapid_simplify(expr);
  return expr;
}

ExprPtr simplify(ExprPtr&& expr_rv, const NamedIndexSpec& named_indices,
                 const CanonicalizeOptions& opts) {
  auto expr = std::move(expr_rv);
  simplify(expr, named_indices, opts);
  return expr;
}

ResultExpr& simplify(ResultExpr& expr, const CanonicalizeOptions& opts) {
  if (!opts.is_internal_form())
    throw Exception(
        "simplify: presentation settings are not supported; canonicalize "
        "explicitly to choose a presentation");
  // every stage reads the contexts
  const auto contexts_in_effect = pin_default_contexts();
  expand(expr.expression());
  rapid_simplify(expr.expression());
  canonicalize(expr, opts);
  rapid_simplify(expr.expression());
  return expr;
}

ResultExpr& simplify(ResultExpr&& expr, const CanonicalizeOptions& opts) {
  return simplify(expr, opts);
}

ExprPtr& non_canon_simplify(ExprPtr& expr) {
  expand(expr);
  rapid_simplify(expr);
  return expr;
}

ResultExpr& non_canon_simplify(ResultExpr& expr) {
  expand(expr);
  rapid_simplify(expr);
  return expr;
}

}  // namespace sequant
