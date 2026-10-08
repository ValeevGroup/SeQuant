#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/complex.hpp>
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
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/range/primitives.hpp>

#include <algorithm>
#include <iostream>
#include <optional>
#include <string>
#include <utility>

namespace sequant {

std::wstring to_latex_align(const ExprPtr& exprptr, size_t max_lines_per_align,
                            size_t max_terms_per_line) {
  std::wstring result = io::latex::to_string(exprptr);
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

ExprPtr& canonicalize(ExprPtr& expr) {
  // the contexts are read at several points of this call and of the hooks
  const auto contexts_in_effect = pin_default_contexts();
  const auto opts = CanonicalizeOptions::default_options();
  if (expr->is_canonical(opts)) return expr;
  const auto contexts_version = current_contexts_version();
  const auto byproduct = expr->canonicalize(opts);
  if (byproduct && byproduct->is<Constant>()) {
    expr = byproduct * expr;
  }
  // a Sum or Product left with a single term is that term, whose canonical
  // form on its own may differ from its form as a summand or factor
  if (auto term = lone_term(expr)) {
    expr = std::move(term);
    return canonicalize(expr);
  }
  expr->mark_canonical(opts, contexts_version);
  return expr;
}

ExprPtr canonicalize(ExprPtr&& expr_rv) {
  canonicalize(expr_rv);
  return std::move(expr_rv);
}

ResultExpr& canonicalize(ResultExpr& expr) {
  expr.expression() = canonicalize(expr.expression());

  return expr;
}

ResultExpr& canonicalize(ResultExpr&& expr) { return canonicalize(expr); }

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
      if (key.conjugated) variable.as<Variable>().conjugate();
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

namespace {

/// whether the default context's registry contains a complex-field base
/// space or named basis instance; a null registry counts as real
bool default_field_is_complex() {
  auto isr = get_default_context().index_basis_registry();
  if (!isr) return false;
  return ranges::any_of(isr->bases(), [&isr](const auto& entry) {
    const IndexBasis& basis = entry.second;
    return basis.field() == Field::Complex &&
           (basis.has_basis_instance() || isr->is_base(basis.space()));
  });
}

/// scopes copies of the default contexts whose canonicalization options are
/// the current ones with @p setting applied
template <typename Setting>
[[nodiscard]] auto scoped_canonicalize_options_with(Setting setting) {
  return set_scoped_modified_default_context([setting](Context& ctx) {
    ctx.set(ctx.canonicalization_options()
                .value_or(CanonicalizeOptions{})
                .copy_and_set(setting));
  });
}

enum class ConjPairEmission {
  ReIm,       // {s, s*} -> 2 Re(s); {s, -s*} -> 2i Im(s)
  DoubleReal  // {s, s*} -> 2 s (caller asserts the sum's value is real)
};

// Merge Re/Im-wrapped c-number summands related by the conjugate identity
// (Re(x*) == Re(x), Im(x*) == -Im(x)): for a fully contracted c-number
// network the adjoint coincides with the conjugate, so wrappers created by
// different simplify passes -- or emitted by the pair fold itself -- may hold
// conjugate-related inners. Bucket them by a canonical representative and
// accumulate scalars; exact cancellations drop out.
/// a Re/Im-wrapped summand taken apart: `kind` is 1 for Re, 2 for Im and 0
/// for a summand that is not a wrapper (or a scaled wrapper)
struct WrapInfo {
  int kind = 0;
  Constant::scalar_type scalar = 1;
  ExprPtr inner;
};

/// @return the wrapper @p sm is, with a real scalar of the wrapped Product
///         hoisted into `scalar` (Re and Im are real-linear, so the two
///         spellings then share one representative)
WrapInfo classify_wrapped(ExprPtr const& sm) {
  auto info = [&sm]() -> WrapInfo {
    if (sm->is<RealPart>()) return {1, 1, sm->as<RealPart>().inner()};
    if (sm->is<ImagPart>()) return {2, 1, sm->as<ImagPart>().inner()};
    if (sm->is<Product>()) {
      auto const& p = sm->as<Product>();
      if (p.factors().size() == 1) {
        auto const& f = p.factor(0);
        if (f->is<RealPart>())
          return {1, p.scalar(), f->as<RealPart>().inner()};
        if (f->is<ImagPart>())
          return {2, p.scalar(), f->as<ImagPart>().inner()};
      }
    }
    return {};
  }();
  if (info.kind != 0) {
    if (info.inner->is<Product>()) {
      auto const& p = info.inner->as<Product>();
      auto const c = p.scalar();
      if (c.imag() == 0 && c.real() != 1) {
        info.inner = detail::strip_scalar(p);
        info.scalar *= c;
      }
    }
  }
  return info;
}

/// @return whether @p sm is scalar-valued c-number content: what the pair
///         fold pairs (see fold_conjugate_pairs_impl)
bool is_closed_cnumber(ExprPtr const& sm) {
  if (!sm->is_cnumber()) return false;
  auto const ext = get_unique_indices(sm);
  return ext.bra.empty() && ext.ket.empty() && ext.aux.empty();
}

template <typename SummandRange>
container::svector<ExprPtr> merge_wrapped_summands(
    SummandRange const& in,
    std::function<ExprPtr(ExprPtr const&)> const& conjugate_op) {
  container::svector<ExprPtr> out;
  struct Bucket {
    int kind;
    ExprPtr rep;
    Constant::scalar_type acc = 0;
  };
  container::map<std::size_t, container::svector<Bucket>> buckets;
  container::svector<std::pair<std::size_t, std::size_t>> order;
  for (auto const& sm : in) {
    auto wi = classify_wrapped(sm);
    if (wi.kind == 0 || !wi.inner->is_cnumber()) {
      out.push_back(sm->clone());
      continue;
    }
    auto ci = canonicalize(wi.inner->clone());
    ExprPtr conj_inner =
        conjugate_op ? conjugate_op(wi.inner) : sequant::adjoint(wi.inner);
    auto cc = canonicalize(conj_inner->clone());
    // the representative is the smaller of the two canonical spellings under
    // Expr::operator<, the criterion fold_conjugate_pairs_impl's
    // `representative` uses: a min over the pair, so it does not depend on the
    // order the sum was written in. When the adjoint hands up a -1 the
    // representative may carry that scalar inside `rep`, and emission spells
    // the summand as e.g. -2 Re[d]; that is value-correct (an anti-Hermitian d
    // with d⁺ = d꙳ has Re d = 0), and hoisting a real scalar out of `ci` and
    // `cc` here is the change to make if a cleaner spelling is wanted
    bool use_conj = *cc < *ci;
    ExprPtr rep = use_conj ? cc : ci;
    auto sc = wi.scalar;
    if (use_conj && wi.kind == 2) sc = -sc;  // Im(x*) = -Im(x)
    auto key = rep->hash_value();
    hash::combine(key, static_cast<std::size_t>(wi.kind));
    auto& vec = buckets[key];
    bool merged = false;
    for (std::size_t b = 0; b != vec.size(); ++b)
      if (vec[b].kind == wi.kind && *vec[b].rep == *rep) {
        vec[b].acc += sc;
        merged = true;
        break;
      }
    if (!merged) {
      vec.push_back(Bucket{wi.kind, rep, sc});
      order.emplace_back(key, vec.size() - 1);
    }
  }
  for (auto const& [key, idx] : order) {
    auto const& b = buckets[key][idx];
    if (b.acc == Constant::scalar_type(0)) continue;
    auto wrapped = b.kind == 1 ? real_part(b.rep->clone())
                               : imaginary_part(b.rep->clone());
    out.push_back(ex<Constant>(b.acc) * wrapped);
  }
  return out;
}

ExprPtr fold_conjugate_pairs_impl(
    ExprPtr const& expr, std::function<ExprPtr(ExprPtr const&)> conjugate_op,
    ConjPairEmission emission) {
  if (!expr || !expr->is<Sum>()) return expr;
  // nothing to fold unless some summand is a wrapper to merge or a
  // scalar-valued summand to pair; a tensor-valued sum (a residual) passes
  // through untouched, without the copies the passes below make
  if (ranges::none_of(expr->as<Sum>().summands(), [](ExprPtr const& sm) {
        return classify_wrapped(sm).kind != 0 || is_closed_cnumber(sm);
      }))
    return expr;
  // cross-summand identity requires meaningful named (external) labels,
  // same reasoning as Sum::canonicalize_impl
  const auto named_labels_matter = scoped_canonicalize_options_with(
      CanonicalizeOptions::IgnoreNamedIndexLabel::No);

  // Pre-merge Re/Im-wrapped summands related by the conjugate identity
  // (Re(x*) == Re(x), Im(x*) == -Im(x)). Wrapped summands are
  // self-conjugate, so the pair fold below leaves them untouched -- but
  // wrappers created by separate simplify passes may hold conjugate-related
  // inners (for a fully contracted c-number network the adjoint coincides
  // with the conjugate), e.g. +c Re(X) and -c Re(X^+) must cancel.
  auto summands_v =
      merge_wrapped_summands(expr->as<Sum>().summands(), conjugate_op);
  auto const& summands = summands_v;
  const std::size_t n = summands.size();
  // the fold applies to scalar-valued summands only. Re/Im of
  // operator-valued content is out of scope here (the operator analogue --
  // anti-Hermitian splitting -- comes with the time-reversal work), and an
  // operator string's adjoint reverses the operators, which is not this
  // fold's elementwise conjugation. A tensor-valued summand (one with
  // external indices) stays out for two reasons: the pairing below compares
  // canonical forms, which identify a summand with its conjugate as a value
  // only when there are no externals to line up (the adjoint of R{a;i} is
  // R{i;a}, a different tensor, so 2 Re would not be the sum's value), and
  // Re/Im are evaluated and exported for scalar results only.
  std::vector<bool> eligible(n);
  for (std::size_t i = 0; i != n; ++i)
    eligible[i] = is_closed_cnumber(summands[i]);
  std::vector<ExprPtr> canon(n), canon_conj(n), canon_negconj(n);
  container::map<std::size_t, container::svector<std::size_t>> buckets;
  for (std::size_t i = 0; i != n; ++i) {
    if (!eligible[i]) continue;
    canon[i] = canonicalize(summands[i]->clone());
    buckets[canon[i]->hash_value()].push_back(i);
    ExprPtr conj;
    if (conjugate_op) {
      conj = conjugate_op(summands[i]);
    } else {
      conj = sequant::adjoint(summands[i]);
    }
    canon_conj[i] = canonicalize(conj->clone());
    if (emission == ConjPairEmission::ReIm)
      canon_negconj[i] = canonicalize(ex<Constant>(-1) * conj);
  }

  // greedy first-match pairing via the hash buckets, verified structurally
  std::vector<bool> consumed(n, false);
  std::vector<int8_t> fold(n, 0);  // 0 = keep as-is, +1 = 2Re, -1 = 2iIm
  std::vector<std::size_t> partner(n, n);  // the other member of the pair
  auto probe = [&](std::size_t i, ExprPtr const& key) -> std::size_t {
    auto it = buckets.find(key->hash_value());
    if (it == buckets.end()) return n;
    for (std::size_t j : it->second)
      if (j != i && !consumed[j] && !fold[j] && *canon[j] == *key) return j;
    return n;
  };
  for (std::size_t i = 0; i != n; ++i) {
    if (!eligible[i] || consumed[i] || fold[i]) continue;
    if (canon[i]->hash_value() == canon_conj[i]->hash_value() &&
        *canon[i] == *canon_conj[i])
      continue;  // self-conjugate (manifestly real): leave untouched
    if (auto j = probe(i, canon_conj[i]); j != n) {
      fold[i] = +1;
      partner[i] = j;
      consumed[j] = true;
      continue;
    }
    if (emission == ConjPairEmission::ReIm) {
      if (auto j = probe(i, canon_negconj[i]); j != n) {
        fold[i] = -1;
        partner[i] = j;
        consumed[j] = true;
      }
    }
  }

  // the representative a folded pair's Re/Im wrapper carries: the smaller of
  // the pair's two canonical forms. Both denote the same wrapped value
  // (`Re s = Re s*` and `Im s = Im(-s*)`), so picking by the expression order
  // is what makes the fold independent of the order the pair was written in
  auto representative = [&](std::size_t i) {
    SEQUANT_ASSERT(partner[i] != n);
    auto const& a = canon[i];
    auto const& b = canon[partner[i]];
    return (*b < *a ? b : a)->clone();
  };

  auto result = std::make_shared<Sum>();
  for (std::size_t i = 0; i != n; ++i) {
    if (consumed[i]) continue;
    if (fold[i] == +1) {
      result->append(emission == ConjPairEmission::ReIm
                         ? ex<Constant>(2) * real_part(representative(i))
                         : ex<Constant>(2) * summands[i]->clone());
    } else if (fold[i] == -1) {
      // s + (-s*) = s - s* = 2i Im(s)
      result->append(ex<Constant>(Constant::scalar_type(0, 2)) *
                     imaginary_part(representative(i)));
    } else {
      result->append(summands[i]->clone());
    }
  }
  auto merged_out = merge_wrapped_summands(result->summands(), conjugate_op);
  auto result2 = std::make_shared<Sum>();
  for (auto& sm : merged_out) result2->append(std::move(sm));
  if (result2->summands().empty()) return ex<Constant>(0);
  if (result2->summands().size() == 1) return result2->summands().front();
  return std::static_pointer_cast<Expr>(result2);
}

}  // namespace

ExprPtr fold_conjugate_pairs(
    ExprPtr const& expr, std::function<ExprPtr(ExprPtr const&)> conjugate_op) {
  return fold_conjugate_pairs_impl(expr, std::move(conjugate_op),
                                   ConjPairEmission::ReIm);
}

ExprPtr fold_conjugate_pairs_of_real_sum(
    ExprPtr const& expr, std::function<ExprPtr(ExprPtr const&)> conjugate_op) {
  return fold_conjugate_pairs_impl(expr, std::move(conjugate_op),
                                   ConjPairEmission::DoubleReal);
}

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

ExprPtr& simplify(ExprPtr& expr) {
  // every stage reads the contexts
  const auto contexts_in_effect = pin_default_contexts();
  expand(expr);
  rapid_simplify(expr);
  canonicalize(expr);
  // complex field: fold conjugate-related summand pairs exactly
  // (A + A* -> 2 Re(A)); in a real field conjugation is trivial and plain
  // canonicalization already merges such pairs. The fold applies only to
  // fully c-number content: an expression still carrying operators is an
  // intermediate of a derivation (Wick consumes it next, and the Wick
  // engine does not ingest RealPart/ImagPart wrappers). Within that, the
  // fold itself pairs scalar-valued summands only, so a tensor-valued sum
  // (a residual) passes through unfolded. The fold runs after
  // the canonicalize pass (so trivially-cancelling spellings are already
  // merged); its own pre-pass canonicalizes existing wrappers' inners, so
  // conjugate-related wrappers from earlier simplify passes merge exactly.
  if (CanonicalizeOptions::default_options().fold_conjugate_pairs ==
          CanonicalizeOptions::FoldConjugatePairs::Yes &&
      default_field_is_complex() && expr->is_cnumber()) {
    expr = fold_conjugate_pairs(expr);
  }
  rapid_simplify(expr);
  return expr;
}

ExprPtr simplify(ExprPtr&& expr_rv) {
  auto expr = std::move(expr_rv);
  simplify(expr);
  return expr;
}

ResultExpr& simplify(ResultExpr& expr) {
  expr.expression() = simplify(expr.expression());

  return expr;
}

ResultExpr& simplify(ResultExpr&& expr) { return simplify(expr); }

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

bool is_hermitian_network(ExprPtr const& expr) {
  SEQUANT_ASSERT(expr && expr->is_cnumber());
  // cross-expression identity requires meaningful named (external) labels
  const auto named_labels_matter = scoped_canonicalize_options_with(
      CanonicalizeOptions::IgnoreNamedIndexLabel::No);
  auto lhs = canonicalize(expr->clone());
  auto rhs = canonicalize(sequant::adjoint(expr));
  return lhs->hash_value() == rhs->hash_value() && *lhs == *rhs;
}

namespace {

/// whether @p idx and every index of its proto-index closure lie in a
/// `K`-closed space
bool is_kclosed(const Index& idx) {
  if (idx.basis().field() != Field::Real) return false;
  return ranges::all_of(idx.proto_indices(),
                        [](const Index& p) { return is_kclosed(p); });
}

/// whether every index @p op acts on lies in a `K`-closed space, the
/// condition for `K O K⁻¹` to be the same operator string
template <Statistics S>
bool acts_on_kclosed_spaces(const NormalOperator<S>& op) {
  return ranges::all_of(op,
                        [](const auto& o) { return is_kclosed(o.index()); });
}

/// @overload for a sequence of normal operators
template <Statistics S>
bool acts_on_kclosed_spaces(const NormalOperatorSequence<S>& seq) {
  return ranges::all_of(
      seq, [](const auto& op) { return acts_on_kclosed_spaces(op); });
}

/// @throw Exception if an operator in @p expr acts on an index space that
///        `K` does not close, so that `K E K⁻¹` has no spelling in the
///        indices at hand
void assert_kclosed_operators(const ExprPtr& expr) {
  std::as_const(*expr).visit(
      [](const ExprPtr& atom) {
        const bool closed =
            atom->is<FNOperator>()
                ? acts_on_kclosed_spaces(atom->as<FNOperator>())
            : atom->is<BNOperator>()
                ? acts_on_kclosed_spaces(atom->as<BNOperator>())
            : atom->is<FNOperatorSeq>()
                ? acts_on_kclosed_spaces(atom->as<FNOperatorSeq>())
            : atom->is<BNOperatorSeq>()
                ? acts_on_kclosed_spaces(atom->as<BNOperatorSeq>())
                : true;
        if (!closed)
          throw Exception(
              "sequant::kconjugate: an operator over a complex basis has no "
              "K-closed index space");
      },
      /*atoms_only=*/true);
}

}  // namespace

ExprPtr kconjugate(const ExprPtr& expr) {
  SEQUANT_ASSERT(expr);
  // K acts on the basis the operators are written in, so an operator string
  // is reproduced only where every index space is K-closed
  if (!expr->is_cnumber()) assert_kclosed_operators(expr);
  auto result = expr->clone();
  const auto sign = result->kconjugate();
  if (sign == 1) return result;
  return ex<Product>(sign, ExprPtrList{std::move(result)});
}

ExprPtr conjugate(const ExprPtr& expr) {
  SEQUANT_ASSERT(expr);
  // the complex conjugate of a value: for a matrix element
  // conj <p|O|q> = <q|O⁺|p>, so on c-number content this is the adjoint (a
  // Product's factors commute, so the adjoint's reversal is not observable)
  if (!expr->is_cnumber())
    throw Exception(
        "sequant::conjugate: an operator has no value to conjugate; "
        "sequant::kconjugate is the conjugation of an operator");
  return sequant::adjoint(expr);
}

}  // namespace sequant
