#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_operators.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/options.hpp>

#include <range/v3/range/primitives.hpp>

#include <cassert>
#include <iostream>
#include <string>
#include <utility>

namespace sequant {

std::wstring to_latex(const ExprPtr& exprptr) { return exprptr->to_latex(); }

std::wstring to_latex_align(const ExprPtr& exprptr, size_t max_lines_per_align,
                            size_t max_terms_per_line) {
  std::wstring result = to_latex(exprptr);
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
      assert(pos != std::wstring::npos);
      result.insert(at, str);
      const auto str_nchar = std::size(str) - 1;  // neglect end-of-string
      pos += str_nchar;
      if (plus_pos != std::wstring::npos) plus_pos += str_nchar;
      if (minus_pos != std::wstring::npos) minus_pos += str_nchar;
      if (pos != plus_pos) assert(plus_pos == result.find(L" + ", plus_pos));
      if (pos != minus_pos) assert(minus_pos == result.find(L" - ", minus_pos));
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

std::wstring to_wolfram(const ExprPtr& exprptr) {
  return exprptr->to_wolfram();
}

std::size_t size(const Expr& expr) { return ranges::size(expr); }

std::size_t size(const ExprPtr& exprptr) {
  if (exprptr) {
    return size(*exprptr);
  }

  return 0;
}

ExprPtr& canonicalize(ExprPtr& expr, CanonicalizeOptions opts) {
  const auto byproduct = expr->canonicalize(opts);
  if (byproduct && byproduct->is<Constant>()) {
    expr = byproduct * expr;
  }
  return expr;
}

ExprPtr canonicalize(ExprPtr&& expr_rv, CanonicalizeOptions opts) {
  const auto byproduct = expr_rv->canonicalize(opts);
  if (byproduct && byproduct->is<Constant>()) {
    expr_rv = byproduct * expr_rv;
  }
  return std::move(expr_rv);
}

ResultExpr& canonicalize(ResultExpr& expr, CanonicalizeOptions opts) {
  expr.expression() = canonicalize(expr.expression(), std::move(opts));

  return expr;
}

ResultExpr& canonicalize(ResultExpr&& expr, CanonicalizeOptions opts) {
  return canonicalize(expr, std::move(opts));
}

struct ExpandVisitor {
  void operator()(ExprPtr& expr) {
    if (Logger::instance().expand)
      std::wcout << "expand_visitor received " << to_latex(expr) << std::endl;
    // apply expand() iteratively until done
    while (expand(expr)) {
      if (Logger::instance().expand)
        std::wcout << "after 1 round of expansion have " << to_latex(expr)
                   << std::endl;
    }
    if (Logger::instance().expand)
      std::wcout << "expansion result = " << to_latex(expr) << std::endl;
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
      std::wcout << "in expand_sum: expr = " << to_latex(expr) << std::endl;
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
                     << " result = " << to_latex(result ? result : expr)
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
                     << to_latex(result ? result : expr) << std::endl;
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
  SimplifyOptions opts;

  RapidSimplifyVisitor(SimplifyOptions opts) : opts(std::move(opts)) {
    opts.method = CanonicalizationMethod::Rapid;
  }

  void operator()(ExprPtr& expr) {
    if (Logger::instance().simplify)
      std::wcout << "rapid_simplify_visitor received " << to_latex(expr)
                 << std::endl;
    // apply simplify() iteratively until done
    while (simplify(expr, opts)) {
      if (Logger::instance().simplify)
        std::wcout << "after 1 round of simplification have " << to_latex(expr)
                   << std::endl;
    }
    if (Logger::instance().simplify)
      std::wcout << "simplification result = " << to_latex(expr) << std::endl;
  }

  /// simplifies a Product by:
  /// - flattening subproducts
  /// - factoring in constants
  /// @param[in,out] expr (shared_ptr to ) a Product
  bool simplify_product(ExprPtr& expr,
                        SimplifyOptions = SimplifyOptions::default_options()) {
    auto& expr_ref = *expr;

    // need to rebuild if any factor is a constant or product
    bool need_to_rebuild = false;
    const auto nsubexpr = size(expr);
    for (std::size_t i = 0; i != nsubexpr; ++i) {
      const auto expr_i_is_product = expr_ref[i]->is<Product>();
      const auto expr_i_is_constant = expr_ref[i]->is<Constant>();
      if (expr_i_is_product || expr_i_is_constant) {
        need_to_rebuild = true;
        break;
      }
    }
    bool expr_changed = false;
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
  bool simplify_sum(ExprPtr& expr,
                    SimplifyOptions = SimplifyOptions::default_options()) {
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

  // @return true if any simplifications were performed
  bool simplify(ExprPtr& expr,
                SimplifyOptions opts = SimplifyOptions::default_options()) {
    if (expr->is<Product>()) {
      return simplify_product(expr, opts);
    } else if (expr->is<Sum>()) {
      return simplify_sum(expr, opts);
    } else
      return false;
  }
};

ExprPtr& rapid_simplify(ExprPtr& expr, SimplifyOptions opts) {
  RapidSimplifyVisitor simplifier{opts};
  expr->visit(simplifier);
  simplifier(expr);
  return expr;
}

ResultExpr& rapid_simplify(ResultExpr& expr, SimplifyOptions opts) {
  expr.expression() = rapid_simplify(expr.expression(), std::move(opts));

  return expr;
}

ResultExpr& rapid_simplify(ResultExpr&& expr, SimplifyOptions opts) {
  return rapid_simplify(expr, std::move(opts));
}

ExprPtr& simplify(ExprPtr& expr, SimplifyOptions opts) {
  expand(expr);
  rapid_simplify(expr, opts);
  canonicalize(expr, opts);
  rapid_simplify(expr, opts);
  return expr;
}

ExprPtr simplify(ExprPtr&& expr_rv, SimplifyOptions opts) {
  auto expr = std::move(expr_rv);
  simplify(expr, opts);
  return expr;
}

ResultExpr& simplify(ResultExpr& expr, SimplifyOptions opts) {
  expr.expression() = simplify(expr.expression(), std::move(opts));

  return expr;
}

ResultExpr& simplify(ResultExpr&& expr, SimplifyOptions opts) {
  return simplify(expr, std::move(opts));
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
