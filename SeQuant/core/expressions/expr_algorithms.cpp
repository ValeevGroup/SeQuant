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
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/range/primitives.hpp>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

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
  SimplifyOptions opts;

  RapidSimplifyVisitor(SimplifyOptions opts) : opts(std::move(opts)) {
    opts.method = CanonicalizationMethod::Rapid;
  }

  void operator()(ExprPtr& expr) {
    if (Logger::instance().simplify)
      std::wcout << "rapid_simplify_visitor received "
                 << io::latex::to_string(expr) << std::endl;
    // apply simplify() iteratively until done
    while (simplify(expr, opts)) {
      if (Logger::instance().simplify)
        std::wcout << "after 1 round of simplification have "
                   << io::latex::to_string(expr) << std::endl;
    }
    if (Logger::instance().simplify)
      std::wcout << "simplification result = " << io::latex::to_string(expr)
                 << std::endl;
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
      // try to flatten Power: mutates Power -> Constant if possible
      if (expr_ref[i]->is<Power>()) Power::flatten(expr_ref[i]);

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

namespace {

/// whether the default context's registry contains a complex-field base
/// space; a null registry counts as real
bool default_field_is_complex() {
  auto isr = get_default_context().index_space_registry();
  if (!isr) return false;
  return ranges::any_of(isr->base_spaces(), [](const IndexSpace& s) {
    return s.field() == Field::Complex;
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
template <typename SummandRange>
container::svector<ExprPtr> merge_wrapped_summands(
    SummandRange const& in, CanonicalizeOptions const& opts,
    std::function<ExprPtr(ExprPtr const&)> const& conjugate_op) {
  struct WrapInfo {
    int kind = 0;
    Constant::scalar_type scalar = 1;
    ExprPtr inner;
  };
  // Re and Im are real-linear (Re(c X) = c Re(X), Im(c X) = c Im(X) for a
  // real c), so a real scalar belongs with the summand's scalar rather than
  // inside the wrapper: the two spellings then share one representative
  auto hoist_real_scalar = [](ExprPtr& e) -> Constant::scalar_type {
    if (!e->is<Product>()) return 1;
    auto const& p = e->as<Product>();
    auto const c = p.scalar();
    if (c.imag() != 0 || c.real() == 1) return 1;
    e = detail::strip_scalar(p);
    return c;
  };
  auto classify = [&hoist_real_scalar](ExprPtr const& sm) -> WrapInfo {
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
    if (info.kind != 0) info.scalar *= hoist_real_scalar(info.inner);
    return info;
  };
  container::svector<ExprPtr> out;
  struct Bucket {
    int kind;
    ExprPtr rep;
    Constant::scalar_type acc = 0;
  };
  container::map<std::size_t, container::svector<Bucket>> buckets;
  container::svector<std::pair<std::size_t, std::size_t>> order;
  for (auto const& sm : in) {
    auto wi = classify(sm);
    if (wi.kind == 0 || !wi.inner->is_cnumber()) {
      out.push_back(sm->clone());
      continue;
    }
    auto ci = canonicalize(wi.inner->clone(), opts);
    ExprPtr conj_inner =
        conjugate_op ? conjugate_op(wi.inner) : sequant::adjoint(wi.inner);
    auto cc = canonicalize(conj_inner->clone(), opts);
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
    ExprPtr const& expr, CanonicalizeOptions opts,
    std::function<ExprPtr(ExprPtr const&)> conjugate_op,
    ConjPairEmission emission) {
  if (!expr || !expr->is<Sum>()) return expr;
  // cross-summand identity requires meaningful named (external) labels,
  // same reasoning as Sum::canonicalize_impl
  opts = opts.copy_and_set(CanonicalizeOptions::IgnoreNamedIndexLabel::No);

  // Pre-merge Re/Im-wrapped summands related by the conjugate identity
  // (Re(x*) == Re(x), Im(x*) == -Im(x)). Wrapped summands are
  // self-conjugate, so the pair fold below leaves them untouched -- but
  // wrappers created by separate simplify passes may hold conjugate-related
  // inners (for a fully contracted c-number network the adjoint coincides
  // with the conjugate), e.g. +c Re(X) and -c Re(X^+) must cancel.
  auto summands_v =
      merge_wrapped_summands(expr->as<Sum>().summands(), opts, conjugate_op);
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
  for (std::size_t i = 0; i != n; ++i) {
    eligible[i] = summands[i]->is_cnumber();
    if (eligible[i]) {
      auto const ext = get_unique_indices(summands[i]);
      eligible[i] = ext.bra.empty() && ext.ket.empty() && ext.aux.empty();
    }
  }
  std::vector<ExprPtr> canon(n), canon_conj(n), canon_negconj(n);
  container::map<std::size_t, container::svector<std::size_t>> buckets;
  for (std::size_t i = 0; i != n; ++i) {
    if (!eligible[i]) continue;
    canon[i] = canonicalize(summands[i]->clone(), opts);
    buckets[canon[i]->hash_value()].push_back(i);
    ExprPtr conj;
    if (conjugate_op) {
      conj = conjugate_op(summands[i]);
    } else {
      conj = sequant::adjoint(summands[i]);
    }
    canon_conj[i] = canonicalize(conj->clone(), opts);
    if (emission == ConjPairEmission::ReIm)
      canon_negconj[i] = canonicalize(ex<Constant>(-1) * conj, opts);
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
  auto merged_out =
      merge_wrapped_summands(result->summands(), opts, conjugate_op);
  auto result2 = std::make_shared<Sum>();
  for (auto& sm : merged_out) result2->append(std::move(sm));
  if (result2->summands().empty()) return ex<Constant>(0);
  if (result2->summands().size() == 1) return result2->summands().front();
  return std::static_pointer_cast<Expr>(result2);
}

}  // namespace

ExprPtr fold_conjugate_pairs(
    ExprPtr const& expr, CanonicalizeOptions opts,
    std::function<ExprPtr(ExprPtr const&)> conjugate_op) {
  return fold_conjugate_pairs_impl(
      expr, std::move(opts), std::move(conjugate_op), ConjPairEmission::ReIm);
}

ExprPtr fold_conjugate_pairs_of_real_sum(
    ExprPtr const& expr, CanonicalizeOptions opts,
    std::function<ExprPtr(ExprPtr const&)> conjugate_op) {
  return fold_conjugate_pairs_impl(expr, std::move(opts),
                                   std::move(conjugate_op),
                                   ConjPairEmission::DoubleReal);
}

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
  if (opts.fold_conjugate_pairs == SimplifyOptions::FoldConjugatePairs::Yes &&
      default_field_is_complex() && expr->is_cnumber()) {
    expr = fold_conjugate_pairs(expr, opts);
  }
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

bool is_hermitian_network(ExprPtr const& expr, CanonicalizeOptions opts) {
  SEQUANT_ASSERT(expr && expr->is_cnumber());
  // cross-expression identity requires meaningful named (external) labels
  opts = opts.copy_and_set(CanonicalizeOptions::IgnoreNamedIndexLabel::No);
  auto lhs = canonicalize(expr->clone(), opts);
  auto rhs = canonicalize(sequant::adjoint(expr), opts);
  return lhs->hash_value() == rhs->hash_value() && *lhs == *rhs;
}

namespace {

/// whether @p idx and every index of its proto-index closure lie in a
/// `K`-closed space
bool is_kclosed(const Index& idx) {
  if (idx.space().field() != Field::Real) return false;
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
