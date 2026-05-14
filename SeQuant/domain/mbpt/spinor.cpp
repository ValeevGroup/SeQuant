//
// Spinor decomposition passes for relativistic 2-component theories.
// See spinor.hpp for design overview.
//

#include <SeQuant/domain/mbpt/spinor.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/index_space_registry.hpp>

#include <algorithm>
#include <stdexcept>

namespace sequant::mbpt {

namespace detail {

/// Construct a new Index in the same IndexSpace type as @p idx, with
/// the @p k Kramers QN bit set (replacing any existing Kramers state)
/// and the label decorated with the matching ⇑/⇓ glyph. Spin bits
/// (alpha/beta) are NOT touched here — Spin and Kramers are mutually
/// exclusive on a single index, but the per-index check is left to
/// the caller (a typical workflow only ever sets one of them).
///
/// Mirrors `detail::make_index_with_spincase` in `spin.cpp` so that
/// SeQuant's IndexSpace registry is reused: if a space with the
/// resulting label and qns already exists, we adopt it.
Index make_index_with_kramers(const Index& idx, mbpt::Kramers k) {
  // sanity check: at most one Kramers annotation already
  SEQUANT_ASSERT(!(idx.label().find(L'⇑') != std::wstring::npos &&
                   idx.label().find(L'⇓') != std::wstring::npos));

  // Strip any existing Kramers bits, then OR in the requested state.
  auto qns = mbpt::kramers_annotation_remove(idx.space().qns()).unIon(k);

  IndexSpace space;
  const auto label =
      mbpt::kramers_annotation_replace(idx.space().base_key(), k);
  if (auto isr = get_default_context().index_space_registry()) {
    auto* space_ptr = isr->retrieve_ptr(label);
    if (space_ptr && space_ptr->type() == idx.space().type() &&
        space_ptr->qns() == qns) {
      space = *space_ptr;
    }
  }
  if (!space) {
    space = IndexSpace{label, idx.space().type(), qns,
                       // assume size does not depend on Kramers state
                       idx.space().approximate_size()};
  }
  auto protoindices = idx.proto_indices();
  for (auto& pidx : protoindices) pidx = make_index_with_kramers(pidx, k);
  return Index{space, idx.ordinal(), protoindices};
}

}  // namespace detail

Index make_kramers_up(const Index& idx) {
  return detail::make_index_with_kramers(idx, mbpt::Kramers::up);
}

Index make_kramers_dn(const Index& idx) {
  return detail::make_index_with_kramers(idx, mbpt::Kramers::down);
}

Index make_kramers_free(const Index& idx) {
  return detail::make_index_with_kramers(idx, mbpt::Kramers::any);
}

// =====================================================================
// Conjugation infrastructure (groundwork for Kramers tracer + quaternion
// decomposer).
// =====================================================================

std::wstring toggle_conj_suffix(std::wstring_view label) {
  if (has_conj_suffix(label)) return std::wstring{label.substr(0, label.size() - 1)};
  return std::wstring{label} + L'*';
}

namespace {

/// Apply conjugation to a single Tensor: rebuild it with toggled label
/// and identical bra/ket/aux/symmetry/etc.
ExprPtr conjugate_tensor(const Tensor& t) {
  auto new_label = toggle_conj_suffix(t.label());
  // Tensor copy ctor + label edit isn't directly exposed; reconstruct.
  return ex<Tensor>(
      new_label,
      bra(container::svector<Index>{t.bra().begin(), t.bra().end()}),
      ket(container::svector<Index>{t.ket().begin(), t.ket().end()}),
      aux(container::svector<Index>{t.aux().begin(), t.aux().end()}),
      t.symmetry(), t.braket_symmetry(), t.column_symmetry());
}

}  // namespace

ExprPtr conjugate(const ExprPtr& expr) {
  if (expr->is<Tensor>()) {
    return conjugate_tensor(expr->as<Tensor>());
  }
  if (expr->is<Constant>()) {
    auto v = expr->as<Constant>().value();
    return ex<Constant>(sequant::conj(v));
  }
  if (expr->is<Variable>()) {
    throw std::runtime_error(
        "conjugate: Variable conjugation not yet supported");
  }
  if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();
    auto p = std::make_shared<Product>();
    p->scale(sequant::conj(prod.scalar()));
    for (auto&& f : prod) p->append(1, conjugate(f));
    return p;
  }
  if (expr->is<Sum>()) {
    auto s = std::make_shared<Sum>();
    for (auto&& summand : *expr) s->append(conjugate(summand));
    return s;
  }
  throw std::runtime_error("conjugate: unsupported Expr type");
}

ExprPtr real_part(const ExprPtr& expr) {
  // Re(z) = (z + z*) / 2 — represented as a literal Sum with a 1/2 prefactor.
  auto s = std::make_shared<Sum>();
  s->append(expr);
  s->append(conjugate(expr));
  return ex<Constant>(rational{1, 2}) * ExprPtr{s};
}

ExprPtr imaginary_part(const ExprPtr& expr) {
  // Im(z) = (z - z*) / (2i) = -i (z - z*) / 2.
  auto s = std::make_shared<Sum>();
  s->append(expr);
  // -1 * conj(expr): scale conj by -1 via a wrapping Product.
  auto neg_conj = ex<Constant>(-1) * conjugate(expr);
  s->append(neg_conj);
  // multiply by -i / 2:
  return ex<Constant>(rational{1, 2}) *
         ex<Constant>(Complex<rational>{0, -1}) * ExprPtr{s};
}

ExprPtr append_kramers(
    const ExprPtr& expr,
    const container::map<Index, Index>& index_replacements) {
  auto add_to_tensor = [&](const Tensor& tensor) {
    auto t = std::make_shared<Tensor>(tensor);
    t->transform_indices(index_replacements);
    return t;
  };
  auto add_to_product = [&](const Product& product) {
    auto p = std::make_shared<Product>();
    p->scale(product.scalar());
    for (auto&& term : product) {
      if (term->is<Tensor>()) {
        p->append(1, add_to_tensor(term->as<Tensor>()));
      } else if (term->is<Constant>() || term->is<Variable>()) {
        p->append(1, term);
      } else {
        throw std::runtime_error(
            "append_kramers: unsupported Expr type in product: " +
            term->type_name());
      }
    }
    return p;
  };

  if (expr->is<Tensor>()) return add_to_tensor(expr->as<Tensor>());
  if (expr->is<Product>()) return add_to_product(expr->as<Product>());
  if (expr->is<Sum>()) {
    auto s = std::make_shared<Sum>();
    for (auto&& summand : *expr) {
      s->append(append_kramers(summand, index_replacements));
    }
    return s;
  }
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;
  throw std::runtime_error("append_kramers: unsupported Expr type");
}

namespace {

/// Walk @p expr and collect the unique Kramers-eligible (i.e.,
/// `Kramers::any`) indices in encounter order. Indices that already
/// have a specific Kramers state (up/down) are skipped.
container::svector<Index> collect_kramers_indices(const ExprPtr& expr) {
  container::svector<Index> seen;
  container::set<Index> seen_set;
  auto consider = [&](const Index& idx) {
    auto qns = idx.space().qns().to_int32();
    if ((qns & mask_v<Kramers>) != static_cast<bitset_t>(Kramers::any)) {
      // Either no Kramers state set or already specialized — skip.
      return;
    }
    if (seen_set.insert(idx).second) seen.push_back(idx);
  };
  expr->visit(
      [&](const ExprPtr& e) {
        if (e->is<Tensor>()) {
          for (const auto& idx : e->as<Tensor>().const_braket()) {
            consider(idx);
          }
        }
      },
      /* recursive = */ true);
  return seen;
}

}  // namespace

ExprPtr kramers_trace(const ExprPtr& expr) {
  // Enumerate the 2^N Kramers configurations of the expression's spinor
  // (Kramers::any) indices, then apply TRS pairing: configs (X, X̄) related
  // by full bar-flip are complex conjugates of each other (per the
  // (-1)^k full-bar identity, with the per-tensor signs cancelling for
  // products in which every barred index appears in an even number of
  // tensors — true for any all-internal-index contraction, MP2 included).
  //
  // For each TRS pair we materialize only the "lex-smaller" configuration
  // X explicitly; the X̄ partner is emitted as `conjugate(T_X)`. Tensors
  // in the conj branch carry a `*` label suffix (per the convention in
  // the conjugate/real_part/imaginary_part scaffolding above) which the
  // evaluator dispatches into a `.conj()` call on the underlying numeric
  // tensor. Term count is unchanged, but the symbolic structure exposes
  // the conjugate-pair relationship for downstream simplification (e.g.,
  // `real_part` collapse in real-scalar callers like SQ_MP2).
  //
  // Standard SeQuant simplifications (dummy renaming, antisym fold, etc.)
  // are NOT applied here — caller invokes `sequant::canonicalize` if it
  // wants those.

  const auto indices = collect_kramers_indices(expr);
  const std::size_t n = indices.size();
  if (n == 0) return expr;
  if (n > 30) {
    throw std::runtime_error(
        "kramers_trace: too many Kramers indices (would emit > 2^30 terms)");
  }

  const std::uint64_t n_configs = 1ull << n;
  const std::uint64_t mask = n_configs - 1;

  // Pass 1: build the "X" term for the lex-smaller half of TRS pairs.
  container::map<std::uint64_t, ExprPtr> x_terms;
  for (std::uint64_t cfg = 0; cfg < n_configs; ++cfg) {
    if (cfg > (mask ^ cfg)) continue;  // partner already (or will be) handled
    container::map<Index, Index> repl;
    for (std::size_t i = 0; i < n; ++i) {
      const bool is_down = (cfg >> i) & 1u;
      Index new_idx = is_down ? make_kramers_dn(indices[i])
                              : make_kramers_up(indices[i]);
      repl[indices[i]] = new_idx;
    }
    x_terms[cfg] = append_kramers(expr, repl);
  }

  // Pass 2: assemble the Sum. The X̄ partner of cfg=k is k^mask; whichever
  // is lex-smaller is the "X" (real) term, the other is `conjugate(X)`.
  auto result = std::make_shared<Sum>();
  for (std::uint64_t cfg = 0; cfg < n_configs; ++cfg) {
    const std::uint64_t partner = mask ^ cfg;
    if (cfg <= partner) {
      result->append(x_terms.at(cfg));
    } else {
      result->append(conjugate(x_terms.at(partner)));
    }
  }
  return result;
}

}  // namespace sequant::mbpt
