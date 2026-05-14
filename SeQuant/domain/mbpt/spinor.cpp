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
  // MVP (Phase 2b/i): enumerate the 2^N Kramers configurations of the
  // expression's spinor (Kramers::any) indices and emit a Sum of substituted
  // copies. Each summand is a Product (or Sum/Tensor) where every Kramers
  // index has been pinned to up/down.
  //
  // Per-tensor TRS canonicalization (the (-1)^k bit-table fold) and
  // cross-tensor simplification land in subsequent commits.
  //
  // The eventual evaluator (LCAOFactory) sees each summand as a concrete
  // (g_block, t_block, ...) contraction; LCAOFactory's own canonical-storage
  // layer (Phase 1.5) ensures we still pay only the orbit-rep integral cost.

  const auto indices = collect_kramers_indices(expr);
  const std::size_t n = indices.size();
  if (n == 0) return expr;
  if (n > 30) {
    throw std::runtime_error(
        "kramers_trace: too many Kramers indices (would emit > 2^30 terms)");
  }

  auto result = std::make_shared<Sum>();
  const std::uint64_t n_configs = 1ull << n;
  for (std::uint64_t cfg = 0; cfg < n_configs; ++cfg) {
    container::map<Index, Index> repl;
    for (std::size_t i = 0; i < n; ++i) {
      const bool is_down = (cfg >> i) & 1u;
      Index new_idx = is_down ? make_kramers_dn(indices[i])
                              : make_kramers_up(indices[i]);
      repl[indices[i]] = new_idx;
    }
    result->append(append_kramers(expr, repl));
  }
  return result;
}

}  // namespace sequant::mbpt
