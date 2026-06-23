#include <SeQuant/domain/mbpt/spinor.hpp>

#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <cstdint>
#include <string_view>
#include <utility>

namespace sequant::mbpt {

namespace {

/// Expands the bra antisymmetry of every tensor named @p label into the
/// NonSymm form, leaving all other tensors untouched. Intended to be run on a
/// Kramers-FREE expression (all indices spin-`any`): `expand_antisymm`'s
/// `ms_conserving_columns` guard then keeps every permutation (no Ms filter),
/// which is what we want — Kramers is not an Ms-conserved label.
ExprPtr expand_label_antisymm(const ExprPtr& expr, std::wstring_view label) {
  if (expr->is<Tensor>()) {
    const auto& t = expr->as<Tensor>();
    return t.label() == label ? expand_antisymm(t) : expr;
  } else if (expr->is<Product>()) {
    const auto& p = expr->as<Product>();
    auto result = std::make_shared<Product>();
    result->scale(p.scalar());
    for (const auto& factor : p) {
      if (factor->is<Tensor>() && factor->as<Tensor>().label() == label) {
        result->append(1, expand_antisymm(factor->as<Tensor>()),
                       Product::Flatten::No);
      } else {
        result->append(1, factor, Product::Flatten::No);
      }
    }
    ExprPtr r = result;
    expand(r);  // distribute the antisymmetrizer-expansion Sum factor
    rapid_simplify(r);
    return r;
  } else if (expr->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (const auto& summand : *expr) {
      result->append(expand_label_antisymm(summand, label));
    }
    return result;
  }
  return expr;
}

}  // namespace

ExprPtr closed_shell_kramers_trace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups) {
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  // Step 0/6: expand the integral `g`'s antisymmetry (Kramers-free, so no Ms
  // filtering), keep the amplitude `t` antisymmetric (the level-1 form).
  ExprPtr traced_input = expand_label_antisymm(expr, L"g");
  canonicalize(traced_input);
  rapid_simplify(traced_input);

  // Steps 1-2: classify indices and build groups (each internal index its own
  // group; external groups appended verbatim). No Ms / rank assumptions.
  const auto all_indices =
      get_used_indices<container::set<Index, Index::LabelCompare>>(
          traced_input);

  container::set<Index> ext_indices;
  for (const auto& group : ext_index_groups) {
    for (const auto& idx : group) {
      Index x = idx;
      x.reset_tag();
      ext_indices.insert(std::move(x));
    }
  }

  using IndexGroup = container::svector<Index>;
  container::svector<IndexGroup> index_groups;
  for (const auto& idx : all_indices) {
    if (ext_indices.find(idx) == ext_indices.end())
      index_groups.emplace_back(IndexGroup(1, idx));
  }
  for (const auto& group : ext_index_groups) {
    index_groups.emplace_back(IndexGroup(group.begin(), group.end()));
  }

  // Step 3: enumerate the 2^n Kramers configurations and fold conjugate
  // (global time-reversal) partners. The config integer assigns bit k to
  // index_groups[k]; bit 0 = Kramers-up (Spin::alpha), bit 1 = down
  // (Spin::beta). The global-T partner is the full one's-complement; for a
  // closed contraction block(comp) = conj(block(cfg)), so the pair sums to
  // 2 Re(block). We keep the lexicographically smaller of each T-pair.
  const std::size_t n = index_groups.size();
  SEQUANT_ASSERT(n <= 62);
  const std::uint64_t nconfigs = pow2(n);
  const std::uint64_t full_mask = nconfigs - 1;

  // Accumulate canonical representatives, merging configurations that
  // `canonicalize` made identical (the particle-interchange / sigma fold) and
  // summing their T-fold multiplicities. RealPart is opaque to
  // `simplify`/`canonicalize`, so we merge explicitly here (exact Expr equality
  // of the canonicalized blocks) rather than relying on the Sum machinery.
  container::svector<std::pair<ExprPtr, std::int64_t>> reps;
  for (std::uint64_t cfg = 0; cfg < nconfigs; ++cfg) {
    const std::uint64_t comp = cfg ^ full_mask;
    if (cfg > comp) continue;  // keep canonical T-rep (cfg < comp for n >= 1)

    container::map<Index, Index> replacements;
    for (std::size_t k = 0; k < n; ++k) {
      const bool down = (cfg >> k) & 1u;
      for (const auto& idx : index_groups[k]) {
        replacements.emplace(idx,
                             down ? make_spinbeta(idx) : make_spinalpha(idx));
      }
    }

    ExprPtr block = append_spin(traced_input, replacements);
    canonicalize(block);  // particle-interchange (sigma) merge within the block
    rapid_simplify(block);

    bool merged = false;
    for (auto& [rep, mult] : reps) {
      if (*rep == *block) {
        ++mult;
        merged = true;
        break;
      }
    }
    if (!merged) reps.emplace_back(block, std::int64_t{1});
  }

  // Emit one representative per orbit; coefficient = 2 (T-fold) x multiplicity.
  auto result = std::make_shared<Sum>();
  for (const auto& [block, mult] : reps) {
    result->append(ex<Constant>(2 * mult) * real_part(block));
  }
  return result;
}

}  // namespace sequant::mbpt
