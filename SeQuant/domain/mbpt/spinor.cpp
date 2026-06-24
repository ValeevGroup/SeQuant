#include <SeQuant/domain/mbpt/spinor.hpp>

#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <bit>
#include <cstdint>
#include <numeric>
#include <optional>
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

/// Returns the first antisymmetrizer (Â) tensor found anywhere in @p e, or
/// nullopt. In a CC residual term Â is a leading factor whose bra/ket are the
/// external (virtual/occupied) index groups.
std::optional<Tensor> find_antisymmetrizer(const ExprPtr& e) {
  if (e->is<Tensor>()) {
    const auto& t = e->as<Tensor>();
    if (t.label() == reserved::antisymm_label()) return t;
    return std::nullopt;
  } else if (e->is<Product>()) {
    for (const auto& f : e->as<Product>()) {
      if (auto r = find_antisymmetrizer(f)) return r;
    }
    return std::nullopt;
  } else if (e->is<Sum>()) {
    for (const auto& s : *e) {
      if (auto r = find_antisymmetrizer(s)) return r;
    }
    return std::nullopt;
  }
  return std::nullopt;
}

/// Removes every antisymmetrizer (Â) tensor from @p e (factor it out / set it
/// aside). The external indices it carried remain on the other tensors; Â is
/// reattached later by the A-expand stage. Handles Product and Sum.
ExprPtr strip_antisymmetrizer(const ExprPtr& e) {
  if (e->is<Product>()) {
    const auto& p = e->as<Product>();
    auto result = std::make_shared<Product>();
    result->scale(p.scalar());
    for (const auto& f : p) {
      if (f->is<Tensor>() &&
          f->as<Tensor>().label() == reserved::antisymm_label())
        continue;
      result->append(1, f, Product::Flatten::No);
    }
    return result;
  } else if (e->is<Sum>()) {
    auto result = std::make_shared<Sum>();
    for (const auto& s : *e) result->append(strip_antisymmetrizer(s));
    return result;
  }
  return e;
}

/// Applies a bit permutation @p perm (bit p -> bit perm[p]) to config @p m.
std::uint64_t apply_bit_perm(const container::svector<std::size_t>& perm,
                             std::uint64_t m, std::size_t n) {
  std::uint64_t out = 0;
  for (std::size_t p = 0; p < n; ++p)
    if ((m >> p) & 1u) out |= (std::uint64_t{1} << perm[p]);
  return out;
}

}  // namespace

container::svector<container::svector<std::uint64_t>> kramers_config_orbits(
    std::size_t n,
    const container::svector<container::svector<std::size_t>>& bit_perms,
    bool use_T) {
  SEQUANT_ASSERT(n <= 62);
  const std::uint64_t full_mask = (std::uint64_t{1} << n) - 1;

  container::svector<container::svector<std::uint64_t>> orbits;
  container::set<std::uint64_t> seen;
  for (std::uint64_t m = 0; m <= full_mask; ++m) {
    if (seen.find(m) != seen.end()) continue;

    // BFS the orbit of m under the generators.
    container::set<std::uint64_t> orbit{m};
    container::svector<std::uint64_t> frontier{m};
    while (!frontier.empty()) {
      container::svector<std::uint64_t> next;
      for (const std::uint64_t x : frontier) {
        container::svector<std::uint64_t> images;
        for (const auto& perm : bit_perms)
          images.push_back(apply_bit_perm(perm, x, n));
        if (use_T) images.push_back((~x) & full_mask);
        for (const std::uint64_t y : images)
          if (orbit.insert(y).second) next.push_back(y);
      }
      frontier = std::move(next);
    }

    // container::set iterates ascending, so front() is the orbit minimum.
    container::svector<std::uint64_t> members(orbit.begin(), orbit.end());
    for (const std::uint64_t x : members) seen.insert(x);
    orbits.push_back(std::move(members));
  }
  return orbits;
}

container::svector<KramersBlock> kramers_external_blocks(
    std::size_t n,
    const container::svector<container::svector<std::size_t>>& antisym_perms,
    bool use_T,
    const container::svector<container::svector<std::size_t>>& symm_perms) {
  SEQUANT_ASSERT(n <= 62);
  const std::uint64_t full_mask = (std::uint64_t{1} << n) - 1;

  using Perm = container::svector<std::size_t>;
  Perm identity(n);
  std::iota(identity.begin(), identity.end(), std::size_t{0});

  // transform from canonical to a config: block(cfg) = sign*[conj]*perm(canon)
  struct Xform {
    int sign;
    bool conj;
    Perm perm;
  };
  // compose: apply generator g (cfg x -> y) on top of transform tx (canon -> x)
  auto compose = [n](const Xform& tx, int g_sign, bool g_conj,
                     const Perm& g_perm) {
    Xform ty;
    ty.sign = tx.sign * g_sign;
    ty.conj = tx.conj ^ g_conj;
    ty.perm.resize(n);
    for (std::size_t k = 0; k < n; ++k) ty.perm[k] = tx.perm[g_perm[k]];
    return ty;
  };

  container::svector<KramersBlock> blocks;
  container::set<std::uint64_t> seen;
  for (std::uint64_t m = 0; m <= full_mask; ++m) {
    if (seen.find(m) != seen.end()) continue;

    // BFS the orbit, tracking each member's transform from the canonical m.
    container::map<std::uint64_t, Xform> xform;
    xform.emplace(m, Xform{+1, false, identity});
    container::svector<std::uint64_t> frontier{m};
    while (!frontier.empty()) {
      container::svector<std::uint64_t> next;
      for (const std::uint64_t x : frontier) {
        const Xform tx =
            xform.at(x);  // copy: emplace below may invalidate refs
        // antisymmetric generators (residual external antisymmetry): sign -1,
        // no conj, permutation = perm
        for (const auto& perm : antisym_perms) {
          const std::uint64_t y = apply_bit_perm(perm, x, n);
          if (xform.find(y) == xform.end()) {
            xform.emplace(y, compose(tx, -1, false, perm));
            next.push_back(y);
          }
        }
        // symmetric generators (e.g. particle interchange sigma on a raw,
        // non-antisymmetrized integral): sign +1, no conj, permutation = perm
        for (const auto& perm : symm_perms) {
          const std::uint64_t y = apply_bit_perm(perm, x, n);
          if (xform.find(y) == xform.end()) {
            xform.emplace(y, compose(tx, +1, false, perm));
            next.push_back(y);
          }
        }
        // global time reversal: complement, conj, sign (-1)^(#down of x)
        if (use_T) {
          const std::uint64_t y = (~x) & full_mask;
          if (xform.find(y) == xform.end()) {
            const int t_sign = (std::popcount(x) & 1) ? -1 : +1;
            xform.emplace(y, compose(tx, t_sign, true, identity));
            next.push_back(y);
          }
        }
      }
      frontier = std::move(next);
    }

    KramersBlock block;
    block.canonical = m;  // orbit-min (m increases, first unseen is the min)
    for (const auto& [cfg, t] : xform) {  // map iterates ascending by cfg
      block.members.push_back({cfg, t.sign, t.conj, t.perm});
      seen.insert(cfg);
    }
    blocks.push_back(std::move(block));
  }
  return blocks;
}

bool has_antisymmetrizer(const ExprPtr& expr) {
  return find_antisymmetrizer(expr).has_value();
}

container::svector<ExprPtr> closed_shell_kramers_CC_trace(const ExprPtr& expr) {
  // Stage 1: factor out the antisymmetrizer Â (kept, not expanded). Its bra/ket
  // are the external virtual/occupied index groups.
  auto A = find_antisymmetrizer(expr);
  SEQUANT_ASSERT(A.has_value() &&
                 "closed_shell_kramers_CC_trace: residual must carry Â");

  // External indices in bit order: bra group (virtuals) then ket group
  // (occupieds). Bit k = Kramers label of ext[k] (0 = up, 1 = down).
  container::svector<Index> ext;
  for (const auto& idx : A->bra()) ext.push_back(idx);
  for (const auto& idx : A->ket()) ext.push_back(idx);
  const std::size_t n_bra = A->bra_rank();
  const std::size_t n = ext.size();

  // Stage 2 (external fold): generators = the S_k adjacent transpositions of
  // each external group (rank-general: the residual's antisymmetry in each
  // external pair/tuple) + global time reversal T. Doubles -> 5 blocks.
  container::svector<container::svector<std::size_t>> gens;
  auto add_group_swaps = [&](std::size_t lo, std::size_t hi) {
    for (std::size_t p = lo; p + 1 < hi; ++p) {
      container::svector<std::size_t> perm(n);
      std::iota(perm.begin(), perm.end(), std::size_t{0});
      std::swap(perm[p], perm[p + 1]);
      gens.push_back(std::move(perm));
    }
  };
  add_group_swaps(0, n_bra);  // virtual external group
  add_group_swaps(n_bra, n);  // occupied external group

  const auto orbits = kramers_config_orbits(n, gens, /*use_T=*/true);

  // Factor Â out (kept symbolic, reattached at the A-expand stage). The
  // external indices remain on the other tensors.
  const ExprPtr inner = strip_antisymmetrizer(expr);

  // Emit one block per external representative.
  container::svector<ExprPtr> blocks;
  blocks.reserve(orbits.size());
  for (const auto& orbit : orbits) {
    const std::uint64_t cfg = orbit.front();  // canonical (orbit-min)

    // Stage 2a: label the external indices for this block.
    container::map<Index, Index> ext_repl;
    for (std::size_t k = 0; k < n; ++k) {
      const bool down = (cfg >> k) & 1u;
      ext_repl.emplace(ext[k],
                       down ? make_spinbeta(ext[k]) : make_spinalpha(ext[k]));
    }
    ExprPtr block = append_spin(inner, ext_repl);

    // Stage 2b (internal fold): the still-unlabeled (spin-any) indices of a
    // term are its contracted/internal ones. Enumerate all their Kramers
    // configurations and sum (no Ms filter — relativistic); canonicalize then
    // performs the particle-interchange (sigma) merge. The internal-T-reach
    // (extra fold in the self-complementary external block) is deferred.
    //
    // The fold MUST be applied per term (summand): different terms of one
    // residual carry different internal-index sets, and a term invariant under
    // a given internal index (e.g. the driver, with no internal index) would be
    // over-counted 2^k if the whole block were folded over the union of all
    // terms' internal indices.
    auto fold_term = [&](const ExprPtr& term) -> ExprPtr {
      const auto used =
          get_used_indices<container::set<Index, Index::LabelCompare>>(term);
      container::svector<Index> internal;
      for (const auto& idx : used)
        if (to_spin(idx.space().qns()) == Spin::any) internal.push_back(idx);
      if (internal.empty()) return term;
      const std::size_t ni = internal.size();
      SEQUANT_ASSERT(ni <= 62);
      const std::uint64_t nconfigs = pow2(ni);
      auto sum = std::make_shared<Sum>();
      for (std::uint64_t ic = 0; ic < nconfigs; ++ic) {
        container::map<Index, Index> int_repl;
        for (std::size_t k = 0; k < ni; ++k) {
          const bool down = (ic >> k) & 1u;
          int_repl.emplace(internal[k], down ? make_spinbeta(internal[k])
                                             : make_spinalpha(internal[k]));
        }
        sum->append(append_spin(term, int_repl));
      }
      return sum;
    };

    if (block->is<Sum>()) {
      auto folded = std::make_shared<Sum>();
      for (const auto& summand : *block) folded->append(fold_term(summand));
      block = folded;
    } else {
      block = fold_term(block);
    }

    canonicalize(block);  // sigma merge + dummy canonicalization
    rapid_simplify(block);
    blocks.push_back(block);
  }
  return blocks;
}

ExprPtr closed_shell_kramers_trace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups,
    bool fold_T, bool expand_g) {
  if (expr->is<Constant>() || expr->is<Variable>()) return expr;

  // Step 0/6: optionally expand the integral `g`'s antisymmetry (Kramers-free,
  // so no Ms filtering), keeping the amplitude `t` antisymmetric (the level-1
  // form). When expand_g is false, `g` stays antisymmetric (ḡ) so the evaluator
  // fetches the factory [as] block directly — the cross-Kramers antisymmetry is
  // then handled inside the integral, not by a raw config sum. (Expanding g is
  // a later optimization stage.)
  ExprPtr traced_input =
      expand_g ? expand_label_antisymm(expr, L"g") : expr->clone();
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
    if (fold_T) {
      const std::uint64_t comp = cfg ^ full_mask;
      if (cfg > comp) continue;  // keep canonical T-rep (cfg < comp for n >= 1)
    }

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

  // Emit one representative per orbit. With the T-fold, coefficient = 2 x
  // multiplicity and each block is RealPart-wrapped (A + A* = 2 Re A). Without
  // it, every configuration is emitted verbatim (coefficient = multiplicity, no
  // RealPart) — the complex sum whose real part the caller takes.
  auto result = std::make_shared<Sum>();
  for (const auto& [block, mult] : reps) {
    if (fold_T)
      result->append(ex<Constant>(2 * mult) * real_part(block));
    else
      result->append(ex<Constant>(mult) * block);
  }
  return result;
}

}  // namespace sequant::mbpt
