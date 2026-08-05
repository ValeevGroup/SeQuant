#include <SeQuant/domain/mbpt/bernoulli.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/range/primitives.hpp>

#include <algorithm>
#include <string>
#include <utility>

// Bernoulli expansion of the unitary-CC similarity-transformed Hamiltonian
// H̄ = e^{−σ} H e^{σ}, σ = T − T† (anti-Hermitian). For UCC the plain BCH
// series does not terminate, because σ mixes excitation and de-excitation.
// This file implements the Bernoulli-number resummation of 10.1063/1.5030344
// that fixes that. H splits as F (Fock) + V (fluctuation potential); every
// operator O splits into O_N (its pure excitation/de-excitation part) and
// O_R = O − O_N.
//
// Equation numbers below are all from 10.1063/1.5030344, Sec. III B:
// superoperator inversion Eqs. (36)-(39); Bernoulli numbers B₁=−1/2, B₂=1/12,
// B₃=0, B₄=−1/720, Eq. (40); the N/R split and the UCC amplitude condition
// V̄_N = 0, above and at Eq. (43); the H̄ recursion, Eq. (44); the sum
// H̄ = Σ_k H̄^k, Eq. (45); H̄⁰..H̄⁴, Eqs. (46)-(50).
//
// The F-cancellation: at a canonical (converged) HF reference, F has no
// occupied-virtual block (Eq. (32)), so F enters H̄ only through H̄¹ (stated
// just below Eq. (50)). Every H̄^k built below relies on this.

namespace {

/// Returns the single residual fermionic NormalOperator carried by @p term, or
/// nullptr when it has none (a pure scalar / fully-contracted term). Every term
/// produced by wick_reduce is either a bare NormalOperator or a Product with at
/// most one NormalOperator factor times tensor coefficients.
const sequant::NormalOperator<sequant::Statistics::FermiDirac>* find_nop(
    const sequant::ExprPtr& term) {
  using namespace sequant;
  if (term.is<NormalOperator<Statistics::FermiDirac>>())
    return &term.as<NormalOperator<Statistics::FermiDirac>>();

  if (term.is<Product>()) {
    const NormalOperator<Statistics::FermiDirac>* found = nullptr;
    for (const auto& f : term.as<Product>().factors())
      if (f.is<NormalOperator<Statistics::FermiDirac>>()) {
        // The one-residual-operator invariant is load-bearing: N/R
        // classification reads this operator alone, so a second one would be
        // silently ignored and misclassify the term.
        SEQUANT_ASSERT(!found &&
                       "find_nop: term carries >1 NormalOperator; wick_reduce "
                       "is expected to leave at most one residual operator");
        found = &f.as<NormalOperator<Statistics::FermiDirac>>();
      }

    return found;
  }
  return nullptr;
}

/// Classifies one block-resolved term as N or R, per the O_N/O_R split above
/// Eq. (43) of 10.1063/1.5030344. A term is N iff its single residual
/// NormalOperator is a pure excitation (all creators pure-unoccupied AND all
/// annihilators pure-occupied) or a pure de-excitation (the reverse), with
/// rank ≤ @p cutoff. A term with no residual NormalOperator is rank-preserving,
/// hence R.
///
/// Rank > @p cutoff falls to R rather than being dropped. The R filter drops a
/// term only because the amplitude condition V̄_N = 0 (Eq. (43)) makes it
/// zero, and that condition holds for rank ≤ N only, since σ here is
/// truncated at rank N. Eq. (43) itself states O_N with no rank limit,
/// because there σ carries every rank.
bool is_N_term(const sequant::ExprPtr& term, std::size_t cutoff) {
  using namespace sequant;
  auto isr = get_default_context().index_space_registry();

  const auto* nop = find_nop(term);
  if (!nop) return false;  // no residual operator => rank-preserving => R

  const auto ncre = ranges::distance(nop->creators());
  const auto nann = ranges::distance(nop->annihilators());
  if (static_cast<std::size_t>(std::max(ncre, nann)) > cutoff) return false;

  auto all_unocc = [&](auto&& ops) {
    return ranges::all_of(ops, [&](const auto& o) {
      return isr->is_pure_unoccupied(o.index().space());
    });
  };
  auto all_occ = [&](auto&& ops) {
    return ranges::all_of(ops, [&](const auto& o) {
      return isr->is_pure_occupied(o.index().space());
    });
  };

  const bool pure_exc =
      all_unocc(nop->creators()) && all_occ(nop->annihilators());
  const bool pure_deexc =
      all_occ(nop->creators()) && all_unocc(nop->annihilators());

  return pure_exc || pure_deexc;
}

}  // namespace

namespace sequant::mbpt::bernoulli {

namespace detail {

ExprPtr wick_reduce(ExprPtr expr) {
  simplify(expr);
  FWickTheorem wick{expr};
  // use_topology defaults to ON, so it must be turned off explicitly. It keeps
  // one representative per symmetry-equivalent contraction class and multiplies
  // by the class size, weight bookkeeping that holds only on the
  // fully-contracted path. On this partial-contraction path it silently
  // rescales the terms carrying a symmetric amplitude pair, and the damage
  // shows up only under projection.
  wick.use_topology(false).full_contractions(false);
  auto result = wick.compute(/*count_only=*/false,
                             /*skip_input_canonicalization=*/true);
  simplify(result);
  return result;
}

ExprPtr wick_commutator(const ExprPtr& A, const ExprPtr& B) {
  // A and B are built independently, so their summed indices are local to each.
  // If both use the same labels (a block-resolved R/N part and sigma both carry
  // a/i), the product A*B fuses two independent summations. That corrupts the
  // contraction. Reindex B to fresh temporaries. Canonicalization restores tidy
  // labels.
  container::map<Index, Index> repl;
  for (const auto& idx : get_used_indices(B))
    repl.emplace(idx, Index::make_tmp_index(idx.space()));
  const auto Bd = repl.empty() ? B : transform_expr(B, repl);
  return wick_reduce(simplify(A * Bd - Bd * A));
}

namespace {

/// Core of expand_to_blocks for input already in wick_reduce'd form. Skipping
/// the reduction is an identity: wick_reduce is idempotent (terms with a single
/// residual NormalOperator admit no further contractions). @p expr is not
/// mutated.
ExprPtr expand_to_blocks_reduced(const ExprPtr& expr) {
  auto isr = get_default_context().index_space_registry();
  const auto& bases = isr->base_spaces();

  auto is_base_space = [&](const IndexSpace& sp) {
    return ranges::any_of(bases, [&](const auto& b) { return b == sp; });
  };

  // Split each general index over the hole and particle base spaces only. A
  // general index also spans the registry's other base spaces (under the SR
  // convention the frozen-core "o" and inactive-virtual "g"), but terms landing
  // in those are annihilated by the single-reference projection onto the
  // hole/particle manifolds, so dropping them changes no projected quantity.
  // This keeps the expansion 2-way per index instead of 4-way, which otherwise
  // compounds across the nested commutators. Falls back to all base spaces if
  // the registry defines no hole/particle split.
  const auto& hole_t = isr->hole_space(/*nulltype_ok=*/true);
  const auto& particle_t = isr->particle_space(/*nulltype_ok=*/true);
  auto physical = [&](const IndexSpace& b) {
    return hole_t.includes(b.type()) || particle_t.includes(b.type());
  };

  auto expand_term = [&](const ExprPtr& term) -> ExprPtr {
    // collect the residual NormalOperator's distinct general (non-base) indices
    const auto* nop = find_nop(term);
    if (!nop)
      return term->clone();  // pure scalar/contraction: nothing to split
    container::svector<Index> gens;
    for (const auto& op : nop->creann()) {
      if (!is_base_space(op.index().space()) &&
          ranges::none_of(gens, [&](const auto& g) { return g == op.index(); }))
        gens.push_back(op.index());
    }
    if (gens.empty()) return term->clone();
    // candidate base spaces per general index: base b is a sub-block of the
    // general space iff its type bits are included and its quantum numbers
    // match (stay within the same spin sector).
    container::svector<container::svector<IndexSpace>> choices;
    for (const auto& g : gens) {
      container::svector<IndexSpace> c, c_all;
      for (const auto& b : bases)
        if (b.qns() == g.space().qns() && g.space().type().includes(b.type())) {
          c_all.push_back(b);
          if (physical(b)) c.push_back(b);
        }
      choices.push_back(c.empty() ? c_all : c);
    }
    // cartesian product of assignments => sum of transformed terms;
    // accumulate via Sum::append (linear) rather than operator+, which
    // deep-copies the accumulated Sum on every call (quadratic)
    auto sum = std::make_shared<Sum>();
    container::svector<std::size_t> idx(gens.size(), 0);
    for (;;) {
      container::map<Index, Index> repl;
      for (std::size_t k = 0; k < gens.size(); ++k)
        // fresh ordinal (not gens[k].ordinal()): reusing the general index's
        // ordinal would collide with any pre-existing definite index of the
        // same base space and ordinal already in the term (e.g. an a/i index
        // from an amplitude in a commutator result), producing a duplicate
        // index. A globally-unique temporary is disjoint by construction;
        // canonicalization restores tidy labels.
        repl.emplace(gens[k], Index::make_tmp_index(choices[k][idx[k]]));
      sum->append(transform_expr(term, repl));
      // increment mixed-radix counter over the assignments
      std::size_t k = 0;
      for (; k < gens.size(); ++k) {
        if (++idx[k] < choices[k].size()) break;
        idx[k] = 0;
      }
      if (k == gens.size()) break;
    }
    return ExprPtr{sum};
  };

  // transform_sum_expr maps in parallel, canonicalizes each result, and
  // accumulates into a HashingAccumulator
  ExprPtr out;
  if (expr.is<Sum>()) {
    out = transform_sum_expr(expr.as<Sum>().summands(), expand_term);
  } else {
    out = expand_term(expr);
  }
  simplify(out);
  return out;
}

/// Keeps only the N terms of block-resolved @p bx (shared tail of N_part and
/// N_part_reduced).
ExprPtr keep_N_terms(const ExprPtr& bx, std::size_t cutoff) {
  if (bx.is<Sum>()) {
    auto out = std::make_shared<Sum>();
    for (const auto& t : bx.as<Sum>())
      if (is_N_term(t, cutoff)) out->append(t);
    return out->empty() ? ex<Constant>(0) : simplify(ExprPtr{out});
  }
  return is_N_term(bx, cutoff) ? bx : ex<Constant>(0);
}

/// N_part for input already in wick_reduce'd form.
ExprPtr N_part_reduced(const ExprPtr& reduced, std::size_t cutoff) {
  return keep_N_terms(expand_to_blocks_reduced(reduced), cutoff);
}

/// R_part for input already in wick_reduce'd form.
ExprPtr R_part_reduced(const ExprPtr& reduced, std::size_t cutoff) {
  return simplify(reduced - N_part_reduced(reduced, cutoff));
}

}  // namespace

/// Identity expansion of every general index into its base sub-blocks (see
/// header): after expansion every residual index is definite, so the N/R
/// classifier can act on it.
ExprPtr expand_to_blocks(const ExprPtr& expr_in) {
  return expand_to_blocks_reduced(wick_reduce(expr_in->clone()));
}

/// N part of @p expr at truncation @p cutoff (see header): block-resolve, then
/// keep only the pure excitation / de-excitation terms.
ExprPtr N_part(const ExprPtr& expr, std::size_t cutoff) {
  return keep_N_terms(expand_to_blocks(expr), cutoff);
}

/// R part of @p expr at truncation @p cutoff (see header): the reduced operator
/// minus its N part. Because expand_to_blocks is an identity
/// (N ⊎ R = expr as operators), R = expr − N holds exactly while expr stays in
/// its compact (general-index) form. Only N is block-resolved. The result
/// equals the fully block-resolved remainder, and the compact expr makes the
/// nested commutators that consume R operate on far fewer terms.
ExprPtr R_part(const ExprPtr& expr, std::size_t cutoff) {
  auto reduced = wick_reduce(expr->clone());
  return R_part_reduced(reduced, cutoff);
}

}  // namespace detail

/// Assembles H̄ order by order (see header), summing H̄⁰..H̄^rank of Eq. (45).
/// Each H̄^k below is a direct transcription of its equation. A subscript R/N on
/// a commutator means "form the commutator, then keep only its R/N part before
/// the next nesting".
ExprPtr hbar(std::size_t N, std::size_t rank, bool skip1) {
  if (rank > 4)
    throw Exception("bernoulli::hbar: only ranks [0,4] are implemented");

  using namespace detail;
  const auto cutoff = N;
  const auto F = op::tensor::F();
  const auto V = op::tensor::h(2);
  const auto T = op::tensor::T(N, skip1);
  const auto sigma = simplify(T - adjoint(T));  // σ = T − T†

  // Every term of H̄^k is a nested commutator [[..[V_{p0},σ]_{f0}..],σ]_{f_k}
  // with a per-level N/R/A partition tag applied after each commutator ('A' =
  // no filter). nest(p0, f) evaluates such a node, memoizing every prefix
  // (key = p0 + tags applied so far): the terms share prefixes both within a
  // rank (the 9 rank-4 terms have only 3 distinct level-1 and 6 level-2 nodes)
  // and across ranks (all four ranks share the same three level-1 nodes), so
  // the memo avoids recomputing them. Reusing a memoized ExprPtr is safe:
  // expression composition deep-copies operands (Product/Sum append clone), so
  // wick_commutator does not mutate its arguments. Commutator outputs are
  // already wick_reduce'd, so the reduced-input N/R filters apply.
  container::map<std::string, ExprPtr> memo;
  auto nest = [&](char p0, const char* f) -> ExprPtr {
    SEQUANT_ASSERT((p0 == 'A' || p0 == 'N' || p0 == 'R') &&
                   "bernoulli::hbar: partition tag must be one of A, N, R");
    // grow `key` in place rather than deriving it from the memo iterator:
    // container::map is a flat_map, whose insertions invalidate iterators
    std::string key{p0};
    auto it = memo.find(key);
    if (it == memo.end()) {
      ExprPtr base = (p0 == 'N')   ? N_part(V, cutoff)
                     : (p0 == 'R') ? R_part(V, cutoff)
                                   : V;
      it = memo.emplace(key, std::move(base)).first;
    }
    ExprPtr op = it->second;
    for (int i = 0; f[i] != '\0'; ++i) {
      SEQUANT_ASSERT((f[i] == 'A' || f[i] == 'N' || f[i] == 'R') &&
                     "bernoulli::hbar: partition tag must be one of A, N, R");
      key += f[i];
      it = memo.find(key);
      if (it == memo.end()) {
        auto cx = wick_commutator(op, sigma);
        ExprPtr filtered = (f[i] == 'R')   ? R_part_reduced(cx, cutoff)
                           : (f[i] == 'N') ? N_part_reduced(cx, cutoff)
                                           : cx;
        it = memo.emplace(key, std::move(filtered)).first;
      }
      op = it->second;
    }
    return op;
  };

  HashingAccumulator acc;
  auto add = [&acc](rational num, const ExprPtr& e) {
    if (e.is<Sum>()) {
      for (const auto& term : e.as<Sum>()) {
        auto scaled = ex<Product>(ExprPtrList{term});
        scaled.as<Product>().scale(num);
        acc.append(std::move(scaled), /*flatten=*/false);
      }
    } else {
      auto scaled = ex<Product>(ExprPtrList{e});
      scaled.as<Product>().scale(num);
      acc.append(std::move(scaled), /*flatten=*/false);
    }
  };

  add(1, simplify(F + V));  // H̄⁰ = F + V   [Eq. (46)]
  if (rank >= 1) {
    // H̄¹ = [F,σ] + ½[V,σ] + ½[V_R,σ]   [Eq. (47)]. F enters H̄ ONLY here
    // (the F-cancellation, stated just below Eq. (50)).
    add(1, wick_commutator(F, sigma));
    add({1, 2}, nest('A', "A"));
    add({1, 2}, nest('R', "A"));
  }
  if (rank >= 2) {
    // H̄² = 1/12[[V_N,σ],σ] + ¼[[V,σ]_R,σ] + ¼[[V_R,σ]_R,σ]   [Eq. (48)]
    add({1, 12}, nest('N', "AA"));
    add({1, 4}, nest('A', "RA"));
    add({1, 4}, nest('R', "RA"));
  }
  if (rank >= 3) {
    // H̄³ = 1/24[[[V_N,σ],σ]_R,σ] + ⅛[[[V,σ]_R,σ]_R,σ] + ⅛[[[V_R,σ]_R,σ]_R,σ]
    //       − 1/24[[[V,σ]_R,σ],σ] − 1/24[[[V_R,σ]_R,σ],σ]   [Eq. (49)]
    add({1, 24}, nest('N', "ARA"));
    add({1, 8}, nest('A', "RRA"));
    add({1, 8}, nest('R', "RRA"));
    add({-1, 24}, nest('A', "RAA"));
    add({-1, 24}, nest('R', "RAA"));
  }
  if (rank >= 4) {
    // H̄⁴ = Eq. (50), the nine order-4 terms produced by the recursion Eq. (44),
    // V̄^{k+1} = σ̂F + X̂⁻¹(σ̂)e^{σ̂}V − Σ_{n≠0} B_n σ̂^n V̄_R^{k}. F is absent here
    // (the F-cancellation). Listed in the paper's order; the outermost tag is
    // always A.
    add({1, 16}, nest('R', "RRRA"));
    add({1, 16}, nest('A', "RRRA"));
    add({1, 48}, nest('N', "ARRA"));
    add({-1, 48}, nest('A', "RARA"));
    add({-1, 48}, nest('R', "RARA"));
    add({-1, 144}, nest('N', "ARAA"));
    add({-1, 48}, nest('A', "RRAA"));
    add({-1, 48}, nest('R', "RRAA"));
    add({-1, 720}, nest('N', "AAAA"));
  }
  auto result = acc.make_expr();
  return simplify(result);
}

}  // namespace sequant::mbpt::bernoulli
