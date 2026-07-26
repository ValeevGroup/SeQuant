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

#include <string>
#include <utility>

// Bernoulli expansion of the unitary-CC similarity-transformed Hamiltonian
// H̄ = e^{−σ} H e^{σ}, σ = T − T† (anti-Hermitian). Because σ mixes excitation
// and de-excitation the plain BCH series does not terminate; the Bernoulli
// expansion rewrites it so that Bernoulli numbers are the expansion
// coefficients, leaving the final truncation at a chosen commutator rank as the
// only approximation. H is split as F (Fock, rank-preserving) + V (fluctuation
// potential), and every operator O is split into O_N (all excitation and
// de-excitation operators) and O_R = O − O_N. At a converged RHF/UHF reference
// two cancellations hold: F survives only in H̄¹, and the higher orders carry
// only R-subscripted inner commutators.
//
// All equation references are to 10.1063/1.5030344 (Sec. III B): superoperator
// inversion Eqs. (36)-(39); Bernoulli numbers B₁=−1/2, B₂=1/12, B₃=0, B₄=−1/720
// Eq. (40); the N/R split and the UCC amplitude condition V̄_N = 0 above and at
// Eq. (43); the iterative recursion for V̄ Eq. (44); the rank-by-rank operators
// H̄⁰..H̄⁴ Eqs. (45)-(50). Cancellation #1 (F enters only H̄¹) is stated just
// below Eq. (50).

namespace {

/// Returns the single residual fermionic NormalOperator carried by @p term, or
/// nullptr when it has none (a pure scalar / fully-contracted term). Every term
/// produced by wick_reduce is either a bare NormalOperator or a Product with
/// exactly one NormalOperator factor times tensor coefficients.
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
                       "is expected to leave exactly one residual operator");
        found = &f.as<NormalOperator<Statistics::FermiDirac>>();
      }
    return found;
  }
  return nullptr;
}

/// Classifies one block-resolved term as N or R (Cancellation #2). A term is N
/// iff its single residual NormalOperator is a pure excitation (all creators
/// pure-unoccupied AND all annihilators pure-occupied) or a pure de-excitation
/// (the mirror), with rank ≤ @p cutoff. A term with no residual NormalOperator
/// is rank-preserving, hence R.
///
/// Rank > @p cutoff falls to R rather than being dropped. The paper's O_N
/// (above Eq. (43): "containing all the excitation operators and de-excitation
/// operators in O") carries no rank cutoff, but this mirrors pdaggerq (`nt_bra
/// > bernoulli_excitation_level -> R`), which is the convention that defines
/// qUCCSD and the one our numbers are validated against. The choice is a real
/// degree of freedom: at cutoff = 2 the HF rank-3 energy is +2.507 mEh, versus
/// +2.428 mEh with no cutoff.
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

/// Operator-valued Wick reduction (see header): reduces a product of
/// normal-ordered operators to a sum of normal-ordered operators, retaining
/// partial contractions so the result is an operator, not a scalar VEV.
ExprPtr wick_reduce(ExprPtr expr) {
  simplify(expr);
  // full_contractions(false) is the whole point: it yields the normal-ordered
  // operator form rather than the scalar VEV. Otherwise mirrors
  // mbpt::tensor::expectation_value_impl. See core/wick.hpp.
  FWickTheorem wick{expr};
  // use_topology MUST be disabled explicitly -- it defaults to ON
  // (wick.hpp: `bool use_topology_ = true`), so merely not asking for it is not
  // enough. It counts one representative per symmetry-equivalent contraction
  // class times a multiplicity; the weight bookkeeping is exercised by the
  // fully-contracted (vacuum-average) path, not by this partial-contraction
  // one. With it on, <0|H̄³|0> stays correct but 16 of the 332 coefficients of
  // <μ|H̄³|0> are rescaled by 2, 1/2, 3, 8/3 or 2/3 -- every one of them a term
  // with a symmetric amplitude pair. With it off, all 332 signatures match
  // pdaggerq exactly. Cost of giving up the optimisation: the rank-3-amplitude
  // derivation goes 6.1 s -> 7.6 s wall.
  wick.use_topology(false).full_contractions(false);
  auto result = wick.compute(/*count_only=*/false,
                             /*skip_input_canonicalization=*/true);
  simplify(result);
  return result;
}

/// Normal-ordered commutator [A, B] = wick_reduce(A·B − B·A) (see header).
ExprPtr wick_commutator(const ExprPtr& A, const ExprPtr& B) {
  // Disjoin B's (bound) indices from A's before forming the product: A and B
  // are independently constructed operator expressions whose summed indices are
  // local to each. If they happen to share labels (e.g. a block-resolved R/N
  // part, which carries definite a/i/o/g indices, commuted with sigma, which
  // also uses a/i), the naive product A*B would identify two independent
  // summations, corrupting the contraction. Reindexing B to globally-fresh
  // temporaries makes the two index sets disjoint; canonicalization restores
  // tidy labels afterward.
  container::map<Index, Index> repl;
  for (const auto& idx : get_used_indices(B))
    repl.emplace(idx, Index::make_tmp_index(idx.space()));
  const auto Bd = repl.empty() ? B : transform_expr(B, repl);
  return wick_reduce(simplify(A * Bd - Bd * A));
}

namespace {

/// Core of expand_to_blocks for input already in wick_reduce'd form (a
/// simplified sum of coefficient × single-NormalOperator terms). Skipping the
/// reduction is an identity: wick_reduce is idempotent (terms with a single
/// residual NormalOperator admit no further contractions). @p expr is not
/// mutated.
ExprPtr expand_to_blocks_reduced(const ExprPtr& expr) {
  auto isr = get_default_context().index_space_registry();
  const auto& bases = isr->base_spaces();

  auto is_base_space = [&](const IndexSpace& sp) {
    return ranges::any_of(bases, [&](const auto& b) { return b == sp; });
  };

  // Only the hole and particle base spaces are physically populated in the
  // single-reference qUCCSD setting; the other base spaces the complete space
  // nominally spans (e.g. the SR "o"/"g" blocks) are empty, so splitting a
  // general index into them only multiplies the term count without changing any
  // projected quantity. Restricting to {hole, particle} keeps the expansion
  // 2-way (occupied/virtual) per index -- exactly what the N/R classifier needs
  // -- instead of splitting every base space, which compounds across the nested
  // commutators. Falls back to all base spaces if the registry defines no
  // hole/particle split.
  const auto& hole_t = isr->hole_space(/*nulltype_ok=*/true);
  const auto& particle_t = isr->particle_space(/*nulltype_ok=*/true);
  auto physical = [&](const IndexSpace& b) {
    return hole_t.includes(b.type()) || particle_t.includes(b.type());
  };

  auto expand_term = [&](const ExprPtr& term) -> ExprPtr {
    // collect the residual NormalOperator's distinct general (non-base) indices
    const auto* nop = find_nop(term);
    if (!nop) return term;  // pure scalar/contraction: nothing to split
    container::svector<Index> gens;
    for (const auto& op : nop->creann()) {
      if (!is_base_space(op.index().space()) &&
          ranges::none_of(gens, [&](const auto& g) { return g == op.index(); }))
        gens.push_back(op.index());
    }
    if (gens.empty()) return term;
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

  ExprPtr out;
  if (expr.is<Sum>()) {
    auto out_sum = std::make_shared<Sum>();
    for (const auto& t : expr.as<Sum>()) out_sum->append(expand_term(t));
    out = out_sum->empty() ? ex<Constant>(0) : ExprPtr{out_sum};
  } else {
    // clone: expand_term may return its argument, which simplify would mutate
    out = expand_term(expr)->clone();
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
/// (N ⊎ R = expr as operators), R = expr − N holds exactly while keeping expr
/// in its compact (general-index) form -- only N is block-resolved. This is
/// verified equivalent to the fully block-resolved remainder: ref_av([R,σ]) is
/// identical either way. Keeping expr compact makes the nested commutators that
/// consume R operate on far fewer terms.
ExprPtr R_part(const ExprPtr& expr, std::size_t cutoff) {
  auto reduced = wick_reduce(expr->clone());
  return R_part_reduced(reduced, cutoff);
}

}  // namespace detail

/// Assembles H̄ order by order (see header), summing H̄⁰..H̄^rank of Eq. (45).
/// Each H̄^k below is transcribed from its equation, verified term-by-term
/// against the published coefficients and N/R subscripts. A subscript R/N on a
/// commutator means "form the commutator, then keep only its R/N part before
/// the next nesting".
ExprPtr hbar(std::size_t N, std::size_t rank, bool skip1) {
  if (rank > 4)
    throw Exception("bernoulli::hbar: only ranks 0..4 are implemented");

  using detail::N_part;
  using detail::N_part_reduced;
  using detail::R_part;
  using detail::R_part_reduced;
  using detail::wick_commutator;
  const auto cutoff = N;
  const auto F = op::tensor::F();
  const auto V = op::tensor::h(2);
  const auto T = op::tensor::T(N, skip1);
  const auto sigma = simplify(T - adjoint(T));  // σ = T − T†
  auto c = [&](rational num, const ExprPtr& e) {
    return ex<Constant>(num) * e;
  };

  // Every term of H̄^k is a nested commutator [[..[V_{p0},σ]_{f0}..],σ]_{f_k}
  // with a per-level N/R/A partition tag applied after each commutator ('A' =
  // no filter). nest(p0, f) evaluates such a node, memoizing every prefix
  // (key = p0 + tags applied so far): the terms share prefixes both within a
  // rank (the 9 rank-4 terms have only 3 distinct level-1 and 6 level-2 nodes)
  // and across ranks (each H̄^k node is a prefix of H̄^{k+1} nodes), so the memo
  // avoids recomputing them. Reusing a memoized ExprPtr is safe: expression
  // composition deep-copies operands (Product/Sum append clone), so
  // wick_commutator does not mutate its arguments. Commutator outputs are
  // already wick_reduce'd, so the reduced-input N/R filters apply.
  container::map<std::string, ExprPtr> memo;
  auto nest = [&](char p0, const char* f) -> ExprPtr {
    std::string key{p0};
    auto it = memo.find(key);
    if (it == memo.end()) {
      ExprPtr base = (p0 == 'N')   ? N_part(V, cutoff)
                     : (p0 == 'R') ? R_part(V, cutoff)
                                   : V;
      it = memo.emplace(std::move(key), std::move(base)).first;
    }
    ExprPtr op = it->second;
    for (int i = 0; f[i] != '\0'; ++i) {
      key = it->first + f[i];
      it = memo.find(key);
      if (it == memo.end()) {
        auto cx = wick_commutator(op, sigma);
        ExprPtr filtered = (f[i] == 'R')   ? R_part_reduced(cx, cutoff)
                           : (f[i] == 'N') ? N_part_reduced(cx, cutoff)
                                           : cx;
        it = memo.emplace(std::move(key), std::move(filtered)).first;
      }
      op = it->second;
    }
    return op;
  };

  // accumulate H̄ contributions via Sum::append (linear) rather than chained
  // operator+ / operator+=, each of which deep-copies the accumulated Sum --
  // quadratic in the (large) term count at high rank
  auto acc = std::make_shared<Sum>();
  acc->append(simplify(F + V));  // H̄⁰ = F + V   [Eq. (46)]
  if (rank >= 1) {
    // H̄¹ = [F,σ] + ½[V,σ] + ½[V_R,σ]   [Eq. (47)]. F enters H̄ ONLY here
    // (Cancellation #1, stated just below Eq. (50)).
    acc->append(wick_commutator(F, sigma));
    acc->append(c({1, 2}, nest('A', "A")));
    acc->append(c({1, 2}, nest('R', "A")));
  }
  if (rank >= 2) {
    // H̄² = 1/12[[V_N,σ],σ] + ¼[[V,σ]_R,σ] + ¼[[V_R,σ]_R,σ]   [Eq. (48)]
    acc->append(c({1, 12}, nest('N', "AA")));
    acc->append(c({1, 4}, nest('A', "RA")));
    acc->append(c({1, 4}, nest('R', "RA")));
  }
  if (rank >= 3) {
    // H̄³ = 1/24[[[V_N,σ],σ]_R,σ] + ⅛[[[V,σ]_R,σ]_R,σ] + ⅛[[[V_R,σ]_R,σ]_R,σ]
    //       − 1/24[[[V,σ]_R,σ],σ] − 1/24[[[V_R,σ]_R,σ],σ]   [Eq. (49)]
    acc->append(c({1, 24}, nest('N', "ARA")));
    acc->append(c({1, 8}, nest('A', "RRA")));
    acc->append(c({1, 8}, nest('R', "RRA")));
    acc->append(c({-1, 24}, nest('A', "RAA")));
    acc->append(c({-1, 24}, nest('R', "RAA")));
  }
  if (rank >= 4) {
    // H̄⁴ = Eq. (50), the nine order-4 terms produced by the recursion Eq. (44),
    // V̄^{k+1} = σ̂F + X̂⁻¹(σ̂)e^{σ̂}V − Σ_{n≠0} B_n σ̂^n V̄_R^{k}. F is absent here
    // (Cancellation #1). Listed in the paper's order; the outermost tag is
    // always A. Verified term-by-term against Eq. (50): coefficients
    // +1/16, +1/16, +1/48, −1/48, −1/48, −1/144, −1/48, −1/48, −1/720 and the
    // N/R subscripts all match as transcribed.
    acc->append(c({1, 16}, nest('R', "RRRA")));
    acc->append(c({1, 16}, nest('A', "RRRA")));
    acc->append(c({1, 48}, nest('N', "ARRA")));
    acc->append(c({-1, 48}, nest('A', "RARA")));
    acc->append(c({-1, 48}, nest('R', "RARA")));
    acc->append(c({-1, 144}, nest('N', "ARAA")));
    acc->append(c({-1, 48}, nest('A', "RRAA")));
    acc->append(c({-1, 48}, nest('R', "RRAA")));
    acc->append(c({-1, 720}, nest('N', "AAAA")));
  }
  ExprPtr result{std::move(acc)};
  return simplify(result);
}

}  // namespace sequant::mbpt::bernoulli
