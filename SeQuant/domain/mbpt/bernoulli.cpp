#include <SeQuant/domain/mbpt/bernoulli.hpp>

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/wick.hpp>

// Bernoulli expansion of the unitary-CC similarity-transformed Hamiltonian
// H̄ = e^{−σ} H e^{σ}, σ = T − T† (anti-Hermitian). Because σ mixes excitation
// and de-excitation the plain BCH series does not terminate; the Bernoulli
// expansion rewrites it so that Bernoulli numbers are the expansion
// coefficients, leaving the final truncation at a chosen commutator rank as the
// only approximation.
//
// This file builds that expansion in three layers: an operator-valued Wick
// reduction (here), the N/R operator split, and the rank-by-rank assembly of
// H̄. All equation references are to 10.1063/1.5030344 (Sec. III B).

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

}  // namespace detail

}  // namespace sequant::mbpt::bernoulli
