#ifndef SEQUANT_CORE_OPTIMIZE_GA_EMIT_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_EMIT_HPP

// Emission: turn a winning Schedule into evaluable expressions.
//
// Shared intermediates are SYNTACTIC, not structural: every keyed array the
// schedule uses at more than one site becomes a NAMED tensor (label "IGA<n>",
// slots = the producing cluster's canonical face order). Its definition is
// emitted exactly once, in dependency order, and every use site (including the
// producer's own) is a leaf tensor IGA<n>{...} whose slots are the use
// cluster's canonical face order. Since key-equal clusters correspond
// axis-wise by zipping their canonical faces (Fitness::correspondences,
// identity automorphism), slot position k of every use leaf denotes the same
// array axis as slot position k of the definition head: the runtime resolves a
// use leaf purely POSITIONALLY and never has to re-derive a correspondence
// between differently-labeled trees.
//
// Rationale: the label-blind structural matching of the runtime eval cache is
// only sound when structural identity implies value identity -- which per-term
// canonical naming guarantees, but fresh-labeled unrolled producer subtrees do
// not. Named intermediates make the sharing explicit in the expression itself.
//
// Keys used at a single site are still inlined (the picked producer's subtree
// with the face renamed through the canonical correspondence); factored sums
// are emitted as (X1' + X2') * V. inline_definitions substitutes the
// definitions back, reproducing the fully-unrolled emission exactly.
//
// WHICH keys get named is not this file's decision: it is
// `CostModel::runtime_amortized`, the same function `Fitness::resolve` charges
// by. Naming a key the objective still charges per replay wastes residency;
// charging once a key emission inlines makes the objective lie. Under
// `CostModel::amortize_persistent` the named class widens from "used at >= 2
// sites" to "used at >= 2 sites OR amplitude-free, slot-bearing and under
// `naming_cap_elems`", on both sides at once.

#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/optimize/ga/fitness.hpp>
// re-exported: named_intermediate_prefix / is_named_intermediate. Consumers
// that need ONLY the predicate should include that header instead -- it is
// standalone, this one is not.
#include <SeQuant/core/optimize/ga/intermediate_label.hpp>

namespace sequant::opt::ga {

struct Emission {
  /// One factorized expression per target, externals named like the target's
  /// first term; shared arrays appear as named leaf tensors (IGA<n>).
  container::svector<ExprPtr> targets;
  /// Definitions of the named intermediates, in dependency (topological)
  /// order. The head tensor's aux slots are the producing cluster's canonical
  /// face order -- the positional layout every use leaf's slots correspond to.
  container::svector<ResultExpr> definitions;
};

/// Named-intermediate emission of \p schedule (see file comment).
Emission emit_named(Fitness const& fitness, Schedule const& schedule);

/// Substitute every named intermediate of \p em back into the targets
/// (bodies inlined recursively, the face renamed positionally, internal
/// indices freshened). The fully-unrolled form; used to verify equivalence.
container::svector<ExprPtr> inline_definitions(Emission const& em);

/// Fully-inlined per-target expressions:
/// inline_definitions(emit_named(fitness, schedule)).
container::svector<ExprPtr> emit(Fitness const& fitness,
                                 Schedule const& schedule);

/// Render-site count per keyed array -- `emit_named`'s pass 1, run on its own.
/// Exposed so the [ga] tests can pin it against `Schedule::uses`: emission
/// counts render sites and `Fitness::resolve` counts L2 demands plus parent
/// slots, and the matched pair is only matched while those two counts agree.
container::map<std::size_t, std::size_t> emission_use_counts(
    Fitness const& fitness, Schedule const& schedule);

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_EMIT_HPP
