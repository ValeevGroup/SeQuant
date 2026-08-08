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

#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/optimize/ga/fitness.hpp>

#include <string_view>

namespace sequant::opt::ga {

/// Label prefix of named GA intermediates ("IGA1", "IGA2", ...).
inline constexpr std::wstring_view named_intermediate_prefix = L"IGA";

/// Whether \p label names a GA intermediate (i.e. starts with
/// named_intermediate_prefix followed by digits).
inline bool is_named_intermediate(std::wstring_view label) {
  auto const& p = named_intermediate_prefix;
  if (label.size() <= p.size() || label.substr(0, p.size()) != p) return false;
  for (auto c : label.substr(p.size()))
    if (c < L'0' || c > L'9') return false;
  return true;
}

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

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_EMIT_HPP
