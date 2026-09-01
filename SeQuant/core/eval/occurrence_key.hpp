#ifndef SEQUANT_EVAL_OCCURRENCE_KEY_HPP
#define SEQUANT_EVAL_OCCURRENCE_KEY_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network/canonicals.hpp>
#include <SeQuant/core/tensor_network/typedefs.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <algorithm>
#include <cstddef>

namespace sequant::eval {

/// \brief The ambient batch-loop modes \p ctx_modes that actually appear on
///        \p node's own result slots.
///
/// \details A mode \c m in \p ctx_modes is included iff it appears among
/// \p node's own slots (\c detail::slot_modes_of, i.e. \c canon_indices() taken
/// as-is). Both sides are plain modes: batch loops are always over plain
/// occ/aux modes, and a composite slot \c a<i,j> contributes only mode \c a
/// (see
/// \c slot_modes_of for why the \c i,j proto pair is NOT exposed). This is the
/// same "which in-scope batch modes live on this node's result" filter the
/// residency meet applies in \c stamp_lifetime_masks / \c stamp_residency_impl.
/// The kind (\c BatchModeType) is never inspected.
///
/// Site-invariant: the result is identical whether computed at a read site
/// (a consumer's \c cache.batch_context()) or a store site (the producer's
/// \c ectx), because it depends only on \p node's own slots and the modes
/// that are actually in scope there -- not on how deep in the batch-context
/// stack those modes sit.
///
/// \param node the (sub-expression) node whose in-scope batched modes are
///        wanted.
/// \param ctx_modes the ambient batch-loop modes in scope at the call site.
/// \return the subset of \p ctx_modes that lives on \p node's own slots.
template <meta::eval_node Node>
TensorNetwork::NamedIndexSet in_scope_batched_on_node(
    Node const& node, container::svector<Index> const& ctx_modes) {
  // node's own result slots. Fully qualified: this file's enclosing namespace
  // is sequant::eval, which (via schedule_dump.hpp, in some inclusion orders)
  // can ALSO have its own sequant::eval::detail -- an unqualified `detail::`
  // here would silently resolve to the wrong one.
  auto const slots = sequant::detail::slot_modes_of(node);

  TensorNetwork::NamedIndexSet named;
  for (auto const& m : ctx_modes)
    if (std::find(slots.begin(), slots.end(), m) != slots.end())
      named.insert(m);
  return named;
}

/// \brief The router occurrence key: canonicalizes \p node's actual
///        sub-expression as a colored \c TensorNetwork with the in-scope
///        batched indices (\c in_scope_batched_on_node) as \c named_indices.
///
/// \details An "occurrence" is a use-site: a place a value is read. Two reads
/// are the SAME occurrence when they read the same value with the same
/// batched-slot structure, up to the value's symmetries. This is computed by
/// treating \p node's leaf tensors as a \c TensorNetwork and canonicalizing
/// its slots (\c TensorNetwork::canonicalize_slots, which colors with
/// \c distinct_named_indices=false): the in-scope batched indices are passed as
/// \c named_indices, so they are colored by space and by their named-ness --
/// distinct from non-batched (anonymous) indices of the same space, and from
/// batched indices of a different space -- but NOT by their specific label
/// (named indices of one space are interchangeable within that class), while
/// every other index is colored by space alone. \c RouterKeyEqual compares only
/// the canonical colored graph, so the key is DAG-GLOBAL: two occurrences of
/// one value whose batched slot binds different physical labels (the g.C legs,
/// i_3 vs i_4) map to the SAME key -- one overlay serves both, and the specific
/// binding is resolved later in home-scope computation (see the design spec
/// doc/dev/specs/2026-08-07-remat-cse-aware-split-design.md, and the
/// "same batched slot, different batched label" test in test_occurrence_key).
/// Symmetry of the resulting key is thus INTRINSIC to the canonicalized colored
/// graph -- no separately-deduced symmetry annotation is needed. This mirrors
/// \c build_subnet_metadata
/// (\c single_term_detail.hpp) exactly, applied to a node's own leaves rather
/// than a DP subset.
///
/// \par Loop coloring (sliced-value canonical-layout design, Task 5). When
/// \p named_index_colors is non-null the in-scope batched (named) indices are
/// colored not merely by space but by the DAG-scope LOOP that slices each --
/// see \c tensor_network::NamedIndexColorMap and \c loop_colored_id
/// (ordered_schedule.hpp), which builds the map from
/// \c compute_sliced_mode_assignment. Two same-space named slots bound to
/// DIFFERENT loops then receive DIFFERENT graph colors (no longer
/// interchangeable), and same-loop bindings still fold. A null (the default,
/// and every unsliced/production consult) or empty map leaves the key
/// BYTE-IDENTICAL to today's space-only named coloring -- the empty-sliced-set
/// reduction the design pins as its #1 non-regression anchor: the null path
/// below is literally the pre-coloring two-argument \c canonicalize_slots call.
///
/// \param node the node whose sub-expression is keyed.
/// \param ctx_modes the ambient batch-loop modes in scope at the call site.
/// \param named_index_colors optional per-named-index loop-color map (default
///        null => today's space-only coloring, byte-identical).
/// \return the canonicalization metadata; compare with \c RouterKeyEqual,
///         hash with \c RouterKeyHash.
template <meta::eval_node Node>
TensorNetwork::SlotCanonicalizationMetadata occurrence_key(
    Node const& node, container::svector<Index> const& ctx_modes,
    tensor_network::NamedIndexColorMap const* named_index_colors = nullptr) {
  container::vector<ExprPtr> leaves;
  auto collect = [&](auto&& self, Node const& n) -> void {
    if (n.leaf()) {
      if (n->is_tensor()) leaves.emplace_back(n->as_tensor().clone());
      return;
    }
    // PRECONDITION: occurrence_key flattens the subtree's leaves into a single
    // TensorNetwork, which is only well-defined for a tensorial (contraction)
    // subtree. A Sum node -- anywhere in the subtree -- is NOT a tensor
    // network: its summands reuse the same dummy labels independently, so
    // unioning their leaves glues unrelated terms together (an index would
    // connect to >1 bra/ket slot, which create_graph then rejects with a
    // cryptic strict-braket abort). A Sum-bearing occurrence needs a structural
    // (op-over-child-identities) key, not a flattened TN (see the hybrid design
    // discussion); until that exists, fail fast with a clear diagnostic here
    // rather than deep inside create_graph.
    SEQUANT_ASSERT(
        !n->is_sum() &&
        "occurrence_key: cannot key a subtree containing a Sum node -- a Sum "
        "is "
        "not a tensor network; only tensorial (contraction) subtrees have a "
        "well-defined occurrence key");
    self(self, n.left());
    self(self, n.right());
  };
  collect(collect, node);

  auto named = in_scope_batched_on_node(node, ctx_modes);
  auto tn = TensorNetwork{leaves};
  // REDUCTION anchor: with no loop colors this is exactly the pre-coloring
  // two-argument call (default named_index_compare, no color map) -- keep it a
  // distinct branch so the null path stays byte-identical rather than relying
  // on `{}` matching the default comparator (it does NOT: an empty
  // named_index_compare falls back to space-only, whereas the two-arg default
  // is default_idxptr_slottype_lesscompare, which also orders by proto-arity
  // and slot type -- see canonicals.hpp).
  if (!named_index_colors)
    return tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels(),
                                 &named);
  return tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels(),
                               &named, default_idxptr_slottype_lesscompare{},
                               named_index_colors);
}

/// \brief Hasher for \c TensorNetwork::SlotCanonicalizationMetadata, for use
///        as a router occurrence-key map's hash functor.
struct RouterKeyHash {
  std::size_t operator()(
      TensorNetwork::SlotCanonicalizationMetadata const& m) const noexcept {
    return m.hash_value();
  }
};

/// \brief Equality for \c TensorNetwork::SlotCanonicalizationMetadata, for
///        use as a router occurrence-key map's equality functor.
struct RouterKeyEqual {
  bool operator()(TensorNetwork::SlotCanonicalizationMetadata const& a,
                  TensorNetwork::SlotCanonicalizationMetadata const& b) const {
    return bliss::ConstGraphCmp::cmp(*a.graph, *b.graph) == 0;
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_OCCURRENCE_KEY_HPP
