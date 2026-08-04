#ifndef SEQUANT_EVAL_OCCURRENCE_KEY_HPP
#define SEQUANT_EVAL_OCCURRENCE_KEY_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <algorithm>
#include <cstddef>

namespace sequant::eval {

/// \brief Batched modes (proto-expanded) that actually appear on \p node's own
///        slots, drawn from the ambient batch-loop modes \p ctx_modes.
///
/// \details A mode \c m in \p ctx_modes is proto-expanded first: if
/// \c m.has_proto_indices() its protoindices stand in for it (a composite
/// batched mode is domain-tied to its proto pair), else \c m stands for
/// itself. An expanded mode is included in the result iff it also appears
/// among \p node's own canonical slots (\c node->canon_indices()),
/// proto-expanded the same way on the slot side -- i.e. iff the mode actually
/// lives on \p node's result. This mirrors \c slot_modes_of /
/// \c ext_modes_of in \c lifetime_mask.hpp exactly (same proto-expansion on
/// both the ambient-mode side and the node-slot side).
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
/// \return the subset of \p ctx_modes (proto-expanded) that lives on \p
///         node's own slots.
template <meta::eval_node Node>
TensorNetwork::NamedIndexSet in_scope_batched_on_node(
    Node const& node, container::svector<Index> const& ctx_modes) {
  // proto-expand the ambient batch-loop modes
  container::svector<Index> ctx_expanded;
  for (auto const& m : ctx_modes) {
    if (m.has_proto_indices())
      for (auto const& p : m.proto_indices()) ctx_expanded.push_back(p);
    else
      ctx_expanded.push_back(m);
  }

  // proto-expand node's own canonical slots
  container::svector<Index> slots;
  for (auto const& s : node->canon_indices()) {
    if (s.has_proto_indices())
      for (auto const& p : s.proto_indices()) slots.push_back(p);
    else
      slots.push_back(s);
  }

  TensorNetwork::NamedIndexSet named;
  for (auto const& ix : ctx_expanded)
    if (std::find(slots.begin(), slots.end(), ix) != slots.end())
      named.insert(ix);
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
/// its slots (\c TensorNetwork::canonicalize_slots): the in-scope batched
/// indices are passed as \c named_indices, so they are colored by space +
/// label (cannot be renamed, their labels/positions are meaningful) while
/// every other index is colored by space alone (freely renamed,
/// interchangeable). Symmetry of the resulting key is thus INTRINSIC to the
/// canonicalized colored graph -- no separately-deduced symmetry annotation
/// is needed. This mirrors \c build_subnet_metadata
/// (\c single_term_detail.hpp) exactly, applied to a node's own leaves rather
/// than a DP subset.
///
/// \param node the node whose sub-expression is keyed.
/// \param ctx_modes the ambient batch-loop modes in scope at the call site.
/// \return the canonicalization metadata; compare with \c RouterKeyEqual,
///         hash with \c RouterKeyHash.
template <meta::eval_node Node>
TensorNetwork::SlotCanonicalizationMetadata occurrence_key(
    Node const& node, container::svector<Index> const& ctx_modes) {
  container::vector<ExprPtr> leaves;
  auto collect = [&](auto&& self, Node const& n) -> void {
    if (n.leaf()) {
      if (n->is_tensor()) leaves.emplace_back(n->as_tensor().clone());
      return;
    }
    self(self, n.left());
    self(self, n.right());
  };
  collect(collect, node);

  auto named = in_scope_batched_on_node(node, ctx_modes);
  auto tn = TensorNetwork{leaves};
  return tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels(),
                               &named);
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
