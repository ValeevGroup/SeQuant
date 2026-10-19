//
// Symbolic half of schedule pinning: a context-independent key for a term, a
// record of the optimizer's output for that term (binary tree over factor
// ordinals plus per-node batch annotations), and the inverse that rebuilds
// the optimizer's output from such a record.
//

#ifndef SEQUANT_CORE_OPTIMIZE_SCHEDULE_PIN_HPP
#define SEQUANT_CORE_OPTIMIZE_SCHEDULE_PIN_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/schedule_document.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>

#include <string>
#include <vector>

namespace sequant::opt::schedule {

/// canonicalize(flatten(term->clone())); the object every other function here
/// works on. \p term itself is not modified.
ExprPtr canonical_term(ExprPtr const& term);

/// toUtf8(v1::to_string(canonical)) -- the identity of a term in the document
std::string term_key(ExprPtr const& canonical);

/// factors of a canonical Product in order (a lone Tensor is one factor; a
/// Constant is none)
std::vector<ExprPtr> canonical_factors(ExprPtr const& canonical);

/// Records the optimizer's output \p optimized for \p canonical (a nested
/// Product over the same factor objects, or the lone Tensor/Constant itself)
/// and its annotations \p nodes (one per tensor x tensor node, left-first
/// post-order, as consumed by \c binarize).
/// \throws sequant::Exception if \p canonical is a mixed product (has a
///         non-Tensor factor), or if \p optimized / \p nodes do not match it
eval::schedule::TermRecord record_term(
    ExprPtr const& canonical, ExprPtr const& optimized,
    container::vector<NodeBatchAnnotation> const& nodes);

/// the optimizer's output rebuilt from a TermRecord
struct PinnedTerm {
  ExprPtr expr;
  container::vector<NodeBatchAnnotation> nodes;
};

/// Rebuild the optimizer's output from a record over \p canonical's factors;
/// the annotations' Index objects are looked up by full_label among the
/// subtree's indices.
/// \throws sequant::Exception naming `rec.key` when the factors differ from
///         \p canonical's, or an axis label is not an index of its node's
///         subtree
PinnedTerm rebuild_pinned_term(eval::schedule::TermRecord const& rec,
                               ExprPtr const& canonical);

}  // namespace sequant::opt::schedule

#endif  // SEQUANT_CORE_OPTIMIZE_SCHEDULE_PIN_HPP
