//
// Symbolic half of schedule pinning: a context-independent key for a term, a
// record of the optimizer's output for that term (binary tree over factor
// ordinals plus per-node batch annotations), and the inverse that rebuilds
// the optimizer's output from such a record.
//

#ifndef SEQUANT_CORE_OPTIMIZE_SCHEDULE_PIN_HPP
#define SEQUANT_CORE_OPTIMIZE_SCHEDULE_PIN_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/schedule_document.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstddef>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sequant::opt::schedule {

/// canonicalize(flatten(term->clone())); the object every other function here
/// works on. \p term itself is not modified.
ExprPtr canonical_term(ExprPtr const& term);

/// v1::to_string(canonical) -- the identity of a term in the document.
/// Canonicalization renames the dummy (contracted) indices only: the external
/// indices keep the caller's labels and are part of the key, so a root must be
/// spelled with the same external labels on every run that pins it. The key
/// also depends on the active Context's canonicalization settings.
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

/// The record of a term the optimizer handles but a pin cannot rebuild (a
/// factor is a Sum or a Variable): key, factors and scalar of \p canonical,
/// no tree, `pinnable` false. Recording lists such a term for completeness;
/// rebuild_pinned_term refuses it.
eval::schedule::TermRecord record_unpinnable_term(ExprPtr const& canonical);

/// the optimizer's output rebuilt from a TermRecord
struct PinnedTerm {
  ExprPtr expr;
  container::vector<NodeBatchAnnotation> nodes;
};

/// Rebuild the optimizer's output from a record over \p canonical's factors;
/// the annotations' Index objects are looked up by
/// eval::schedule::index_label among the subtree's indices (protos included).
/// \throws sequant::Exception naming `rec.key` when the factors differ from
///         \p canonical's, or an axis label is not an index of its node's
///         subtree
PinnedTerm rebuild_pinned_term(eval::schedule::TermRecord const& rec,
                               ExprPtr const& canonical);

/// Replays one RootRecord through \ref sequant::optimize_result: hand it to
/// the optimizer via OptimizeOptions::pinned. Every summand looks its record
/// up by key; a miss is remembered (and the summand left unoptimized), and
/// \ref finish reports every miss, unused record and duplicate key at once.
/// \ref take and \ref bind are thread-safe (summands are optimized in
/// parallel).
class PinSession {
 public:
  explicit PinSession(eval::schedule::RootRecord const& root);

  PinSession(PinSession const&) = delete;
  PinSession& operator=(PinSession const&) = delete;

  /// thread-safe; returns nullptr on a miss (recorded) and records a duplicate
  /// hit (a key taken more than once)
  eval::schedule::TermRecord const* take(std::string const& key);

  /// position of \p key in the document order, for callers assembling their
  /// own Sum
  std::optional<std::size_t> document_position(std::string const& key) const;

  /// \throws sequant::Exception listing "missing: ...", "unused: ...",
  ///         "duplicate: ..." (all keys) if any key was missed, any record was
  ///         not taken, or any key was taken twice (or occurs twice in the
  ///         record)
  /// \post the session is finished: a session is single-use, and \ref take
  ///       and \ref bind throw from now on (calling finish again is fine)
  void finish() const;

  /// Called by the optimizer: remembers that \p optimized is the summand it
  /// produced for \p key, so the Sum can be assembled in document order.
  /// A repeated address (freed by a per-summand caller and reused) is
  /// overwritten: the newest binding wins. Thread-safe.
  void bind(Expr const* optimized, std::string key);

  /// the key \ref bind associated with \p optimized, if any. Only meaningful
  /// while that result is alive (the optimizer reads it only for the summands
  /// of the whole Sum it is assembling); a freed address may since have been
  /// reused.
  std::optional<std::string> key_of(Expr const* optimized) const;

 private:
  eval::schedule::RootRecord root_;
  std::map<std::string, std::size_t> position_;  // key -> first position
  std::vector<std::string> record_duplicates_;   // keys repeated in root_
  mutable std::mutex mutex_;
  std::vector<std::size_t> uses_;  // per record, times taken
  std::vector<std::string> missing_, duplicates_;
  std::unordered_map<Expr const*, std::string> keys_;
  mutable bool finished_ = false;
};

/// Collects the record_term results of an optimizer run and emits them as a
/// RootRecord: hand it to the optimizer via OptimizeOptions::recording.
/// \ref add and \ref relocate are thread-safe.
///
/// Records are kept in call order and never dropped. Resolution by pointer
/// (\ref finish(std::vector<Expr const*> const&) and \ref relocate) requires
/// the optimized results to be alive, since a freed address can be reused by a
/// later result (the newest record for an address wins); call-order
/// \ref finish() does not, so per-summand callers may drop each result.
class RecordSession {
 public:
  explicit RecordSession(std::string root_label);

  RecordSession(RecordSession const&) = delete;
  RecordSession& operator=(RecordSession const&) = delete;

  /// thread-safe; appends \p rec, and keys it by the optimized summand's
  /// Expr const* (same key as term_batch_axes) for pointer resolution
  void add(Expr const* optimized, eval::schedule::TermRecord rec);

  /// records in \p order (final Sum summand pointers)
  /// \throws sequant::Exception on a pointer with no record, or if two
  ///         records share a key
  eval::schedule::RootRecord finish(
      std::vector<Expr const*> const& order) const;

  /// records in call order (per-summand callers)
  /// \throws sequant::Exception if two records share a key
  eval::schedule::RootRecord finish() const;

  /// Called by the optimizer when it replaces summands by copies (opt::reorder
  /// clones each summand into the final Sum): re-keys the record of each
  /// `first` pointer onto the matching `second` pointer. Pointers with no
  /// record are skipped. Thread-safe.
  void relocate(std::vector<std::pair<Expr const*, Expr const*>> const& moves);

 private:
  std::string label_;
  mutable std::mutex mutex_;
  std::vector<eval::schedule::TermRecord> records_;     // in call order
  std::unordered_map<Expr const*, std::size_t> index_;  // into records_
};

/// a recorded root: its record and the optimized expression it was taken from
struct RecordedRoot {
  eval::schedule::RootRecord record;
  ExprPtr optimized;
};

/// Runs \ref sequant::optimize_result on \p root_sum with a RecordSession
/// labeled \p label set as \p opts.recording, and returns the record (terms in
/// the order of the optimized Sum's summands) with the optimized expression.
/// A root that is not a Sum is recorded as a single term. The Sum must be
/// simplified: no two summands may have the same canonical key.
/// \throws sequant::Exception listing the duplicate keys otherwise
RecordedRoot record_root(std::string label, ExprPtr const& root_sum,
                         OptimizeOptions opts);

namespace detail {

/// the non-template part of \ref fingerprint_of: \p representatives has one
/// Index per IndexBasis::base_key() present in the forest, \p
/// volatile_leaf_labels the labels of the volatile leaf tensors
eval::schedule::Fingerprint make_fingerprint(
    std::vector<std::string> root_labels, OptimizeOptions const& opts,
    std::vector<std::string> const& mode_order,
    std::map<std::string, Index> const& representatives,
    std::set<std::string> const& volatile_leaf_labels);

}  // namespace detail

/// Everything beyond the per-term trees that shapes the fused DAG: the
/// batchability predicates evaluated on one representative Index per basis
/// present in \p forest (keyed by IndexBasis::base_key(), the key the
/// engine's loops and the DAG record's mode spellings carry, which is the
/// basis name when the basis is named), the volatile leaf labels, the scalar
/// BatchPolicy and OptimizeOptions settings that steer structure, \p mode_order
/// and the active Context (vacuum, basis, strict bra-ket symmetry, cardinal
/// tensor labels, canonicalization method and the named-index flag, the
/// settings a term key depends on). Size-dependent settings are deliberately
/// absent (see \ref planner_record).
///
/// Two limits: a predicate is evaluated on one index per basis, so a policy
/// that also reads the index itself (its proto indices, an attribute that
/// varies within a basis) is fingerprinted by that one index alone and a
/// difference there surfaces as a different DAG rather than a mismatch; and
/// the tensor canonicalizers and index comparers the Context registers are
/// not spelled, so a difference there surfaces as a key miss, which
/// PinSession::finish names.
template <meta::eval_node_range R>
eval::schedule::Fingerprint fingerprint_of(
    R const& forest, std::vector<std::string> root_labels,
    OptimizeOptions const& opts, std::vector<std::string> const& mode_order) {
  std::map<std::string, Index> representatives;
  std::set<std::string> volatile_leaf_labels;
  auto const& is_volatile = opts.batch_policy.is_volatile_leaf;
  auto walk = [&](auto const& self, auto const& node) -> void {
    if (!node.leaf()) {
      self(self, node.left());
      self(self, node.right());
      return;
    }
    if (!node->is_tensor()) return;
    for (auto const& ix : node->canon_indices())
      representatives.try_emplace(ix.basis().base_key(), ix);
    auto const& t = node->as_tensor();
    if (is_volatile && is_volatile(t))
      volatile_leaf_labels.insert(std::string(t.label()));
  };
  for (auto const& root : forest) walk(walk, root);
  return detail::make_fingerprint(std::move(root_labels), opts, mode_order,
                                  representatives, volatile_leaf_labels);
}

/// \throws sequant::Exception listing every field in which \p loaded and
///         \p live differ, one "fingerprint.<field>: <a> vs <b>" per line
void check_fingerprint(eval::schedule::Fingerprint const& loaded,
                       eval::schedule::Fingerprint const& live);

/// the optimizer settings the fingerprint excludes (they only steer the
/// single-term optimizer), stringified; for the document's run.planner
std::map<std::string, std::string> planner_record(OptimizeOptions const& opts);

}  // namespace sequant::opt::schedule

#endif  // SEQUANT_CORE_OPTIMIZE_SCHEDULE_PIN_HPP
