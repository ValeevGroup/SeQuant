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
#include <SeQuant/core/optimize/options.hpp>

#include <cstddef>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
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
  void finish() const;

  /// Called by the optimizer: remembers that \p optimized is the summand it
  /// produced for \p key, so the Sum can be assembled in document order.
  /// Thread-safe.
  void bind(Expr const* optimized, std::string key);

  /// the key \ref bind associated with \p optimized, if any
  std::optional<std::string> key_of(Expr const* optimized) const;

 private:
  eval::schedule::RootRecord root_;
  std::map<std::string, std::size_t> position_;  // key -> first position
  std::vector<std::string> record_duplicates_;   // keys repeated in root_
  mutable std::mutex mutex_;
  std::vector<std::size_t> uses_;  // per record, times taken
  std::vector<std::string> missing_, duplicates_;
  std::unordered_map<Expr const*, std::string> keys_;
};

/// Collects the record_term results of an optimizer run and emits them as a
/// RootRecord: hand it to the optimizer via OptimizeOptions::recording.
/// \ref add and \ref relocate are thread-safe.
class RecordSession {
 public:
  explicit RecordSession(std::string root_label);

  RecordSession(RecordSession const&) = delete;
  RecordSession& operator=(RecordSession const&) = delete;

  /// thread-safe; keyed by the optimized summand's Expr const* (same key as
  /// term_batch_axes)
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
  std::vector<std::pair<Expr const*, eval::schedule::TermRecord>> records_;
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
/// A root that is not a Sum is recorded as a single term.
RecordedRoot record_root(std::string label, ExprPtr const& root_sum,
                         OptimizeOptions opts);

}  // namespace sequant::opt::schedule

#endif  // SEQUANT_CORE_OPTIMIZE_SCHEDULE_PIN_HPP
