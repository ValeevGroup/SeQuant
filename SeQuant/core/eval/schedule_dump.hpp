//
// Per-term schedule-record emitter (JSON) for the batched-evaluation
// visualizer. Two producers share this schema:
//   * IR side (this header): schedule_ir_json() walks an annotated eval tree
//     and emits the factorizer's INTENDED schedule -- node_slice_mask loops,
//     sliced/contracted modes, order-aware gate, effective (replay) count.
//   * Runtime side (eval.hpp): the evaluator emits the SAME node-keyed schema
//     with the ACTUAL lifetimes (life/max_life), hoist placement, replay, and
//     footprint, so the two can be diffed node-by-node.
//
// The schema is deliberately minimal JSON (no external dependency): a per-term
// object { term_id, kind, root } whose nodes carry the fields below and a
// `children` array (binary: [left,right] for the IR tree; runtime nodes may
// carry richer children). Index labels and expr identities are emitted RAW
// (SeQuant full_label / serialization); the renderer formats them.
//

#ifndef SEQUANT_CORE_EVAL_SCHEDULE_DUMP_HPP
#define SEQUANT_CORE_EVAL_SCHEDULE_DUMP_HPP

#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <ostream>
#include <sstream>
#include <string>

namespace sequant::eval {

namespace detail {

/// DIAGNOSTIC: the flops the LAST product op computed (cm->flops), stashed by
/// DryRunOps::prod and read by the Build-event emitter in eval.hpp so each
/// per-build event can carry the SAME cost the replay tallied -- letting a
/// caller sum per DISTINCT node (by hash) with the replay's own flop method,
/// no second flop path. Thread-local; single-threaded dry-run replay.
inline double& last_op_flops() noexcept {
  static thread_local double f = 0.0;
  return f;
}

/// DIAGNOSTIC: the roofline exec-cost estimate (cm->exec_cost) for the LAST
/// product op, stashed by DryRunOps::prod alongside \c last_op_flops() and
/// read by the Build-event choke in eval.hpp so each per-build tally entry
/// carries the SAME time estimate the replay computed. Thread-local;
/// single-threaded dry-run replay.
inline double& last_op_exec() noexcept {
  static thread_local double f = 0.0;
  return f;
}

inline std::string sched_json_escape(std::string const& s) {
  std::string o;
  o.reserve(s.size() + 8);
  for (char c : s) {
    if (c == '"' || c == '\\') {
      o += '\\';
      o += c;
    } else if (c == '\n') {
      o += "\\n";
    } else {
      o += c;
    }
  }
  return o;
}

inline char const* sched_mode_name(BatchModeType m) {
  return m == BatchModeType::External ? "External" : "Contracted";
}

template <typename Range>
void sched_json_index_array(Range const& r, std::ostream& os) {
  os << '[';
  bool first = true;
  for (auto const& ix : r) {
    if (!first) os << ',';
    first = false;
    os << '"' << sched_json_escape(sequant::toUtf8(ix.full_label())) << '"';
  }
  os << ']';
}

}  // namespace detail

// NOTE: the former `cost_op_signature` (a full-label string keying the
// avoidable-recompute rollup) was REMOVED: a signature is dummy-label- and
// batch-slice-dependent, so per-block and alpha-renamed builds of one value
// split across keys and the rollup mis-attributed them. Nodes are now joined to
// cost_profile's per-node recompute by their TOPOLOGICAL hash (\c
// EvalExpr::hash_value), the same identity the CacheManager dedups on.

/// Emit one node of the IR schedule record. \p n is a FullBinaryNode whose
/// value type exposes the batch-annotation accessors (EvalExpr / EvalExprTA).
template <typename Node>
void schedule_ir_json_node(Node const& n, std::ostream& os) {
  using detail::sched_json_escape;
  using detail::sched_json_index_array;
  using detail::sched_mode_name;

  std::string ident;
  try {
    ident = sequant::toUtf8(sequant::io::serialization::to_string(n->expr()));
  } catch (...) {
    ident.clear();
  }

  // hash = the cache/DAG identity: canonically-equal subexpressions share it,
  // so the renderer collapses the forest into the DAG the CacheManager sees
  // (and the runtime overlay joins by the same key). Emitted as a string to
  // survive 64-bit precision through JSON/JS.
  os << "{\"hash\":\"" << n->hash_value() << "\",\"expr\":\""
     << sched_json_escape(ident) << "\",\"result\":";
  sched_json_index_array(n->canon_indices(), os);
  os << ",\"leaf\":" << (n.leaf() ? "true" : "false")
     << ",\"order_aware\":" << (n->batch_order_aware() ? "true" : "false")
     << ",\"effective_count\":" << n->batch_effective_count()
     << ",\"node_slice_mask\":[";
  {
    bool first = true;
    for (auto const& [ix, knd] : n->node_slice_mask()) {
      if (!first) os << ',';
      first = false;
      os << "{\"index\":\""
         << sched_json_escape(sequant::toUtf8(ix.full_label()))
         << "\",\"mode\":\"" << sched_mode_name(knd) << "\"}";
    }
  }
  os << "],\"sliced_modes\":";
  sched_json_index_array(n->sliced_modes(), os);
  // Nodes are joined to cost_profile's per-node recompute by their topological
  // hash (emitted above), not by a signature: a signature is dummy-label- and
  // batch-slice-dependent and mis-folds per-block / alpha-renamed builds.
  os << ",\"children\":[";
  if (!n.leaf()) {
    schedule_ir_json_node(n.left(), os);
    os << ',';
    schedule_ir_json_node(n.right(), os);
  }
  os << "]}";
}

/// Emit the per-term IR schedule record as a single JSON object.
template <typename Node>
std::string schedule_ir_json(Node const& root, std::string const& term_id) {
  std::ostringstream os;
  os << "{\"term_id\":\"" << detail::sched_json_escape(term_id)
     << "\",\"kind\":\"ir\",\"root\":";
  schedule_ir_json_node(root, os);
  os << '}';
  return os.str();
}

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_SCHEDULE_DUMP_HPP
