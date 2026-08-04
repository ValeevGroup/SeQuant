//
// Per-term schedule-record emitter (JSON) for the batched-evaluation
// visualizer. Two producers share this schema:
//   * IR side (this header): schedule_ir_json() walks an annotated eval tree
//     and emits the factorizer's INTENDED schedule -- batched_here loops,
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

/// Canonical per-op signature for the avoidable-recompute JOIN: the result
/// index labels, then the two operands' label lists in a sorted
/// (order-independent) pair, so a binary contraction gets the SAME key
/// whichever operand binarize placed left. `full_label` (proto indices
/// included) keeps composite (proto-carrying) indices distinct. ONE definition,
/// two producers: the dry-run cost tally (\c DryRunOps::prod) keys its per-node
/// avoidable by it, and the runtime Build event (\c eval.hpp) plus the IR node
/// record below stamp each node with it -- so the visualizer maps a DAG node
/// (by hash) to cost_profile's per-node avoidable without reconstructing the
/// signature in the renderer, and without a second, drifting definition.
template <typename OutRange, typename LhsRange, typename RhsRange>
std::string cost_op_signature(OutRange const& out, LhsRange const& lhs,
                              RhsRange const& rhs) {
  auto join = [](auto const& r) {
    std::string s;
    for (auto const& ix : r) s += sequant::toUtf8(ix.full_label()) + ',';
    return s;
  };
  std::string const l = join(lhs);
  std::string const r = join(rhs);
  std::string sig = join(out);
  sig += '<';
  sig += (l <= r ? l : r);
  sig += '*';
  sig += (l <= r ? r : l);
  return sig;
}

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
     << ",\"batched_here\":[";
  {
    bool first = true;
    for (auto const& [ix, knd] : n->batched_here()) {
      if (!first) os << ',';
      first = false;
      os << "{\"index\":\""
         << sched_json_escape(sequant::toUtf8(ix.full_label()))
         << "\",\"mode\":\"" << sched_mode_name(knd) << "\"}";
    }
  }
  os << "],\"sliced_modes\":";
  sched_json_index_array(n->sliced_modes(), os);
  // sig = the avoidable-recompute join key (result + sorted operand pair). Only
  // internal (binary contraction) nodes have one; it matches the runtime Build
  // event's sig and cost_profile's per-node label, so the renderer maps this
  // node to cost_profile's avoidable count/exec by hash->sig.
  if (!n.leaf()) {
    os << ",\"sig\":\""
       << sched_json_escape(cost_op_signature(n->canon_indices(),
                                              n.left()->canon_indices(),
                                              n.right()->canon_indices()))
       << "\"";
  }
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
