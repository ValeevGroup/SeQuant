//
// The read-only DAG section of the schedule document: the fused values,
// loop blocks, loop instances and operand edges of a batched-evaluation
// schedule, recorded from the (RichSchedule, LegalitySchedule,
// OrderedSchedule) triple with run-independent value ids.
//

#ifndef SEQUANT_EVAL_SCHEDULE_DAG_RECORD_HPP
#define SEQUANT_EVAL_SCHEDULE_DAG_RECORD_HPP

#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/eval/cell_table_builder.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/schedule_document.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace sequant::eval::schedule {

/// "<label>:<base_key>" spelling used for every mode list in the DAG record
inline std::string mode_spelling(Index const& ix) {
  return index_label(ix) + ":" + toUtf8(ix.basis().base_key());
}

namespace detail {

inline std::string hex_string(std::size_t v) {
  char buf[2 + 2 * sizeof(std::size_t) + 1];
  std::snprintf(buf, sizeof buf, "0x%zx", v);
  return buf;
}

template <typename Indices>
std::vector<std::string> spell_modes(Indices const& indices) {
  std::vector<std::string> out;
  for (Index const& ix : indices) out.push_back(mode_spelling(ix));
  return out;
}

inline std::string kind_name(BatchModeType k) {
  return k == BatchModeType::External ? "external" : "contracted";
}

inline std::string role_name(LoopRole r) {
  switch (r) {
    case LoopRole::LoopLocal:
      return "LoopLocal";
    case LoopRole::Reduction:
      return "Reduction";
    case LoopRole::LoopCarried:
      return "LoopCarried";
    case LoopRole::LoopInvariant:
      return "LoopInvariant";
  }
  throw Exception("record_dag: unknown LoopRole");
}

inline std::string output_kind_name(OutputKind k) {
  switch (k) {
    case OutputKind::Transient:
      return "Transient";
    case OutputKind::AccumulateSum:
      return "AccumulateSum";
    case OutputKind::AccumulateScatter:
      return "AccumulateScatter";
  }
  throw Exception("record_dag: unknown OutputKind");
}

/// Where a forest node sits in the document's naming: the occurrence's
/// `term` and `node` and the id it spells.
///  - a node of a term's pinned tree: (term key, left-first post-order
///    ordinal over leaves and products), id "<key>#<ordinal>";
///  - the node binarize adds to apply a term's scalar: (term key, number of
///    tree positions + 1 -- it is the term's top, after every tree node and
///    the scalar), id "<key>#scale"; its right operand, the Constant leaf
///    holding the scalar, is (term key, number of tree positions), id
///    "<key>#scalar" (spec 6.1 names no id for it; it is a DAG value too);
///  - an Adjoint node (a conjugated factor) takes one tree ordinal as a leaf
///    would; its bare-label operand is id "<key>#<ordinal>/0" and its
///    Constant{1} sentinel "<key>#<ordinal>/1", both at the same (term, node);
///  - the k-th partial sum of a root's left-fold Sum: (root label, k), id
///    "<label>+<k>".
///
/// A value's home block is the block holding its BuildStep; a value with no
/// BuildStep anywhere (escape-only) is homed at the innermost block whose
/// outputs list it; a leaf is homed at the root block (0). An edge is sliced
/// when the consumer's home block is nested deeper than the producer's (it
/// reads a slice of a value homed outside its own block). A value's
/// `persistent` is the engine's own flag, \c TableCell::persistent of its
/// production cell at its home block in the cell table \c build_cell_table
/// derives from the schedule (a Build cell or, for a leaf, its Leaf cell;
/// an escape-only value takes its Assemble cell there; a value with no cell is
/// never stored, so false).
struct NodePlace {
  std::string term;
  std::size_t node = 0;
  std::string id;
};

/// "<label>[<sliced modes>]" for a node: its tensor (or variable) label, or
/// its scalar for a constant, with the node_slice_mask() modes in brackets
template <typename Node>
std::string node_label(Node const& n) {
  std::string out;
  if (n->is_tensor())
    out = toUtf8(n->as_tensor().label());
  else if (n->is_variable())
    out = toUtf8(n->as_variable().label());
  else
    out = toUtf8(io::serialization::v1::to_string(n->expr()));
  auto const& mask = n->node_slice_mask();
  if (!mask.empty()) {
    out += "[";
    for (std::size_t k = 0; k < mask.size(); ++k) {
      if (k) out += ",";
      out += index_label(mask[k].first);
    }
    out += "]";
  }
  return out;
}

/// "<result>[...] = <left>[...] <op> <right>[...]", or the label alone on a
/// leaf
template <typename Node>
std::string node_expr(Node const& n) {
  if (n.leaf()) return node_label(n);
  auto const op = n->op_type();
  std::string const sym = op && *op == EvalOp::Sum ? " + " : " * ";
  return node_label(n) + " = " + node_label(n.left()) + sym +
         node_label(n.right());
}

}  // namespace detail

/// Value ids: smallest occurrence "<term key>#<ordinal>", "<term key>#scale",
/// "<root label>+<k>". `roots` supplies term keys in forest order
/// (roots[r].terms[t].key for root r, summand t).
///
/// \details Static points are matched to forest nodes by reproducing \c
/// compute_dag_boulevard's numbering: one post-order (left, right, node)
/// over every tree in forest order, leaves included, one counter across
/// trees. This needs nothing beyond the forest itself and resolves every
/// occurrence (each occurrence is one point), whereas going through
/// \c build_value_key_node_map / \c resolve_value_node would give one node
/// per value and lose which term each other occurrence sits in. The count is
/// cross-checked against \c RichSchedule::num_points.
///
/// \param is_volatile the node-level volatility predicate the evaluation
///        runs with (see \c evaluate_ordered_schedule), which decides
///        persistence; empty means no leaf is volatile
/// \throws sequant::Exception when the forest and \p roots disagree (tree or
///         summand count, a Sum spine node that is not an in-place Sum), on
///         an Adjoint node (not addressable by a pinned tree), or when a
///         value, block step or legality record cannot be resolved
template <meta::eval_node_range R>
DagRecord record_dag(
    RichSchedule const& rich, LegalitySchedule const& legality,
    OrderedSchedule const& ordered, R const& forest,
    std::vector<RootRecord> const& roots,
    std::function<bool(std::ranges::range_value_t<R> const&)> const&
        is_volatile = {}) {
  // nodes are keyed by address across two passes over the forest
  static_assert(
      std::is_lvalue_reference_v<std::ranges::range_reference_t<R const>>,
      "record_dag needs a range that yields node references");
  using Node = std::ranges::range_value_t<R>;
  using detail::NodePlace;
  using detail::spell_modes;

  // ---- 1. name every forest node by its position
  std::unordered_map<Node const*, NodePlace> place;

  // a term's subtree: an optional #scale top, then the pinned tree
  auto const classify_term = [&](Node const& top, std::string const& key) {
    Node const* tree = &top;
    bool const scaled = !top.leaf() && top->op_type() &&
                        *top->op_type() == EvalOp::Product &&
                        top.right().leaf() && top.right()->is_constant();
    if (scaled) tree = &top.left();
    std::size_t ordinal = 0;
    auto const walk = [&](auto const& self, Node const& n) -> void {
      if (!n.leaf()) {
        if (n->op_type() && *n->op_type() == EvalOp::Adjoint) {
          // the pinned (conjugated) factor: one ordinal, as a leaf; its
          // bare-label operand and Constant{1} sentinel are sub-positions
          std::string const id = key + "#" + std::to_string(ordinal);
          place[&n.left()] = {key, ordinal, id + "/0"};
          place[&n.right()] = {key, ordinal, id + "/1"};
        } else {
          self(self, n.left());
          self(self, n.right());
        }
      }
      place[&n] = {key, ordinal, key + "#" + std::to_string(ordinal)};
      ++ordinal;
    };
    walk(walk, *tree);
    if (scaled) {
      place[&top.right()] = {key, ordinal, key + "#scalar"};
      place[&top] = {key, ordinal + 1, key + "#scale"};
    }
  };

  std::size_t r = 0;
  for (auto const& tree : forest) {
    if (r == roots.size())
      throw Exception("record_dag: the forest has more trees than roots (" +
                      std::to_string(roots.size()) + ")");
    auto const& root = roots[r++];
    std::size_t const nterms = root.terms.size();
    if (nterms == 0)
      throw Exception("record_dag: root " + root.label + " has no terms");
    // a Sum root is a left fold ((s1 + s2) + s3) ...: walk its left spine
    std::vector<Node const*> summands(nterms);
    Node const* n = &tree;
    for (std::size_t k = nterms - 1; k >= 1; --k) {
      if (n->leaf() || !(*n)->op_type() || *(*n)->op_type() != EvalOp::Sum ||
          !(*n)->accumulate_in_place())
        throw Exception("record_dag: root " + root.label + " has " +
                        std::to_string(nterms) +
                        " terms but its tree is not a left-fold Sum of as "
                        "many summands");
      place[n] = {root.label, k, root.label + "+" + std::to_string(k)};
      summands[k] = &n->right();
      n = &n->left();
    }
    summands[0] = n;
    for (std::size_t t = 0; t < nterms; ++t)
      classify_term(*summands[t], root.terms[t].key);
  }
  if (r != roots.size())
    throw Exception("record_dag: " + std::to_string(roots.size()) +
                    " roots but the forest has " + std::to_string(r) +
                    " trees");

  // ---- 2. static points, as compute_dag_boulevard numbers them
  std::vector<Node const*> at_point;
  for (auto const& tree : forest) {
    // iterative post-order (a Sum spine is as deep as the term count)
    std::vector<std::pair<Node const*, bool>> stack{{&tree, false}};
    while (!stack.empty()) {
      auto [n, expanded] = stack.back();
      stack.pop_back();
      if (expanded || n->leaf()) {
        at_point.push_back(n);
        continue;
      }
      stack.push_back({n, true});
      stack.push_back({&n->right(), false});
      stack.push_back({&n->left(), false});
    }
  }
  if (at_point.size() != rich.num_points)
    throw Exception("record_dag: the forest has " +
                    std::to_string(at_point.size()) +
                    " nodes but the schedule has " +
                    std::to_string(rich.num_points) + " static points");
  auto const place_of = [&](std::size_t point) -> NodePlace const& {
    if (point >= at_point.size())
      throw Exception("record_dag: static point " + std::to_string(point) +
                      " is past the forest");
    auto const it = place.find(at_point[point]);
    if (it == place.end())
      throw Exception("record_dag: static point " + std::to_string(point) +
                      " is not a node of any root's terms");
    return it->second;
  };

  // ---- 3. value ids: the smallest occurrence id
  std::size_t const ncells = rich.cells.size();
  std::vector<std::string> id_of(ncells);
  for (std::size_t v = 0; v < ncells; ++v) {
    auto const& c = rich.cells[v];
    if (c.occurrences.empty())
      throw Exception("record_dag: value " + std::to_string(v) +
                      " has no occurrence to name it by");
    for (auto const& occ : c.occurrences) {
      auto const& id = place_of(occ.point).id;  // bounds-checks the point
      // identity cross-check: the node at this point is this value's node
      if (std::size_t const h = (*at_point[occ.point])->hash_value();
          h != c.hash)
        throw Exception(
            "record_dag: static point " + std::to_string(occ.point) +
            " holds node " + detail::hex_string(h) + " but value " +
            std::to_string(v) + " is node " + detail::hex_string(c.hash));
      if (id_of[v].empty() || id < id_of[v]) id_of[v] = id;
    }
  }
  auto const vid_name = [&](std::size_t v) -> std::string const& {
    if (v >= ncells)
      throw Exception("record_dag: value id " + std::to_string(v) +
                      " is not a value of the schedule");
    return id_of[v];
  };

  DagRecord dag;

  // ---- 4. blocks, preorder over ordered.root
  // value -> the block holding its BuildStep, and (for escape-only values)
  // the innermost block whose outputs list it: a child block's outputs are
  // seen before its ancestors' (they are read after its steps), so the first
  // escape site recorded is the innermost. Kept apart so that no outputs
  // entry ever overrides a build block.
  std::vector<std::optional<std::size_t>> build_block(ncells),
      escape_block(ncells);
  std::vector<std::size_t> nest_depth;  // per block id: number of ancestors
  std::vector<CellScope> scope_of;      // per block id: its cell-table scope
  auto const visit = [&](auto const& self, ScopeBlock const& b,
                         std::optional<std::size_t> parent) -> std::size_t {
    std::size_t const id = dag.blocks.size();
    dag.blocks.emplace_back();
    nest_depth.push_back(parent ? nest_depth[*parent] + 1 : 0);
    scope_of.push_back(parent ? scope_of[*parent] : CellScope{});
    if (parent)
      scope_of.back().path.push_back({b.level.key(), b.latitude_ordinal});
    {
      BlockRecord& rec = dag.blocks.back();
      rec.id = id;
      rec.parent = parent;
      if (parent) {
        rec.axis = mode_spelling(b.axis);
        rec.kind = detail::kind_name(b.kind);
      } else {
        rec.kind = "root";
      }
      rec.depth = b.level.depth;
      rec.slot = b.level.loop_slot;
      rec.altitude = b.level.altitude_ordinal;
      rec.pass = b.level.latitude_ordinal;
    }
    std::vector<StepRecord> steps;
    for (Step const& s : b.steps) {
      if (auto const* build = std::get_if<BuildStep>(&s.value)) {
        steps.push_back({vid_name(build->value_id), std::nullopt});
        build_block[build->value_id] = id;
      } else {
        steps.push_back(
            {std::nullopt, self(self, std::get<ScopeBlock>(s.value), id)});
      }
    }
    BlockRecord& rec = dag.blocks[id];  // the vector may have grown
    rec.steps = std::move(steps);
    for (auto const& [v, kind] : b.outputs) {
      rec.outputs.emplace_back(vid_name(v), detail::output_kind_name(kind));
      if (!escape_block[v]) escape_block[v] = id;
    }
    return id;
  };
  visit(visit, ordered.root, std::nullopt);
  // second pass: escape-only values take their innermost escape block
  std::vector<std::size_t> home_of(ncells, 0);
  for (std::size_t v = 0; v < ncells; ++v) {
    if (build_block[v])
      home_of[v] = *build_block[v];
    else if (escape_block[v])
      home_of[v] = *escape_block[v];
    else if (!rich.cells[v].is_leaf)
      throw Exception("record_dag: value " + id_of[v] +
                      " is not a leaf and no block builds or outputs it");
  }

  // ---- 5. persistence: the engine's flag, read off the cell table the
  // ordered executor builds from this schedule. The flag needs no lives, so
  // the batch count is a constant.
  auto const resolve = [&](std::size_t vid) -> Node const& {
    return *at_point[rich.cells.at(vid).occurrences.front().point];
  };
  SlicedModeAssignment const sliced_modes =
      compute_sliced_mode_assignment(ordered, rich);
  CellTableInputs table_inputs = sequant::eval::detail::make_cell_table_inputs(
      ordered, rich, sliced_modes, resolve, is_volatile);
  table_inputs.n_batches_of = [](LoopKey const&) { return std::size_t{1}; };
  CellTable const table = build_cell_table(table_inputs);
  std::vector<std::optional<bool>> produced(ncells), assembled(ncells);
  for (TableCell const& cell : table.cells) {
    std::size_t const v = cell.value_id;
    if (v >= ncells)
      throw Exception("record_dag: the cell table names value " +
                      std::to_string(v) + " past the schedule");
    if (!(cell.scope == scope_of[home_of[v]])) continue;
    (cell.production.kind == ProductionKind::Assemble ? assembled
                                                      : produced)[v] =
        cell.persistent;
  }

  // ---- 6. values
  std::unordered_map<std::size_t, CellLegality const*> legality_of;
  for (auto const& cl : legality.cells) legality_of[cl.hash] = &cl;
  std::map<std::size_t, std::vector<std::size_t>> by_hash;  // node id
  for (std::size_t v = 0; v < ncells; ++v)
    by_hash[rich.cells[v].hash].push_back(v);

  for (std::size_t v = 0; v < ncells; ++v) {
    auto const& c = rich.cells[v];
    ValueRecord rec;
    rec.id = id_of[v];
    rec.expr = detail::node_expr(*at_point.at(c.occurrences.front().point));
    rec.hash = detail::hex_string(c.hash);
    rec.key = detail::hex_string(value_key_of(c));
    rec.leaf = c.is_leaf;
    rec.result_modes = spell_modes(c.carried);
    rec.home_modes = spell_modes(c.home_modes);
    rec.enclosing_modes = spell_modes(c.enclosing_modes);
    rec.divergent_modes = spell_modes(c.divergent_modes);
    rec.home_block = home_of[v];
    rec.home_depth = c.home_depth;
    rec.persistent = produced[v] ? *produced[v] : assembled[v].value_or(false);
    rec.first_use = c.first_use;
    rec.last_use = c.last_use;
    auto const lit = legality_of.find(value_key_of(c));
    if (lit == legality_of.end())
      throw Exception("record_dag: value " + rec.id +
                      " has no legality record");
    auto const& cl = *lit->second;
    rec.build_site = spell_modes(cl.build_site);
    rec.home_floor = spell_modes(cl.home_floor);
    rec.forced_split_axes = spell_modes(cl.forced_split_axes);
    for (auto const& ac : cl.per_axis)
      rec.per_axis.push_back(
          {mode_spelling(ac.axis), detail::role_name(ac.role)});
    for (auto const& occ : c.occurrences) {
      auto const& p = place_of(occ.point);
      OccurrenceRecord o;
      o.term = p.term;
      o.node = p.node;
      // the consumer's position when it is in the same term (or, for a
      // partial sum, the same root's spine); none at a root or across
      if (occ.consumer_point != occ.point) {
        auto const& cp = place_of(occ.consumer_point);
        if (cp.term == p.term) o.consumer = cp.node;
      }
      o.home = spell_modes(occ.home);
      for (auto const& e : occ.ectx)
        o.enclosing.push_back(mode_spelling(e.first));
      for (auto const& [ix, kind] : occ.opens)
        o.opens.push_back({mode_spelling(ix), kind == BatchModeType::External
                                                  ? ModeKind::External
                                                  : ModeKind::Contracted});
      rec.occurrences.push_back(std::move(o));
    }
    for (std::size_t s : by_hash.at(c.hash))
      if (s != v) rec.siblings.push_back(id_of[s]);
    std::ranges::sort(rec.siblings);
    dag.values.push_back(std::move(rec));
  }

  // ---- 7. loops and loop order
  auto const loop_rec = [&](std::pair<std::wstring, int> const& k) {
    if (k.second < 0)
      throw Exception("record_dag: loop " + toUtf8(k.first) + " has no slot (" +
                      std::to_string(k.second) + ")");
    auto const it = rich.loop_kind.find(k);
    return LoopRecord{toUtf8(k.first), static_cast<std::size_t>(k.second),
                      it == rich.loop_kind.end()
                          ? std::string{}
                          : detail::kind_name(it->second)};
  };
  for (auto const& [k, kind] : rich.loop_kind) dag.loops.push_back(loop_rec(k));
  for (auto const& [pair, witness] : rich.loop_order)
    dag.loop_order.emplace_back(loop_rec(pair.first), loop_rec(pair.second),
                                witness);

  // ---- 8. operand edges; sliced = the consumer's home block is nested
  // deeper than the producer's (it reads a slice of a value homed outside its
  // own block). Nesting depth, not DagScopeLevel::depth (a loop-group
  // identity: a root child block can carry depth 2).
  std::set<std::pair<std::size_t, std::size_t>> seen;
  for (std::size_t to = 0; to < ncells; ++to) {
    auto const it = ordered.operand_vids.find(to);
    if (it == ordered.operand_vids.end()) continue;
    for (std::size_t from : it->second) {
      if (!seen.insert({from, to}).second) continue;
      dag.edges.push_back(
          {vid_name(from), id_of[to],
           nest_depth[home_of[to]] > nest_depth[home_of.at(from)]});
    }
  }
  return dag;
}

}  // namespace sequant::eval::schedule

#endif  // SEQUANT_EVAL_SCHEDULE_DAG_RECORD_HPP
