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

inline std::string role_name(OccurrenceRole r) {
  switch (r) {
    case OccurrenceRole::Lockstep:
      return "Lockstep";
    case OccurrenceRole::Assembled:
      return "Assembled";
  }
  throw Exception("record_dag: unknown OccurrenceRole");
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
///  - a unary node (an Adjoint, a conjugated factor, or a Power, a
///    TensorPower factor) takes one tree ordinal as a leaf would; its operand
///    is id "<key>#<ordinal>/0" and its Constant sentinel (Adjoint's
///    Constant{1}, Power's exponent) "<key>#<ordinal>/1", both at the same
///    (term, node); when the operand is a subtree (an inlined equation), its
///    other nodes are "<operand id>.<j>" in left-first post-order, a unary
///    node among them taking one j with its own "/0" and "/1" below, all at
///    the same (term, node);
///  - the k-th partial sum of a root's left-fold Sum: (root label, k), id
///    "<label>+<k>".
///
/// A value's home block is the block holding its BuildStep; a value with no
/// BuildStep anywhere (escape-only) is homed at the innermost block whose
/// outputs list it; a leaf is homed at the root block (0). An edge is sliced
/// when the consumer's home block is nested deeper than the producer's (it
/// reads a slice of a value homed outside its own block). A value's
/// `persistent` is the engine's own flag, \c TableCell::persistent of the
/// value's outermost cell in the cell table \c build_cell_table derives from
/// the schedule: the complete form the rest of the DAG consumes, which is the
/// Assemble cell at the parent of the loop block that produces it per batch,
/// else its Build or Leaf cell. A leaf no cell reads is never stored, so
/// false.
struct NodePlace {
  std::string term;
  std::size_t node = 0;
  std::string id;
};

/// a node's bare label: its tensor (or variable) label, or its scalar for a
/// constant
template <typename Node>
std::string bare_label(Node const& n) {
  if (n->is_tensor()) return toUtf8(n->as_tensor().label());
  if (n->is_variable()) return toUtf8(n->as_variable().label());
  return toUtf8(io::serialization::v1::to_string(n->expr()));
}

/// UTF-8 "\u03a3", the partial-sum marker of a reduced loop axis in `expr`
inline constexpr char const* kPartialMarker = "\xCE\xA3";

/// "<label>[<form>]": a value in one of its forms: \p slices, the carried
/// modes a read or its own home loops slice, and \p partials, the loop axes
/// its node reduces when the form is a per-batch partial sum, each marked
/// with \ref kPartialMarker; bare when neither applies, which is the complete
/// value. Mode spellings ("label:basis") are shown by label.
inline std::string form_label(std::string const& label,
                              std::vector<std::string> const& slices,
                              std::vector<std::string> const& partials) {
  std::string brackets;
  auto const add = [&](std::string const& mode) {
    brackets += (brackets.empty() ? "" : ",") + mode.substr(0, mode.find(':'));
  };
  for (auto const& m : slices) add(m);
  for (auto const& m : partials) add(kPartialMarker + m);
  return brackets.empty() ? label : label + "[" + brackets + "]";
}

}  // namespace detail

/// Value ids: the value's first occurrence in forest order, "<term key>#<
/// ordinal>", "<term key>#scale" or "<root label>+<k>". That occurrence is
/// also the engine's frame for the value's cell (\c home_modes_in_cell_frame:
/// "the rich cell's own frame -- the first occurrence's labels"), so the id's
/// term is the term whose index labels `expr`, the mode lists and the legality
/// fields of the record are spelled in; every other occurrence has its own
/// labels, recorded per occurrence. `roots` supplies term keys in forest order
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
///         summand count, a Sum spine node that is not an in-place Sum), or
///         when a value, block step or legality record cannot be resolved
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
    auto const is_unary_node = [](Node const& n) {
      return !n.leaf() && n->op_type() && is_unary(*n->op_type());
    };
    // a sub-position `id` at the ordinal and what hangs below it: a unary
    // node's operand and sentinel at id/0 and id/1, the other nodes of a
    // subtree (an inlined operand) at id.<j> in left-first post-order, a
    // unary node among them counting as one
    auto const place_sub = [&](auto const& self, Node const& n,
                               std::string const& id) -> void {
      place[&n] = {key, ordinal, id};
      if (n.leaf()) return;
      if (is_unary_node(n)) {
        self(self, n.left(), id + "/0");
        self(self, n.right(), id + "/1");
        return;
      }
      std::size_t j = 0;
      auto const below = [&](auto const& below_self, Node const& m) -> void {
        if (!m.leaf() && !is_unary_node(m)) {
          below_self(below_self, m.left());
          below_self(below_self, m.right());
          place[&m] = {key, ordinal, id + "." + std::to_string(j++)};
        } else {
          self(self, m, id + "." + std::to_string(j++));
        }
      };
      below(below, n.left());
      below(below, n.right());
    };
    auto const walk = [&](auto const& self, Node const& n) -> void {
      if (!n.leaf()) {
        if (is_unary_node(n)) {
          // the pinned factor (conjugated or raised): one ordinal, as a leaf;
          // its operand and Constant sentinel are sub-positions
          std::string const id = key + "#" + std::to_string(ordinal);
          place_sub(place_sub, n.left(), id + "/0");
          place_sub(place_sub, n.right(), id + "/1");
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

  // ---- 3. value ids: the first occurrence, whose labels are the cell's frame
  std::size_t const ncells = rich.cells.size();
  std::vector<std::string> id_of(ncells);
  // the short name of a term's values: the root label, with the term index
  // when the root has several terms
  std::unordered_map<std::string, std::string> short_of_term;
  for (auto const& root : roots)
    for (std::size_t t = 0; t < root.terms.size(); ++t)
      short_of_term.emplace(root.terms[t].key,
                            root.terms.size() > 1
                                ? root.label + "/" + std::to_string(t)
                                : root.label);
  std::vector<std::string> name_of(ncells);
  std::unordered_map<std::string, std::size_t> value_of_occurrence;  // by id
  for (std::size_t v = 0; v < ncells; ++v) {
    auto const& c = rich.cells[v];
    if (c.occurrences.empty())
      throw Exception("record_dag: value " + std::to_string(v) +
                      " has no occurrence to name it by");
    for (auto const& occ : c.occurrences) {
      value_of_occurrence.emplace(place_of(occ.point).id, v);  // bounds-checks
      // identity cross-check: the node at this point is this value's node
      if (std::size_t const h = (*at_point[occ.point])->hash_value();
          h != c.hash)
        throw Exception(
            "record_dag: static point " + std::to_string(occ.point) +
            " holds node " + detail::hex_string(h) + " but value " +
            std::to_string(v) + " is node " + detail::hex_string(c.hash));
    }
    auto const& first = place_of(c.occurrences.front().point);
    id_of[v] = first.id;
    auto const sit = short_of_term.find(first.term);
    name_of[v] = sit == short_of_term.end()
                     ? first.id
                     : sit->second + first.id.substr(first.term.size());
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
  auto const visit = [&](auto const& self, ScopeBlock const& b,
                         std::optional<std::size_t> parent) -> std::size_t {
    std::size_t const id = dag.blocks.size();
    dag.blocks.emplace_back();
    nest_depth.push_back(parent ? nest_depth[*parent] + 1 : 0);
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
  // ordered executor builds from this schedule, on the value's outermost
  // cell (the fewest enclosing loop instances; an Assemble before a Build at
  // the same depth): the complete form the rest of the DAG consumes. The flag
  // needs no lives, so the batch count is a constant.
  auto const resolve = [&](std::size_t vid) -> Node const& {
    return *at_point[rich.cells.at(vid).occurrences.front().point];
  };
  SlicedModeAssignment const sliced_modes =
      compute_sliced_mode_assignment(ordered, rich);
  CellTableInputs table_inputs = sequant::eval::detail::make_cell_table_inputs(
      ordered, rich, sliced_modes, resolve, is_volatile);
  table_inputs.n_batches_of = [](LoopKey const&) { return std::size_t{1}; };
  CellTable const table = build_cell_table(table_inputs);
  struct Outermost {
    std::size_t depth = 0;
    bool assemble = false;
    bool persistent = false;
    bool any = false;
  };
  std::vector<Outermost> outermost(ncells);
  // the value's production form: its deepest Build (or Leaf) cell
  std::vector<TableCell const*> production(ncells, nullptr);
  for (TableCell const& cell : table.cells) {
    std::size_t const v = cell.value_id;
    if (v >= ncells)
      throw Exception("record_dag: the cell table names value " +
                      std::to_string(v) + " past the schedule");
    std::size_t const depth = cell.scope.path.size();
    bool const assemble = cell.production.kind == ProductionKind::Assemble;
    if (!assemble &&
        (!production[v] || production[v]->scope.path.size() < depth))
      production[v] = &cell;
    auto& best = outermost[v];
    if (best.any && (best.depth < depth ||
                     (best.depth == depth && best.assemble && !assemble)))
      continue;
    best = {depth, assemble, cell.persistent, true};
  }
  for (std::size_t v = 0; v < ncells; ++v)
    if (!outermost[v].any && !rich.cells[v].is_leaf)
      throw Exception("record_dag: value " + id_of[v] +
                      " has no cell in the cell table");

  // ---- 6. values
  std::unordered_map<std::size_t, CellLegality const*> legality_of;
  for (auto const& cl : legality.cells) legality_of[cl.hash] = &cl;
  std::map<std::size_t, std::vector<std::size_t>> by_hash;  // node id
  for (std::size_t v = 0; v < ncells; ++v)
    by_hash[rich.cells[v].hash].push_back(v);

  // Forms (detail::form_label), spelled in the cell frame (the carried list,
  // so the labels agree with result_modes). A value's own form is its
  // production: the carried axes its production cell is sliced on and, for
  // every slice-mask axis the value does not carry (the node reduces it),
  // the per-batch partial over it. An operand's form is how the value reads
  // it, from the cell table: the source cell's own slicing plus the read's,
  // and the partials when the source is the operand's per-batch Build cell
  // rather than its assembled form.
  auto const node_of = [&](std::size_t w) -> Node const& {
    return *at_point[rich.cells[w].occurrences.front().point];
  };
  auto const reduced_of = [&](std::size_t w) {
    auto const carried = spell_modes(rich.cells[w].carried);
    std::vector<std::string> out;
    for (auto const& [ix, kind] : node_of(w)->node_slice_mask())
      if (auto const m = mode_spelling(ix);
          std::ranges::find(carried, m) == carried.end())
        out.push_back(m);
    return out;
  };
  auto const sliced_labels = [&](std::size_t w, auto const& positions,
                                 std::vector<std::string>& out) {
    auto const carried = spell_modes(rich.cells[w].carried);
    for (auto const& [pos, key] : positions)
      if (pos < carried.size() &&
          std::ranges::find(out, carried[pos]) == out.end())
        out.push_back(carried[pos]);
  };
  // as an operand: a leaf by its tensor label, an intermediate by its name
  auto const label_of = [&](std::size_t w) {
    return rich.cells[w].is_leaf ? detail::bare_label(node_of(w)) : name_of[w];
  };
  // the left-hand side of expr: the name in the form the value is produced in
  auto const own_form = [&](std::size_t w) {
    std::vector<std::string> slices;
    if (production[w]) sliced_labels(w, production[w]->sliced, slices);
    return detail::form_label(name_of[w], slices, reduced_of(w));
  };
  // (consumer value, operand value) -> the read's source cell and slicing
  struct ReadForm {
    CellId source;
    container::svector<std::pair<std::size_t, LoopKey>> const* slice;
  };
  std::map<std::pair<std::size_t, std::size_t>, ReadForm> read_form;
  for (Read const& r : table.reads) {
    if (r.consumer >= table.cells.size() || r.source >= table.cells.size())
      continue;
    read_form.try_emplace(
        {table.cells[r.consumer].value_id, r.operand_value_id},
        ReadForm{r.source, &r.slice});
  }
  // an operand of value v's node: the value at that child, as v reads it
  auto const operand_of = [&](std::size_t v, Node const& child) {
    auto const pit = place.find(&child);
    auto const wit = pit == place.end()
                         ? value_of_occurrence.end()
                         : value_of_occurrence.find(pit->second.id);
    if (wit == value_of_occurrence.end()) return detail::bare_label(child);
    std::size_t const w = wit->second;
    auto const rit = read_form.find({v, w});
    if (rit == read_form.end()) return label_of(w);
    TableCell const& source = table.cells[rit->second.source];
    std::vector<std::string> slices;
    sliced_labels(w, source.sliced, slices);
    sliced_labels(w, *rit->second.slice, slices);
    return detail::form_label(label_of(w), slices,
                              source.production.kind == ProductionKind::Build
                                  ? reduced_of(w)
                                  : std::vector<std::string>{});
  };

  for (std::size_t v = 0; v < ncells; ++v) {
    auto const& c = rich.cells[v];
    ValueRecord rec;
    rec.id = id_of[v];
    rec.name = name_of[v];
    {
      auto const& node = node_of(v);
      rec.expr = own_form(v) + " = ";
      if (node.leaf()) {
        rec.expr += detail::bare_label(node);
      } else {
        auto const op = node->op_type();
        rec.expr += operand_of(v, node.left());
        if (op == EvalOp::Adjoint)
          rec.expr += "⁺";
        else
          rec.expr += (op == EvalOp::Sum     ? " + "
                       : op == EvalOp::Power ? " ^ "
                                             : " * ") +
                      operand_of(v, node.right());
      }
    }
    rec.hash = detail::hex_string(c.hash);
    rec.key = detail::hex_string(value_key_of(c));
    rec.leaf = c.is_leaf;
    rec.result_modes = spell_modes(c.carried);
    rec.home_modes = spell_modes(c.home_modes);
    rec.enclosing_modes = spell_modes(c.enclosing_modes);
    rec.divergent_modes = spell_modes(c.divergent_modes);
    rec.home_block = home_of[v];
    rec.home_depth = c.home_depth;
    rec.persistent = outermost[v].persistent;
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
      // a per-occurrence role is on a carried position of the value; spell
      // the position in this occurrence's own frame, as home and enclosing
      // are (ac.axis is in the first occurrence's)
      for (auto const& ac : cl.per_axis)
        for (auto const& por : ac.per_occurrence) {
          if (por.point != occ.point) continue;
          auto const pit =
              std::find(c.carried.begin(), c.carried.end(), ac.axis);
          std::size_t const pos =
              static_cast<std::size_t>(pit - c.carried.begin());
          if (pit == c.carried.end() || pos >= occ.carried.size())
            throw Exception("record_dag: value " + rec.id +
                            " has an occurrence role on " +
                            mode_spelling(ac.axis) +
                            ", which is not a carried position of the value");
          o.roles.push_back(
              {mode_spelling(occ.carried[pos]), detail::role_name(por.role)});
        }
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
