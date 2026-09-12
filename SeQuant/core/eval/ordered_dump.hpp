#ifndef SEQUANT_CORE_EVAL_ORDERED_DUMP_HPP
#define SEQUANT_CORE_EVAL_ORDERED_DUMP_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

///
/// \file ordered_dump.hpp
/// \brief The ordered (DAG) evaluator's environment-gated diagnostic dumps,
/// in ONE place.
///
/// \details Every dump below is OFF unless its environment variable is set;
/// when it is unset the call site is a single \c std::getenv test and the
/// evaluated result is byte-identical. Nothing here participates in
/// evaluation -- these functions only print.
///
/// The dumps take their operands as template parameters rather than including
/// the schedule/table/executor headers, so this header sits BELOW all of them
/// and can be included from any of them without an include cycle.
///
/// ENVIRONMENT VARIABLES (the complete list of what SeQuant reads)
///
/// Ordered-evaluator dumps (this header):
///
/// | Variable                  | Prints                                      |
/// |---------------------------|---------------------------------------------|
/// | SEQUANT_DUMP_SCHEDULE     | `[levels]`, `[sched]`, `[sched-collapse]`,  |
/// |                           | `[sched-materialize]`, `[sched-nest]` and   |
/// |                           | `[sched-tree]` -- how the ordered schedule  |
/// |                           | was built and the block tree it became      |
/// |                           | (plus `build_ordered_schedule`'s own        |
/// |                           | `dump_reject` rejection report, which reads |
/// |                           | that builder's local frame and so stays     |
/// |                           | there).                                     |
/// | SEQUANT_DUMP_USEINDUCED   | `[useinduced]` -- each use-induced slicing  |
/// |                           | fact (value, mode, position, loop,          |
/// |                           | consumer) recorded on an operand fetch.     |
/// | SEQUANT_DUMP_OPENS        | `[opens]` -- the batch loops each node      |
/// |                           | opens, with its carried modes.              |
/// | SEQUANT_DUMP_CELLS_OF     | `[cells]` / `[occs]` for the listed value   |
/// |                           | ids (`=<vid>[,<vid>...]`): every cell of    |
/// |                           | each value, and every occurrence layout and |
/// |                           | table read of it.                           |
/// | SEQUANT_DUMP_ROOT_NORMS   | `[build-norm]` / `[root-norm]` -- the norm  |
/// |                           | of every root-scope build and of every      |
/// |                           | forest root at the combine, for comparing   |
/// |                           | term values across schedules.               |
/// | SEQUANT_UT_READ_DIAG      | `[READ]` -- every table-driven operand      |
/// |                           | read, with consumer, source and slices.     |
/// | SEQUANT_UT_BLOCK_DIAG     | `[BLOCK]` / `[cell-registry]` -- per-batch  |
/// |                           | build and assemble steps, block skips, and  |
/// |                           | the end-of-call registry residency.         |
/// | SEQUANT_SCHED_DUMP        | `ORDERED_RUN_BLOCK` / `SCHEDULE_RUN_GROUP`  |
/// |                           | JSON lines for the schedule visualizer      |
/// |                           | (emitted by the executors themselves).      |
///
/// Cache / runtime knobs documented with their own facility (cache_manager.hpp,
/// runtime.hpp, optimize/): SEQUANT_NUM_THREADS, SEQUANT_FACTORIZER_DEBUG,
/// SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING, SEQUANT_UT_ACCESS_CLOCK,
/// SEQUANT_UT_LOOKUP_METER, SEQUANT_UT_PHASE, SEQUANT_UT_EVALIMPL,
/// SEQUANT_UT_BUILD_METER, SEQUANT_UT_BUILD_DUMP, SEQUANT_UT_DEFUSE_METER,
/// SEQUANT_UT_DEFUSE_DUMP, SEQUANT_UT_PEAK_COMPOSE, SEQUANT_REBUILD_TRACE,
/// SEQUANT_UT_STRICT_FILL_ONCE, SEQUANT_SYNC_STATS.
///
/// SEQUANT_UT_FORCE_SYNC and SEQUANT_UT_PROD_TR name per-op diagnostics in \c
/// apply_one_op_traced / \c Result::fence that the LIBRARY never reads itself:
/// the gate is evaluated at the CALL SITE, by the consuming backend (mpqc's
/// dry-run and wet drivers). They are listed here because the library's own
/// doxygen names them (\c eval.hpp, \c result.hpp), not because SeQuant reads
/// them.
///
/// Everything else prefixed SEQUANT_UT_ is a unit-test knob read by the tests
/// themselves, not by library code.
///

namespace sequant::eval::detail {

/// \brief Is dump \p var enabled? (Any value, including the empty string,
/// enables it; only an unset variable disables it.)
[[nodiscard]] inline bool dump_enabled(char const* var) noexcept {
  return std::getenv(var) != nullptr;
}

/// \brief The raw value of dump variable \p var, or nullptr when unset.
[[nodiscard]] inline char const* dump_value(char const* var) noexcept {
  return std::getenv(var);
}

// ---------------------------------------------------------------------------
// SEQUANT_DUMP_SCHEDULE
// ---------------------------------------------------------------------------

/// \brief `[levels]`: value \p v's pass was bumped by its operand \p o, which
/// is either carried across a loop or read inside the loop that reduces it.
inline void dump_level_bump(std::size_t v, std::size_t o, bool carried,
                            int level) {
  std::wcerr << L"[levels] v" << v << L" bumped by operand v" << o
             << (carried ? L" (carried)" : L" (reduction read inside)")
             << L" -> level " << level << L"\n";
}

/// \brief `[sched]`: the pairwise loop order the rich schedule witnesses, and
/// the per-instance loop chain (depth -> space#slot) realized from it.
template <typename LoopOrder, typename Types, typename Slots>
void dump_loop_chain(LoopOrder const& loop_order, Types const& types,
                     Slots const& type_slot) {
  for (auto const& [pair, witness] : loop_order)
    std::wcerr << L"[sched] loop_order " << pair.first.first << L"#"
               << pair.first.second << L" > " << pair.second.first << L"#"
               << pair.second.second << L" by v" << witness << L"\n";
  std::wcerr << L"[sched] per-instance loop chain: ";
  for (std::size_t d = 0; d < types.size(); ++d)
    std::wcerr << L"d" << d << L"=" << types[d].space().base_key() << L"#slot"
               << type_slot[d] << L" ";
  std::wcerr << L"\n";
}

/// \brief `[sched-collapse]`: a value with more than one non-local mode of ONE
/// space has its distinct per-instance escapes collapsed to fewer escapes (one
/// per depth == one per space). \p nonlocal is that value's non-LoopLocal axes
/// with a one-character role tag each, \p escapes its (depth, kind) list.
template <typename Nonlocal, typename Escapes>
void dump_sched_collapse(std::size_t hash, Nonlocal const& nonlocal,
                         Escapes const& escapes) {
  std::wcerr << L"[sched-collapse] hash=" << hash << L" nonlocal={";
  for (auto const& [ix, role] : nonlocal)
    std::wcerr << ix.full_label() << L":" << role << L" ";
  std::wcerr << L"} -> " << nonlocal.size() << L" modes COLLAPSE to "
             << escapes.size() << L" escapes(depth:kind)={";
  for (auto const& [d, k] : escapes)
    std::wcerr << d << L":"
               << (k == std::decay_t<decltype(k)>::AccumulateScatter
                       ? L"Scatter"
                       : L"Sum")
               << L" ";
  std::wcerr << L"}\n";
}

/// \brief `[sched-materialize]`: the members built AND escaped across a forced
/// loop split (their in-nest readers take the per-batch cell, the other pass
/// the assembled form).
template <typename Ids>
void dump_sched_materialize(Ids const& ids) {
  std::wcerr << L"[sched-materialize] " << ids.size()
             << L" member(s) built AND escaped across the forced split:";
  for (std::size_t v : ids) std::wcerr << L" v" << v;
  std::wcerr << L"\n";
}

/// \brief `[sched-nest]`: per nest, its outermost depth and the passes emitted
/// there, each with the number of builds it carries. \p count(cluster, pass)
/// returns that build count.
template <typename NestPasses, typename ClusterMin, typename CountFn>
void dump_sched_nests(NestPasses const& nest_passes,
                      ClusterMin const& cluster_min, CountFn&& count) {
  for (auto const& [c, ps] : nest_passes) {
    std::wcerr << L"[sched-nest] outermost depth " << cluster_min.at(c)
               << L" passes={";
    for (int k : ps) std::wcerr << k << L":" << count(c, k) << L" ";
    std::wcerr << L"}\n";
  }
}

/// \brief `[sched-tree]`: the block tree of a built schedule -- per block its
/// depth, axis, kind (contracted/external), own BuildStep count, child block
/// count and escaped outputs (value:sum|scatter) -- so a production run shows
/// which batch loops were realized and of which kind.
template <typename Block>
void dump_schedule_tree(Block const& b, int depth) {
  std::size_t builds = 0, children = 0;
  for (auto const& st : b.steps)
    std::visit(
        [&](auto const& alt) {
          if constexpr (requires { alt.steps; })
            ++children;
          else
            ++builds;
        },
        st.value);
  std::wcerr << L"[sched-tree] " << std::wstring(2 * depth, L' ') << L"depth="
             << depth << L" axis=" << (b.axis ? b.axis.full_label() : L"<root>")
             << L" kind="
             << (b.kind == BatchModeType::External ? L"external"
                                                   : L"contracted")
             << L" builds=" << builds << L" children=" << children << L" outs=";
  for (auto const& [ovid, okind] : b.outputs)
    std::wcerr << ovid << L":"
               << (okind == std::decay_t<decltype(okind)>::AccumulateSum
                       ? L"sum"
                   : okind == std::decay_t<decltype(okind)>::AccumulateScatter
                       ? L"scatter"
                       : L"other")
               << L" ";
  std::wcerr << L"\n";
  for (auto const& st : b.steps)
    std::visit(
        [&](auto const& alt) {
          if constexpr (requires { alt.steps; })
            dump_schedule_tree(alt, depth + 1);
        },
        st.value);
}

// ---------------------------------------------------------------------------
// SEQUANT_DUMP_USEINDUCED
// ---------------------------------------------------------------------------

/// \brief `[useinduced]`: value \p value_h is sliced at position \p pos on
/// loop \p loop because its consumer \p consumer_h is sliced on mode \p mode
/// there (Layer 2 use-induced slicing).
inline void dump_use_induced(std::size_t value_h, Index const& mode,
                             std::size_t pos, std::size_t loop,
                             std::size_t consumer_h) {
  std::cerr << "[useinduced] value_h=" << value_h
            << " M=" << toUtf8(mode.full_label()) << " pos=" << pos
            << " loop=" << loop << " consumer_h=" << consumer_h << "\n";
}

// ---------------------------------------------------------------------------
// SEQUANT_DUMP_OPENS
// ---------------------------------------------------------------------------

/// \brief `[opens]`: the batch loops node \p n opens at its own node (the
/// group nest structure the factorizer emits) with its carried modes -- which
/// same-space modes open at ONE node (a multi-loop group) vs at different
/// nodes (separate groups).
template <typename Node>
void dump_opens(Node const& n) {
  std::wcerr << L"[opens] hash=" << n->hash_value() << L" opened_here={";
  for (auto const& [ix, kind] : n->batch_loops_opened_here())
    std::wcerr << ix.full_label() << L":" << ix.space().base_key() << L":"
               << (kind == BatchModeType::Contracted ? L"C" : L"E") << L" ";
  std::wcerr << L"} carried={";
  for (auto const& c : n->canon_indices()) std::wcerr << c.full_label() << L" ";
  std::wcerr << L"}\n";
}

// ---------------------------------------------------------------------------
// SEQUANT_DUMP_CELLS_OF
// ---------------------------------------------------------------------------

/// \brief Parse SEQUANT_DUMP_CELLS_OF's `<vid>[,<vid>...]` value.
[[nodiscard]] inline std::set<std::size_t> dump_cells_of_values(
    char const* spec) {
  std::set<std::size_t> want;
  std::istringstream toks{spec};
  for (std::string tok; std::getline(toks, tok, ',');)
    if (!tok.empty()) want.insert(std::stoul(tok));
  return want;
}

/// \brief `[cells]`: every cell of the wanted values -- kind, scope path,
/// sliced positions, partial_over, scatter map, carried modes, and the
/// produce-if-absent / persistent flags.
template <typename Table, typename Rich>
void dump_cells_of(Table const& table, Rich const& rich,
                   std::set<std::size_t> const& want) {
  for (std::size_t cid = 0; cid < table.cells.size(); ++cid) {
    auto const& c = table.cells[cid];
    if (!want.count(c.value_id)) continue;
    using Kind = std::decay_t<decltype(c.production.kind)>;
    std::cerr << "[cells] value " << c.value_id << " cell#" << cid << " "
              << (c.production.kind == Kind::Build      ? "Build"
                  : c.production.kind == Kind::Assemble ? "Assemble"
                                                        : "Leaf")
              << " path={";
    for (auto const& [k, lat] : c.scope.path)
      std::cerr << "d" << k.depth << "#" << k.loop_slot << " ";
    std::cerr << "} sliced={";
    for (auto const& [p, k] : c.sliced)
      std::cerr << "pos" << p << "@d" << k.depth << "#" << k.loop_slot << " ";
    std::cerr << "} partial_over={";
    for (auto const& k : c.partial_over)
      std::cerr << "d" << k.depth << "#" << k.loop_slot << " ";
    std::cerr << "} scatter={";
    for (auto const& [p, k] : c.production.scatter_map)
      std::cerr << "pos" << p << "@d" << k.depth << "#" << k.loop_slot << " ";
    std::cerr << "} carried={";
    for (auto const& ix : rich.cells[c.value_id].carried)
      std::cerr << toUtf8(std::wstring(ix.full_label())) << ' ';
    std::cerr << "} pia=" << c.produce_if_absent
              << " persistent=" << c.persistent << "\n";
  }
}

/// \brief `[occs]`: per wanted value, every occurrence's OWN index order (its
/// array layout) with its consumer, and every table read of the value with its
/// slices -- the two must agree position-for-position on ONE layout.
template <typename Table, typename Rich>
void dump_occs_of(Table const& table, Rich const& rich,
                  std::set<std::size_t> const& want) {
  std::map<std::size_t, std::pair<std::size_t, std::size_t>> point_to;
  for (std::size_t vid = 0; vid < rich.cells.size(); ++vid)
    for (std::size_t o = 0; o < rich.cells[vid].occurrences.size(); ++o)
      point_to[rich.cells[vid].occurrences[o].point] = {vid, o};
  auto const labels = [](container::svector<Index> const& c) {
    std::string s;
    for (auto const& ix : c) s += toUtf8(std::wstring(ix.full_label())) + " ";
    return s;
  };
  for (std::size_t vid : want) {
    if (vid >= rich.cells.size()) continue;
    auto const& vc = rich.cells[vid];
    std::cerr << "[occs] value " << vid << " cell carried={"
              << labels(vc.carried) << "}\n";
    for (auto const& occ : vc.occurrences) {
      std::cerr << "  occurrence carried={" << labels(occ.carried)
                << "} loop_slot={";
      for (int sl : occ.loop_slot) std::cerr << sl << " ";
      std::cerr << "} consumer=";
      if (auto it = point_to.find(occ.consumer_point);
          it != point_to.end() && occ.consumer_point != occ.point) {
        auto const& co =
            rich.cells[it->second.first].occurrences[it->second.second];
        std::cerr << it->second.first << " carried={" << labels(co.carried)
                  << "}";
      } else {
        std::cerr << "(root)";
      }
      std::cerr << "\n";
    }
    for (auto const& r : table.reads) {
      if (table.cells[r.source].value_id != vid) continue;
      std::cerr << "  read by consumer value "
                << table.cells[r.consumer].value_id << " cell#" << r.consumer
                << " slices={";
      for (auto const& [p, key] : r.slice)
        std::cerr << "pos" << p << "@d" << key.depth << "#" << key.loop_slot
                  << " ";
      std::cerr << "}\n";
    }
  }
}

// ---------------------------------------------------------------------------
// SEQUANT_UT_READ_DIAG
// ---------------------------------------------------------------------------

/// \brief `[READ]`: one table-driven operand read -- its consumer cell, the
/// value, the source cell and the declared slices (each resolved against the
/// in-scope batch context where it names an enclosing loop).
template <typename Cell, typename Read, typename Ctx, typename SameKey>
void dump_read(std::size_t consumer, Cell const& ccell, std::size_t vid,
               Read const& r, Cell const& src_cell, Ctx const& ctx,
               SameKey&& same_key) {
  using Kind = std::decay_t<decltype(ccell.production.kind)>;
  std::cerr << "[READ] consumer=" << consumer
            << " consumer_vid=" << ccell.value_id << " consumer_kind="
            << (ccell.production.kind == Kind::Build      ? "Build"
                : ccell.production.kind == Kind::Assemble ? "Assemble"
                                                          : "Other")
            << " consumer_depth=" << ccell.scope.path.size() << " value=" << vid
            << " source=" << r.source << " src_kind="
            << (src_cell.production.kind == Kind::Build      ? "Build"
                : src_cell.production.kind == Kind::Assemble ? "Assemble"
                                                             : "Leaf")
            << " src_depth=" << src_cell.scope.path.size()
            << " slices=" << r.slice.size();
  for (auto const& [pos, key] : r.slice) {
    std::cerr << " [pos" << pos << "@d" << key.depth << "#" << key.loop_slot;
    for (auto const& e : ctx)
      if (same_key(e.level.key(), key))
        std::cerr << "=" << e.range.first << ".." << e.range.second;
    std::cerr << "]";
  }
  std::cerr << " src_sliced={";
  for (auto const& [pos, key] : src_cell.sliced)
    std::cerr << "pos" << pos << "@d" << key.depth << "#" << key.loop_slot
              << " ";
  std::cerr << "}" << std::endl;
}

// ---------------------------------------------------------------------------
// SEQUANT_UT_BLOCK_DIAG
// ---------------------------------------------------------------------------

/// \brief `[BLOCK] ... SKIPPED WHOLE`: every production of this block (its
/// steps', its descendants' and its outputs') was already resident.
inline void dump_block_skipped(Index const& axis, std::size_t depth,
                               int loop_slot) {
  std::cerr << "[BLOCK] axis=" << toUtf8(axis.full_label())
            << " depth=" << depth << " slot=" << loop_slot
            << " SKIPPED WHOLE (every production resident)" << std::endl;
}

/// \brief `[BLOCK] ... BUILD`: one per-batch BuildStep of this block.
inline void dump_block_build(Index const& axis, std::size_t lo, std::size_t hi,
                             std::size_t vid, std::size_t cell,
                             std::size_t hash) {
  std::cerr << "[BLOCK] axis=" << toUtf8(axis.full_label()) << " batch=[" << lo
            << "," << hi << ") BUILD vid=" << vid << " cell=" << cell
            << " hash=" << hash << std::endl;
}

/// \brief `[BLOCK] ... ASSEMBLE`: one per-batch fold of this block's output.
inline void dump_block_assemble(Index const& axis, std::size_t lo,
                                std::size_t hi, std::size_t vid,
                                std::size_t cell, std::size_t src,
                                bool is_sum) {
  std::cerr << "[BLOCK] axis=" << toUtf8(axis.full_label()) << " batch=[" << lo
            << "," << hi << ") ASSEMBLE vid=" << vid << " cell=" << cell
            << " src=" << src << " kind=" << (is_sum ? "SUM" : "SCATTER")
            << std::endl;
}

/// \brief `[cell-registry] cache-halt skips`: how many cells the cache-halt
/// skip set covers.
inline void dump_cache_halt_skips(std::size_t n_skip, std::size_t n_cells) {
  std::cerr << "[cell-registry] cache-halt skips " << n_skip << " of "
            << n_cells << " cells" << std::endl;
}

/// \brief `[cell-registry] residency`: the end-of-call registry residency
/// split into live / persistent / root bytes.
inline void dump_registry_residency(std::size_t live, std::size_t persistent,
                                    std::size_t roots) {
  std::cerr << "[cell-registry] residency live=" << live
            << " persistent=" << persistent << " roots=" << roots << std::endl;
}

/// \brief `[cell-registry] resolver served`: how many operand reads the
/// table-driven resolver served this call.
inline void dump_resolver_served(std::size_t served) {
  std::cerr << "[cell-registry] resolver served " << served << " reads"
            << std::endl;
}

// ---------------------------------------------------------------------------
// SEQUANT_DUMP_ROOT_NORMS
// ---------------------------------------------------------------------------

/// \brief The norm of \p value, or -2.0 when it is absent or cannot be normed
/// (the dump must never throw out of a diagnostic).
template <typename Ptr>
[[nodiscard]] double dump_norm2(Ptr const& value) noexcept {
  try {
    if (value) return value->norm2();
  } catch (...) {
  }
  return -2.0;
}

/// \brief \c dump_norm2 of whatever \p fetch returns, with the FETCH ITSELF
/// inside the guard: a registry lookup (\c CellRegistry::peek) can throw on a
/// bad cell id, and computing the argument at the call site would put that
/// throw outside \c dump_norm2's own try -- i.e. a diagnostic that aborts the
/// run it is diagnosing.
template <typename F>
[[nodiscard]] double dump_norm2_of(F&& fetch) noexcept {
  try {
    return dump_norm2(fetch());
  } catch (...) {
  }
  return -2.0;
}

/// \brief `[build-norm]`: one root-scope build, by value hash.
inline void dump_build_norm(std::size_t vid, std::size_t hash, double norm) {
  std::cerr << "[build-norm] vid=" << vid << " hash=" << hash
            << " norm=" << std::setprecision(17) << norm << "\n";
}

/// \brief `[root-norm]`: one forest root's value at the combine, by node hash.
inline void dump_root_norm(std::size_t i, std::size_t hash, std::size_t vid,
                           std::size_t cell, int phase, double norm) {
  std::cerr << "[root-norm] i=" << i << " hash=" << hash << " vid=" << vid
            << " cell=" << cell << " phase=" << phase
            << " norm=" << std::setprecision(17) << norm << "\n";
}

}  // namespace sequant::eval::detail

#endif  // SEQUANT_CORE_EVAL_ORDERED_DUMP_HPP
