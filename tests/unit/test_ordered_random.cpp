// Randomized checks of the ordered schedule builder and the release plan:
// random value graphs (loop-local, reduced and carried modes on random loop
// chains, random operands, assembled reads) are built into schedules and
// cell tables, and
//  - the builder never fails to order the steps of a block level: the pass
//    levels are what keep those orderable (see forced_split_levels);
//  - the release plan never lets a cell go while a read of it is still to
//    come, and leaves nothing behind: a simulator of the executor's walk
//    (transcribed from ordered_executor.hpp, on the real CellRegistry and the
//    executor's own skip predicates) runs each legitimate table with and
//    without the plan and compares what the two did.
// Each case is seeded, so a failure prints the graph that produced it.

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/cell_registry.hpp>
#include <SeQuant/core/eval/cell_table_builder.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/release_plan.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using sequant::Index;
using sequant::ResultPtr;
using namespace sequant::eval;
namespace ed = sequant::eval::detail;

// A value of a random graph: the loop slots it is loop-local on, reduced
// over and carried on, its operands, and (slot, reader) pairs for a reader
// that takes it assembled over that slot.
struct V {
  std::vector<int> local, red, carried;
  std::vector<std::size_t> ops;
  std::vector<std::pair<int, std::size_t>> asmb;
  bool leaf = false;
};

Index ix_of(int s) { return Index{L"i_" + std::to_wstring(s + 1)}; }

std::string spec(std::vector<V> const& vs) {
  std::string s;
  for (std::size_t v = 0; v < vs.size(); ++v) {
    s += " v" + std::to_string(v) + "{L";
    for (int x : vs[v].local) s += std::to_string(x);
    s += " R";
    for (int x : vs[v].red) s += std::to_string(x);
    s += " C";
    for (int x : vs[v].carried) s += std::to_string(x);
    if (vs[v].leaf) s += " leaf";
    s += " <-";
    for (auto o : vs[v].ops) s += " " + std::to_string(o);
    for (auto const& [as, rd] : vs[v].asmb)
      s += " asm(" + std::to_string(as) + ",v" + std::to_string(rd) + ")";
    s += "}";
  }
  return s;
}

// The rich schedule and legality of a graph, as the boulevard and the
// legality pass would record them: one occurrence per reader, every mode
// home-sliced on its own slot, reduced slots stamped on every occurrence.
void build(std::vector<V> const& vs, RichSchedule& rich,
           LegalitySchedule& legality) {
  std::size_t const n = vs.size();
  std::vector<std::vector<std::size_t>> readers(n);
  for (std::size_t v = 0; v < n; ++v)
    for (std::size_t o : vs[v].ops) readers[o].push_back(v);
  auto const pt = [](std::size_t v, std::size_t k) {
    return 100 * (v + 1) + k;
  };
  for (std::size_t v = 0; v < n; ++v) {
    std::vector<int> all = vs[v].local;
    all.insert(all.end(), vs[v].carried.begin(), vs[v].carried.end());
    std::sort(all.begin(), all.end());
    ValueCell vc{};
    vc.value_id = v;
    vc.hash = 9000 + v;
    vc.first_use = v;
    vc.last_use = v;
    vc.is_leaf = vs[v].leaf;
    std::vector<std::pair<std::size_t, std::size_t>> occs;
    if (readers[v].empty())
      occs.push_back({pt(v, 0), pt(v, 0)});
    else
      for (std::size_t k = 0; k < readers[v].size(); ++k)
        occs.push_back({pt(v, k), pt(readers[v][k], 0)});
    for (auto const& [p, cp] : occs) {
      OccurrenceRec o{};
      o.point = p;
      o.consumer_point = cp;
      for (int s : all) {
        o.carried.push_back(ix_of(s));
        o.home.push_back(ix_of(s));
        o.loop_slot.push_back(s);
      }
      for (int s : vs[v].red) o.reduced_slot.push_back({ix_of(s), s});
      vc.occurrences.push_back(std::move(o));
    }
    vc.carried = vc.occurrences.front().carried;
    rich.cells.push_back(std::move(vc));

    CellLegality cl;
    cl.hash = 9000 + v;
    for (int s : vs[v].local) {
      AxisClass ac;
      ac.axis = ix_of(s);
      ac.role = LoopRole::LoopLocal;
      for (auto const& [as, rd] : vs[v].asmb)
        if (as == s)
          for (std::size_t k = 0; k < readers[v].size(); ++k)
            if (readers[v][k] == rd)
              ac.per_occurrence.push_back(
                  {pt(v, k), OccurrenceRole::Assembled});
      cl.per_axis.push_back(ac);
    }
    for (int s : vs[v].carried) {
      AxisClass ac;
      ac.axis = ix_of(s);
      ac.role = LoopRole::LoopCarried;
      cl.per_axis.push_back(ac);
      cl.forced_split_axes.push_back(ix_of(s));
    }
    for (int s : vs[v].red) {
      AxisClass ac;
      ac.axis = ix_of(s);
      ac.role = LoopRole::Reduction;
      cl.per_axis.push_back(ac);
    }
    legality.cells.push_back(cl);
  }
  // loop kinds, as the opens record them
  for (auto const& vc : rich.cells)
    for (auto const& occ : vc.occurrences) {
      for (std::size_t p = 0; p < occ.loop_slot.size(); ++p)
        rich.loop_kind.emplace(
            std::make_pair(std::wstring{occ.carried[p].basis().base_key()},
                           occ.loop_slot[p]),
            sequant::BatchModeType::External);
      for (auto const& [m, slot] : occ.reduced_slot)
        rich.loop_kind.emplace(
            std::make_pair(std::wstring{m.basis().base_key()}, slot),
            sequant::BatchModeType::Contracted);
    }
}

// A random graph. mode bits: 1 = one loop chain (else one to three),
// 4 = bigger shapes, 8 = a value spanning each chain, 16 = leaves,
// 32 = reduce every unread loop-local value, 64 = many root values.
std::vector<V> random_graph(std::mt19937_64& rng, int mode) {
  auto const rnd = [&](std::uint64_t n) { return rng() % n; };
  bool const big = mode & 4;
  std::size_t const nch = (mode & 1) ? 1 : 1 + rnd(3);
  std::vector<std::vector<int>> chains(nch);
  int next = 0;
  for (auto& ch : chains) {
    int const len = 1 + static_cast<int>(rnd(big ? 5 : 3));
    for (int k = 0; k < len && next < 8; ++k) ch.push_back(next++);
    if (ch.empty()) ch.push_back(next++);
  }
  std::size_t const n = big ? 6 + rnd(8) : 4 + rnd(6);
  std::vector<V> vs(n);
  std::uint64_t const rootw = (mode & 64) ? 2 * nch + 1 : 1;
  for (std::size_t v = 0; v < n; ++v) {
    V& x = vs[v];
    std::uint64_t const where = rnd(nch * 4 + rootw);
    if (where < nch * 4) {
      auto const& ch = chains[where % nch];
      int const L = static_cast<int>(ch.size());
      int const d = static_cast<int>(rnd(L));
      std::uint64_t const kind = rnd(15);
      if (kind >= 12) {  // escape on d, loop-local on deeper levels too
        (kind == 12 ? x.red : x.carried).push_back(ch[d]);
        for (int k = 0; k < L; ++k)
          if (k != d && (k > d ? rnd(3) != 0 : rnd(3) == 0))
            x.local.push_back(ch[k]);
      } else if (kind < 6) {  // all-local prefix 0..d
        for (int k = 0; k <= d; ++k)
          if (k == d || rnd(8) != 0) x.local.push_back(ch[k]);
      } else if (kind < 9) {  // reduce d0..d1, local prefix
        int const d0 = static_cast<int>(rnd(d + 1));
        for (int k = 0; k < d0; ++k)
          if (rnd(8) != 0) x.local.push_back(ch[k]);
        int const d1 = d0 + static_cast<int>(rnd(L - d0));
        for (int k = d0; k <= d1; ++k) x.red.push_back(ch[k]);
      } else if (kind < 11) {  // carried on d, local on a subset of the rest
        x.carried.push_back(ch[d]);
        for (int k = 0; k < L; ++k)
          if (k != d && rnd(2)) x.local.push_back(ch[k]);
      } else {  // carried on d, local prefix above, reduction below
        for (int k = 0; k < d; ++k) x.local.push_back(ch[k]);
        x.carried.push_back(ch[d]);
        for (int k = d + 1; k < L; ++k)
          if (rnd(2)) x.red.push_back(ch[k]);
      }
    }
    if (v > 0) {
      std::uint64_t const r = rnd(24);
      std::uint64_t const nops = r < 3 ? 0 : r < 16 ? 1 : r < 22 ? 2 : 3;
      for (std::uint64_t k = 0; k < nops; ++k) {
        std::size_t const o =
            rnd(2) ? rnd(v) : v - 1 - rnd(std::min<std::size_t>(v, 3));
        if (std::find(x.ops.begin(), x.ops.end(), o) == x.ops.end())
          x.ops.push_back(o);
      }
    }
    if ((mode & 16) && x.ops.empty() && x.red.empty() && x.carried.empty() &&
        rnd(3) == 0)
      x.leaf = true;
  }
  if (mode & 8)
    for (auto const& ch : chains) {
      V s;
      s.local = ch;
      vs.push_back(s);
    }
  if (mode & 32) {
    std::size_t const m = vs.size();
    std::vector<bool> read(m, false);
    for (auto const& x : vs)
      for (auto o : x.ops) read[o] = true;
    for (std::size_t v = 0; v < m; ++v)
      if (!read[v] && !vs[v].local.empty()) {
        V t;
        t.red = vs[v].local;
        t.ops.push_back(v);
        vs.push_back(t);
      }
  }
  // a reader not sliced on an operand's loop-local slot reads it assembled
  for (std::size_t v = 0; v < vs.size(); ++v)
    for (std::size_t o : vs[v].ops)
      for (int s : vs[o].local) {
        auto const has = [&](std::vector<int> const& c) {
          return std::find(c.begin(), c.end(), s) != c.end();
        };
        if (!(has(vs[v].local) || has(vs[v].red))) vs[o].asmb.push_back({s, v});
      }
  return vs;
}

sequant::BatchPolicy random_policy() {
  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](Index const& ix) {
    return ix.basis().base_key() == L"i";
  };
  return policy;
}

std::uint64_t mix(std::uint64_t a, std::uint64_t b, std::uint64_t c = 0,
                  std::uint64_t d = 0) {
  std::uint64_t h = 0x9E3779B97F4A7C15ull;
  for (std::uint64_t x : {a, b, c, d}) {
    h ^= x + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
    h *= 0xBF58476D1CE4E5B9ull;
    h ^= h >> 31;
  }
  return h;
}

// The executor's walk over a cell table on the real CellRegistry with dummy
// results: run_ordered_schedule_pre_results, run_ordered_contracted_block
// and compute_cell transcribed, the skip predicates the executor's own. With
// use_plan the release plan's releases and takes are applied; without, a
// cell is only ever emptied by clear_bound_to, which is the reference every
// production and read count is compared against.
struct Sim {
  CellTable const& t;
  OrderedSchedule const& sched;
  ReleasePlan const& plan;
  bool use_plan;
  std::uint64_t seed;
  bool custom_on, empty_on;
  std::vector<std::size_t> const& roots;  // forest-root value ids
  PersistentValueStore& store;
  std::vector<std::vector<std::size_t>> const& reads_of;
  std::shared_ptr<sequant::eval::dryrun::CostModel const> cm;

  std::vector<std::string> errs = {};
  std::vector<std::size_t> prod = {}, leaf_fetch = {};
  std::size_t reads = 0, takes = 0;
  std::optional<CellRegistry> reg = {};
  std::string where = {};

  ResultPtr fresh() {
    return std::make_shared<sequant::eval::dryrun::ResultDryRun>(
        sequant::container::svector<Index>{}, cm);
  }
  void err(std::string const& m) {
    if (errs.size() < 6) errs.push_back(m + " @" + where);
  }
  std::size_t nbatches(LoopKey const& k) const {
    return 1 + mix(seed, k.depth, k.loop_slot, 11) % 3;
  }
  bool custom(CellId c) const { return custom_on && mix(seed, c, 99) % 7 == 0; }

  // compute_cell: one production of cell c.
  ResultPtr compute(CellId c) {
    ++prod[c];
    if (custom(c)) return fresh();  // the custom-evaluator seam: no reads
    std::vector<ResultPtr> legs;
    for (std::size_t r : reads_of[c]) {
      CellId const src = t.reads[r].source;
      bool const leaf = t.cells[src].production.kind == ProductionKind::Leaf;
      if (!reg->peek(src)) {
        if (!leaf) {
          err("read of absent " + ed::cell_str(t, src) + " by " +
              ed::cell_str(t, c) + " read#" + std::to_string(r));
          continue;
        }
        ++leaf_fetch[src];
        reg->set(src, fresh());  // record_leaf
      }
      bool const last = use_plan && plan.last_read[r];
      ++reads;
      if (last) {
        ResultPtr v = reg->take(src);
        ++takes;
        // the reader is now the sole holder, besides its own earlier legs
        long others = 0;
        for (auto const& l : legs) others += (l == v);
        if (!t.cells[src].persistent && v.use_count() != 1 + others)
          err("take of " + ed::cell_str(t, src) + " not sole-held");
        legs.push_back(std::move(v));
      } else {
        legs.push_back(reg->read(src));
      }
    }
    return fresh();
  }

  void set(CellId c, ResultPtr v) {
    try {
      reg->set(c, std::move(v));
    } catch (std::exception const& e) {
      err(std::string("set throws: ") + e.what());
    }
  }

  void run_block(ScopeBlock const& block, CellScope const& parent_scope,
                 sequant::container::vector<char> const& skip) {
    CellScope block_scope = parent_scope;
    block_scope.path.push_back({block.level.key(), block.latitude_ordinal});
    sequant::container::vector<char> local_skip_storage;
    sequant::container::vector<char> const* skip_p = &skip;
    {
      bool seeded = false;
      for (CellId c = 0; c < t.cells.size(); ++c) {
        if (skip[c] ||
            !ed::ordered_visit_skip_seedable(t.cells[c], parent_scope))
          continue;
        if (!reg->peek(c)) continue;
        if (!seeded) {
          local_skip_storage = skip;
          seeded = true;
        }
        local_skip_storage[c] = 1;
      }
      if (seeded) {
        local_skip_storage = ed::ordered_skip_closure(
            t, std::move(local_skip_storage), plan.kept);
        skip_p = &local_skip_storage;
      }
    }
    sequant::container::vector<char> const& vskip = *skip_p;
    if (ed::ordered_block_fully_skipped(*reg, vskip, block, parent_scope))
      return;
    std::unordered_set<std::size_t> built_here;
    for (Step const& s : block.steps)
      if (auto const* b = std::get_if<BuildStep>(&s.value))
        built_here.insert(b->value_id);
    std::vector<CellId> out_cells(block.outputs.size(), 0);
    std::vector<char> out_skip(block.outputs.size(), 0);
    for (std::size_t k = 0; k != block.outputs.size(); ++k) {
      auto const a =
          reg->assemble_cell_at(block.outputs[k].first, parent_scope);
      if (!a) {
        err("no Assemble cell for an output");
        return;
      }
      out_cells[k] = *a;
      out_skip[k] = vskip[*a] ||
                    (t.cells[*a].produce_if_absent && reg->peek(*a) != nullptr);
    }
    std::size_t const nb = nbatches(block.level.key());
    std::vector<ResultPtr> acc(block.outputs.size());
    // some batches are empty, never the one at index keep_b
    std::size_t const keep_b =
        mix(seed, block.level.depth, block.level.loop_slot, 5) % nb;
    for (std::size_t b = 0; b < nb; ++b) {
      if (empty_on && b != keep_b &&
          mix(seed, block.level.depth * 131 + block.level.loop_slot, b, 7) %
                  4 ==
              0)
        continue;
      reg->clear_bound_to(block.level.key());
      for (std::size_t at = 0; at != block.steps.size(); ++at) {
        Step const& step = block.steps[at];
        ed::ReleaseAfter const release_after{*reg, plan, block, at};
        if (auto const* build = std::get_if<BuildStep>(&step.value)) {
          auto const bc = reg->build_cell_at(build->value_id, block_scope);
          if (!bc) {
            err("no Build cell for a step");
            continue;
          }
          bool const pia_resident =
              t.cells[*bc].produce_if_absent && reg->peek(*bc) != nullptr;
          if (pia_resident || vskip[*bc]) continue;
          set(*bc, compute(*bc));
        } else {
          run_block(std::get<ScopeBlock>(step.value), block_scope, vskip);
        }
      }
      for (std::size_t k = 0; k != block.outputs.size(); ++k) {
        ed::ReleaseAfter const release_after{*reg, plan, block,
                                             block.steps.size() + k};
        if (out_skip[k]) continue;
        auto const vid = block.outputs[k].first;
        CellId const src = t.cells[out_cells[k]].production.source;
        TableCell const& s = t.cells[src];
        if (s.production.kind == ProductionKind::Build &&
            s.scope == block_scope && !built_here.count(vid))
          set(src, compute(src));
        bool const takes_it = use_plan && plan.last_assemble[out_cells[k]] &&
                              !t.cells[src].persistent;
        if (!reg->peek(src)) {
          err("Assemble read of absent " + ed::cell_str(t, src) + " by " +
              ed::cell_str(t, out_cells[k]));
          continue;
        }
        ++reads;
        ResultPtr part = takes_it ? reg->take(src) : reg->read(src);
        if (takes_it) {
          ++takes;
          if (part.use_count() != 1)
            err("Assemble take of " + ed::cell_str(t, src) + " not sole-held");
        }
        if (!acc[k]) acc[k] = fresh();
      }
    }
    for (std::size_t k = 0; k != block.outputs.size(); ++k) {
      if (out_skip[k]) continue;
      if (!acc[k]) {
        err("zero batches realized");
        continue;
      }
      set(out_cells[k], std::move(acc[k]));
    }
  }

  void evaluate() {
    prod.assign(t.cells.size(), 0);
    leaf_fetch.assign(t.cells.size(), 0);
    CellRegistryHooks hooks;
    hooks.persistent = &store;
    hooks.hash_of = [](std::size_t vid) { return vid; };
    hooks.strict_fill_once = true;
    reg.emplace(t, std::move(hooks));
    reg->seed_persistent();
    auto const skip = ed::ordered_cache_halt_skip(t, *reg, plan.kept);
    for (std::size_t at = 0; at != sched.root.steps.size(); ++at) {
      Step const& step = sched.root.steps[at];
      ed::ReleaseAfter const release_after{*reg, plan, sched.root, at};
      where = "root step " + std::to_string(at);
      if (auto const* build = std::get_if<BuildStep>(&step.value)) {
        auto const rc = reg->build_cell_at(build->value_id, CellScope{});
        if (!rc) {
          err("no root Build cell");
          continue;
        }
        if (skip[*rc] || (t.cells[*rc].produce_if_absent && reg->peek(*rc)))
          continue;
        set(*rc, compute(*rc));
      } else {
        run_block(std::get<ScopeBlock>(step.value), CellScope{}, skip);
      }
    }
    where = "end";
    std::unordered_set<CellId> root_cells;
    for (std::size_t vid : roots) {
      auto const c = reg->cell_of(vid, CellScope{});
      if (!c) {
        err("forest root v" + std::to_string(vid) + " has no root form");
        continue;
      }
      root_cells.insert(*c);
      if (!reg->peek(*c))
        err("forest root " + ed::cell_str(t, *c) + " holds no result");
    }
    reg->for_each_live([&](CellId c, ResultPtr const&) {
      if (!t.cells[c].persistent && !root_cells.count(c))
        err("left over: " + ed::cell_str(t, c));
    });
  }
};

// Three consecutive evaluations of one table on one persistent store (cold,
// warm, and with half the store lost), with and without the plan. Returns
// the discrepancies.
std::vector<std::string> simulate(CellTable const& t,
                                  OrderedSchedule const& sched,
                                  std::vector<std::size_t> const& roots,
                                  std::uint64_t seed, bool custom_on,
                                  bool empty_on) {
  std::vector<std::string> out;
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      sequant::eval::dryrun::SizeRegime{});
  sequant::container::svector<CellId> keep;
  {
    CellRegistry probe(t);
    for (std::size_t vid : roots)
      if (auto const c = probe.cell_of(vid, CellScope{})) keep.push_back(*c);
  }
  ReleasePlan const plan = build_release_plan(t, sched.root, keep);
  ReleasePlan none;
  none.last_read.assign(t.reads.size(), 0);
  none.last_assemble.assign(t.cells.size(), 0);
  none.kept = plan.kept;
  for (auto const& [block, slots] : plan.after)
    none.after[block].resize(slots.size());
  std::vector<std::vector<std::size_t>> reads_of(t.cells.size());
  for (std::size_t r = 0; r < t.reads.size(); ++r)
    reads_of[t.reads[r].consumer].push_back(r);
  PersistentValueStore store_p, store_n;
  for (int ev = 0; ev < 3; ++ev) {
    if (ev == 2) {  // the store lost a random subset
      std::vector<std::size_t> gone;
      store_p.for_each([&](std::size_t h, ResultPtr const&) {
        if (mix(seed, h, 4242) % 2) gone.push_back(h);
      });
      for (auto h : gone) {
        store_p.erase(h);
        store_n.erase(h);
      }
    }
    Sim p{t,        sched, plan,    true,     seed, custom_on,
          empty_on, roots, store_p, reads_of, cm};
    Sim n{t,        sched, none,    false,    seed, custom_on,
          empty_on, roots, store_n, reads_of, cm};
    p.evaluate();
    n.evaluate();
    std::string const ev_s = " [evaluation " + std::to_string(ev) + "]";
    for (auto const& e : p.errs) out.push_back("with the plan: " + e + ev_s);
    for (auto const& e : n.errs)
      if (e.find("left over") == std::string::npos)
        out.push_back("without the plan: " + e + ev_s);
    for (CellId c = 0; c < t.cells.size(); ++c) {
      if (p.prod[c] != n.prod[c])
        out.push_back("productions of " + ed::cell_str(t, c) + " differ: " +
                      std::to_string(p.prod[c]) + " with the plan, " +
                      std::to_string(n.prod[c]) + " without" + ev_s);
      if (p.leaf_fetch[c] != n.leaf_fetch[c])
        out.push_back("leaf fetches of " + ed::cell_str(t, c) + " differ" +
                      ev_s);
    }
    if (p.reads != n.reads)
      out.push_back("read counts differ: " + std::to_string(p.reads) +
                    " with the plan, " + std::to_string(n.reads) + " without" +
                    ev_s);
  }
  return out;
}

// The validated cell table of a built schedule, or nullopt when the table
// is not legitimate (violations, or unresolved sliced positions): the
// generator does not guarantee a consistent legality record.
std::optional<CellTable> legitimate_table(
    RichSchedule const& rich, OrderedSchedule const& sched,
    ed::OrderedScheduleDepGraph const& g) {
  auto const sma = compute_sliced_mode_assignment(sched, rich);
  CellTableInputs in;
  in.ordered = &sched;
  in.rich = &rich;
  in.sliced = &sma;
  in.sliced_modes_of = [&](std::size_t vid) {
    return rich.cells[vid].occurrences.front().home;
  };
  in.volatile_of = [](std::size_t) { return false; };
  in.operands_of = [&](std::size_t vid) {
    auto const it = g.depends_on.find(vid);
    return it == g.depends_on.end() ? sequant::container::svector<std::size_t>{}
                                    : it->second;
  };
  auto table = build_cell_table(in);
  if (!validate_cell_table(table, sched.root).empty() ||
      !table.unresolved.empty())
    return std::nullopt;
  return table;
}

#ifdef SEQUANT_SKIP_LONG_TESTS
constexpr std::size_t n_graphs_per_mode = 60;
#else
constexpr std::size_t n_graphs_per_mode = 1500;
#endif
constexpr int modes[] = {0, 1, 4 | 8, 16 | 32, 64, 4 | 64};

}  // namespace

TEST_CASE("random schedules: every block level is orderable",
          "[ordered-schedule][random]") {
  auto const policy = random_policy();
  std::size_t built = 0, refused = 0;
  for (int mode : modes) {
    std::mt19937_64 rng(1000 + mode);
    for (std::size_t it = 0; it < n_graphs_per_mode; ++it) {
      auto const vs = random_graph(rng, mode);
      INFO("mode " << mode << " graph " << it << ":" << spec(vs));
      RichSchedule rich;
      LegalitySchedule legality;
      build(vs, rich, legality);
      try {
        auto const sched =
            build_ordered_schedule(rich, legality, policy, {L"i"});
        CHECK(well_formed(sched));
        ++built;
      } catch (sequant::Exception const& e) {
        // The generator can ask for the unsupported (a reader outside its
        // operand's nest, say), which the builder refuses by name. It must
        // never fail to order a level: the pass levels keep every level
        // orderable.
        std::string const what = e.what();
        CHECK(what.find("Impossible dependencies") == std::string::npos);
        ++refused;
      }
    }
  }
  CHECK(built > refused);
}

TEST_CASE("random schedules: the release plan never releases early",
          "[ordered-schedule][release-plan][random]") {
  auto const policy = random_policy();
  std::size_t simulated = 0;
  for (int mode : modes) {
    std::mt19937_64 rng(2000 + mode);
    for (std::size_t it = 0; it < n_graphs_per_mode; ++it) {
      auto const vs = random_graph(rng, mode);
      INFO("mode " << mode << " graph " << it << ":" << spec(vs));
      RichSchedule rich;
      LegalitySchedule legality;
      build(vs, rich, legality);
      std::optional<OrderedSchedule> sched;
      try {
        sched = build_ordered_schedule(rich, legality, policy, {L"i"});
      } catch (sequant::Exception const&) {
        continue;
      }
      auto const g = ed::ordered_schedule_dep_graph(rich);
      auto const table = legitimate_table(rich, *sched, g);
      if (!table) continue;
      // forest roots: the non-leaf values nothing reads
      std::vector<char> is_read(rich.cells.size(), 0);
      for (auto const& [v, ops] : g.depends_on)
        for (auto o : ops) is_read[o] = 1;
      std::vector<std::size_t> roots;
      for (std::size_t v = 0; v < rich.cells.size(); ++v)
        if (!is_read[v] && !rich.cells[v].is_leaf) roots.push_back(v);
      // three batch-count draws; the custom-evaluator seam and empty batches
      // on the later two
      for (std::uint64_t k = 0; k < 3; ++k) {
        auto const problems =
            simulate(*table, *sched, roots, 7919 * it + k, k >= 1, k == 2);
        for (auto const& p : problems) UNSCOPED_INFO(p);
        CHECK(problems.empty());
      }
      ++simulated;
    }
  }
  CHECK(simulated > 0);
}
