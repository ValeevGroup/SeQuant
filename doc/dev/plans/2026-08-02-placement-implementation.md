# Home-scope placement — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the placement-as-register-allocation design
(`doc/dev/specs/2026-08-02-home-scope-placement-design.md`), replacing the
caching-veto framing with correct cell placement under a peak budget.

**Architecture:** Five phases, each a working, testable increment; each later
phase gets its own detailed plan once the prior lands. Phase 1 (this plan,
detailed) corrects the peak *measurement* to the true co-resident sum — the
oracle everything else is validated against; Phases 2-5 are a roadmap here.

**Tech Stack:** C++20, SeQuant dry-run backend (`cost_profile()`), the real eval
loop (`eval.hpp`), `cache_manager.hpp`; Catch2 unit tests; Boost rational (for
`W` in a later phase).

## Global Constraints

- **No en-dashes** (U+2013) anywhere — the repo has a `forbid-en-dashes` commit
  hook; use `--` or `-`. Em-dashes are fine.
- **clang-format** touched C/C++ with `/opt/homebrew/opt/llvm/bin/clang-format
  --style=file -i` before committing.
- **No `Co-Authored-By` trailers**; keep dependency-pin bumps in their own commit.
- **Perfect-CSE default must stay byte-identical**: any keying/router default
  must reproduce today's behavior when no split/placement decision is made.
- **The replay is the ground-truth peak oracle**; the static model (Phase 3+)
  must be validated against it, never the reverse.

---

## Roadmap (phase = future plan unless detailed below)

- **Phase 1 (detailed here) — peak oracle.** Correct `cost_profile()`'s
  `peak_bytes` from `max(scratch_hwmark, cache_hwmark)` to the instant-resolved
  **sum** of co-resident residency across the live scope chain (spec §7c / O3b).
  Fixes the documented lower-bound under-count; nothing else depends on the
  router, so it stands alone.
- **Phase 2 — home_scope + router seed** (spec §7a, §7d / O1, O5). Compute
  `home_scope(value) = deepest scope over (sliced_modes ∪ demoted_external_modes)`
  and materialize the explicit `{value, use-site} -> (home, split=0)` router as
  the perfect-CSE seed; reads route via the router and reuse the existing
  `(use - home) INTERSECT carried` slicer instead of the `parent_` walk +
  `hops`. Byte-identical (seed == perfect CSE), but replaces the split authority.
- **Phase 3 — static peak profile** (spec §7c / O3a). The weighted-interval
  sweep (`[first-use, last-use]` per cell x footprint) over the router+schedule,
  validated to equal the Phase-1 replay oracle.
- **Phase 4 — the O2 greedy** (spec §7b). Spill loop: shrink/evict cells alive at
  the binding peak point by ΔPeak/ΔRecompute against the Phase-3 profile,
  emitting a router with non-trivial homes/splits.
- **Phase 5 — feedback** (spec §7b, O6). Detect-and-report the two infeasible
  modes (factorization-inherent vs. re-batch-needed); optional re-batch hint.

---

## Phase 1: correct the replay peak to the co-resident sum

**Why.** `cost_profile()`'s `peak_bytes` is
`max(profile.peak_bytes, peak.load(), cache.working_set_hwmark())`
(`cost_profile.hpp` replay loop) — the max of two *independent* high-watermarks
(the batched-scratch PeakSink and the outer cache), so when a persistent
cross-term cell co-resides with a batched-inner transient the two are additive
and the max under-reports (the field's own doc calls it a lower bound). The true
peak is `max over instants of (sum of all live residency across the scope
chain)`. Note `max(hwmark_a, hwmark_b)` is a lower bound and `hwmark_a + hwmark_b`
is an upper bound (the two peaks may occur at different instants); the correct
value is the sum taken *per instant*, then maxed — which requires summing the
scope chain at each `note_working_set` point, not combining two end-of-run
hwmarks.

**Files:**
- Modify: `SeQuant/core/eval/cache_manager.hpp` (add per-cache current-residency
  + scope-chain sum; feed the chain sum into the hwmark)
- Modify: `SeQuant/core/eval/eval.hpp` (the `note_working_set` call sites pass
  the chain-summed working set)
- Modify: `SeQuant/core/eval/backends/dryrun/cost_profile.hpp` (`peak_bytes` doc;
  the fold no longer needs the separate `working_set_hwmark()` term once the
  chain sum subsumes it)
- Test: `tests/unit/test_cache_manager.cpp` (chain-residency unit), and
  `tests/unit/test_eval_dryrun.cpp` (the `[dryrun][peak]` / `[cost_profile]`
  co-resident assertion + witness re-baseline)

**Interfaces:**
- Produces: `size_t CacheManager::current_residency() const` — sum over *alive*
  entries of `entry.size_in_bytes()` in this cache only.
- Produces: `size_t CacheManager::chain_residency() const` — `current_residency()
  + (parent_ ? parent_->chain_residency() : 0)`.
- Consumes (Phase 3): the corrected `CostProfile::peak_bytes` as the oracle.

- [ ] **Step 0: Read the current peak path.** Read `cache_manager.hpp` lines
  ~185-305 (`working_set_hwmark_`, `note_working_set`, `reset`) and
  `eval.hpp:576-585, 1867-1874` (the `note_working_set(hwmark)` sites; `hwmark =
  log::bytes(cache, post).value`) and `cost_profile.hpp` replay fold
  (`profile.peak_bytes = std::max({..., peak.load(),
  cache.working_set_hwmark()})`). Confirm: each `CacheManager` tracks *its own*
  hwmark; the scratch PeakSink folds the scratch hwmark; the outer cache hwmark
  is separate; they are combined by `max` at the end.

- [ ] **Step 1: Write the failing test for `current_residency()`.**
  In `tests/unit/test_cache_manager.cpp`, in a new `TEST_CASE("cache_manager
  residency", "[cache_manager]")`: build a manager over two nodes, `store()` data
  into both, and assert
  `man.current_residency() == man.entry_size_in_bytes(k0) +
  man.entry_size_in_bytes(k1)` while both are alive, and that it drops as entries
  are drained (mirroring the existing `alive and entry_size_in_bytes` section).

- [ ] **Step 2: Run it, verify it fails** (`current_residency` undefined).
  Run: `ninja -C cmake-build-release unit_tests-sequant -j6 && \
  ./cmake-build-release/tests/unit/unit_tests-sequant "cache_manager residency"`.
  Expected: compile error / FAIL.

- [ ] **Step 3: Implement `current_residency()` + `chain_residency()`** in
  `cache_manager.hpp` (public):
  ```cpp
  [[nodiscard]] size_t current_residency() const noexcept {
    size_t s = 0;
    for (auto const& [k, e] : cache_map_)
      if (e.alive()) s += e.size_in_bytes();
    return s;
  }
  [[nodiscard]] size_t chain_residency() const noexcept {
    return current_residency() + (parent_ ? parent_->chain_residency() : 0);
  }
  ```

- [ ] **Step 4: Run the test, verify it passes.** Same command as Step 2.
  Expected: PASS.

- [ ] **Step 5: Commit.**
  ```bash
  git add SeQuant/core/eval/cache_manager.hpp tests/unit/test_cache_manager.cpp
  git commit -m "eval: add CacheManager current/chain residency accessors"
  ```

- [ ] **Step 6: Thread the chain sum into the per-op hwmark ONLY.** The fix is
  local to the `note_working_set` call sites; **leave both folds unchanged** (see
  Step 7 for why). In `eval.hpp`, at *every* `cache.note_working_set(hwmark)` site
  (576-585 and the others found in Step 0), add the ancestors' live residency to
  the local `hwmark` before the call, so the value recorded into that cache's
  `working_set_hwmark_` is the instant co-resident total:
  `hwmark += (cache.parent() ? cache.parent()->chain_residency() : 0);`
  For the OUTER cache (no parent) this adds 0 -> its hwmark is unchanged. For a
  SCRATCH cache (parent = outer) this adds the outer's current residency at that
  instant -> the scratch's `working_set_hwmark_` becomes
  `max over the scratch's life of (scratch_now + outer_now)` = the true
  co-resident peak. Do NOT touch the scatter-branch fold at 1867-1874: it folds
  `bs.cache.working_set_hwmark()` (the scratch's MAX over its life) into the
  PeakSink, and that max is now already chain-inclusive from these per-op calls;
  folding `chain_residency()` (an end-of-scatter instant) there instead would
  miss intermediate peaks -- keep `working_set_hwmark()`.

- [ ] **Step 7: Leave the `cost_profile()` fold as-is (verify only).** Do NOT
  change `profile.peak_bytes = std::max({profile.peak_bytes, peak.load(),
  double(cache.working_set_hwmark())})`. It is already correct once Step 6 lands:
  `peak.load()` (the PeakSink) now carries the scratch co-resident peak, and the
  `cache.working_set_hwmark()` term is the outer-alone peak, still needed for the
  no-batching case (where `peak.load()` is 0 because no scatter fires). The `max`
  of the two is the true peak in both regimes -- dropping the outer term would
  under-report the unbatched case. Only UPDATE the `peak_bytes` doc comment
  (lines ~189-201): it is now the true co-resident peak (the per-op hwmark sums
  the live scope chain), not the old `max`-of-independent-hwmarks lower bound.

- [ ] **Step 8: Write the failing co-resident test.** In `test_eval_dryrun.cpp`,
  adapt the `[dryrun][cost_profile]` fixture (or the `[dryrun][peak]` one) into a
  case where a persistent cross-term intermediate is alive while a batched-inner
  transient is alive: assert the new `cp.peak_bytes` equals the *sum* of the two
  live footprints at the co-resident instant, which is strictly greater than the
  old `max` of the two separate hwmarks. (Reuse the existing giant-term forest;
  the persistent gC-class intermediate co-resides with the batched W scratch.)
  Compute the expected sum from the two nodes' `memsize` at their realized
  extents, as the `[dryrun][cache]` test already sizes the giant.

- [ ] **Step 9: Run; verify PASS after Steps 6-7** (the summed peak is now
  reported). Run the `[dryrun]` filter.

- [ ] **Step 10: Re-baseline the witnesses.** The corrected (higher) peak shifts
  the documented `[.][dryrun-occ-veto]` / `[.][dryrun-extmode-avoidable]` peak
  figures (e.g. `~443.6 GB`, `~6026 GB`). Run each at nterms=55, record the new
  peaks, and update the doc figures + the `documented-RED research target (~X GB)`
  comments. These are diagnostics, not gates (the `peak < 100` CHECK stays RED);
  do **not** invent numbers — measure. Note in each that the peak rose because it
  is now the true co-resident sum, not the old lower bound.

- [ ] **Step 11: Full suite green.** `ctest --test-dir cmake-build-release` for
  the eval/cache/dryrun tests (`[cache_manager]`, `[lifetime_mask]`, `[dryrun]`,
  `[eval]`); the `[.]` witnesses' *acceptance* gates (role stamps, `abs < 0.02`)
  must stay green, only the documented-RED peak/avoidable targets move.

- [ ] **Step 12: clang-format + commit.**
  ```bash
  /opt/homebrew/opt/llvm/bin/clang-format --style=file -i \
    SeQuant/core/eval/eval.hpp SeQuant/core/eval/cache_manager.hpp \
    SeQuant/core/eval/backends/dryrun/cost_profile.hpp tests/unit/test_eval_dryrun.cpp
  git add -A && git commit -m "eval: peak_bytes = true co-resident sum, not max(scratch,cache)"
  ```

## Self-review checklist (run after Phase 1)

- The `max(a,b) -> per-instant sum` change is correct: the summed value is
  recorded *at each op* (an instant), then maxed, so it is neither the
  lower-bound `max(hwmarks)` nor the upper-bound `sum(hwmarks)`.
- No `[.]` acceptance gate (`contracted_occ == 0`, `external > 0`,
  `abs(ext-con) < 0.02`) moved; only documented-RED peak targets re-baselined
  from *measurement*.
- Perfect-CSE / OFF path unchanged: with no batching every scratch is absent, so
  `chain_residency()` == the single cache's residency == today's value.

## Notes for later phases (not tasks yet)

- Phase 2's router key is `{value, use-site}`; the use-site is the eval-tree node
  occurrence (pointer/position) — O5b / O1 residuals.
- Phase 3's sweep is validated to equal this phase's replay `peak_bytes` on the
  same forests.
- Phase 4 needs Phase 3's per-point live set (the sweep's argmax) for the spill
  candidates.

---

## Phase 2 -- Home-scope placement: the router as a READ+STORE override seam

### Goal
Land an explicit placement router as an override seam that BOTH the Enter-stage read and the `place_at_this_level` store consult: `router: occurrence_key -> optional<(home_target, split_index=0)>`. In Phase 2 the router is seeded EMPTY, so both sides derive home exactly as today and behavior is byte-identical. The occurrence key is a proper, batching-aware, symmetry-intrinsic identity (a canonicalized colored sub-network), NOT the raw eval-node identity. Overrides are populated later (Phase 3/4); Phase 2 delivers the correctly-keyed seam and proves it with one hand-injected relocation.

### Architecture
Values are canonical `EvalExprDryRun` tree nodes; a per-scope `CacheManager` (the cell store) is keyed by structural hash and chained by `parent_` along the batch-loop nest. Phase 2 adds `PlacementRouter`, keyed by a `TensorNetwork::SlotCanonicalizationMetadata` computed from the node's sub-expression with the in-scope batched indices as `named_indices`. The read side (`access_at`/`slice_to_use`) and the store side (`place_at_this_level`'s `rl+1` level choice) each consult the router; an absent override falls to today's derivation. Empty router => both sides derive as today => byte-identical.

### Tech stack
C++20, header-only eval core under `SeQuant/core/eval/`, `TensorNetwork` (tensor_network/v2), range-v3, Catch2 tests under `tests/unit/`, CMake + Ninja. The dry-run backend (`sequant::eval::dryrun`) is the zero-data replay both SeQuant tests and MPQC drive through `cost_profile()`.

### Global Constraints (binding)
- No en-dashes (U+2013) anywhere in code, comments, commit messages, or this plan; use `--`.
- Format every touched file with `/opt/homebrew/opt/llvm/bin/clang-format --style=file -i <file>` before each commit.
- No `Co-Authored-By` trailers on commits.
- Any DP / witness / cost_profile run uses the RELEASE build in `cmake-build-release`; cap `ninja -j 6`.
- Acceptance is byte-identical: EMPTY router => `cost_profile()` equals the baseline on `peak_bytes`, `dryrun_flops`, `avoidable_flops`, `dryrun_n_ops` on `[dryrun-occ-veto]` + `[dryrun-extmode-avoidable]`, and the full `[eval]`/`[cache_manager]` suites stay green.
- OUT OF SCOPE (do not touch): `BatchPolicy::is_batchable_contracted_index`, `BatchPolicy::is_batchable_index()`.

### File Structure

Created:
- `SeQuant/core/eval/occurrence_key.hpp` -- `occurrence_key(node, batched)` (sub-expression -> `TensorNetwork` -> `canonicalize_slots`) plus the local `RouterKeyHash`/`RouterKeyEqual` functors. Header-only.
- `SeQuant/core/eval/placement_router.hpp` -- `PlacementRouter` (keyed by `SlotCanonicalizationMetadata`), `HomeTarget`, `home_depth()` resolver.
- `tests/unit/test_occurrence_key.cpp` and `tests/unit/test_placement_router.cpp` -- registered in `tests/unit/CMakeLists.txt`.

Modified:
- `SeQuant/core/eval/cache_manager.hpp` -- non-owning `PlacementRouter const*` field + setter + chain-walking getter (default null, inherited from `parent_`); `access_at_hops(key, hops)`.
- `SeQuant/core/eval/eval.hpp` -- Enter-stage read seam (714-730) and `place_at_this_level` store seam (1730-1739).
- `SeQuant/core/eval/backends/dryrun/cost_profile.hpp` -- optional caller-supplied `PlacementRouter const*` (default null) attached to the replay cache.
- `tests/unit/test_eval_dryrun.cpp` -- empty-router byte-identity guardrail + the positive relocation demonstration.

---

### Code recon confirmed (with line numbers)

Sub-expression -> TensorNetwork reconstruction (the load-bearing find for the key):
- `to_expr(node)` rebuilds the symbolic `ExprPtr` of an eval node's subtree (`eval_expr.hpp:637-670`); a leaf returns `evxpr.expr()`, an internal node returns a nested `Product`/`Sum`.
- The cleaner, `build_subnet_metadata`-matching path is to gather the subtree's TENSOR leaves and build the network from them: a leaf's tensor is `node->as_tensor()` (`eval_expr.hpp:210`, guarded by `node->is_tensor()` at `:160`). `TensorNetwork` accepts an `ExprPtrRange` of tensors (`v2.hpp:223`), exactly as `build_subnet_metadata` does (`single_term_detail.hpp:592-597`: `ts_expr` of cloned `Tensor`s -> `TensorNetwork{ts_expr}`).
- `TensorNetwork::canonicalize_slots(cardinal_tensor_labels, &named_indices)` returns a `SlotCanonicalizationMetadata` (`v2.hpp:258-318`); its `hash_value()` (`v2.hpp:282`) and `graph->cmp` (bliss) are the hash/equality, wrapped by `SubNetHash`/`SubNetEqual` (`single_term_detail.hpp:413-426`). Named-index coloring: passing ONLY the batched indices makes them position-meaningful (colored by space+label) and renames every other index (colored by space alone), so non-batched external/contracted slots are ignored and symmetry is intrinsic to the colored graph. Cardinal labels: `TensorCanonicalizer::cardinal_tensor_labels()` (`single_term_detail.hpp:599`).
- NON-trivial note: an eval-node subtree can contain scalar (`Variable`) leaves; those are skipped (only `is_tensor()` leaves join the network). Proto/composite leaves (`a1<i1,i2>`) reconstruct with their proto indices intact through `as_tensor()`/`clone()`, so PNO composite coloring is preserved (made an explicit test).

Scope / depth (UNCHANGED, reused):
- Scope == a `CacheManager` in the `parent_` chain (`cache_manager.hpp:197`); run/term cache is the chain root. Each realized loop pushes one outermost-first `batch_context` entry (`cache_manager.hpp:83-84`) and creates a scratch with `parent_` set (`eval.hpp:1723`). Depth == `batch_context().size()`.
- Read: `access_at` returns `hops` (links crossed) (`cache_manager.hpp:365-371`); `slice_to_use` slices the innermost `hops` loops filtered by `index_position(nd,axis)` (`eval.hpp:607-622`); `home_depth = use_depth - hops`.
- Store: `place_at_this_level` picks `rl` = deepest `ectx` index whose mode is in `sliced_modes() UNION contracted_modes()` (`in_union`, `eval.hpp:1724-1729`; `rl` walk, `1730-1739`), `rl == -1` => chain root; then walks up to that level and `target->store(d, ...)` (`1763-1775`). Collect/hoist gate (`residency_all_outer`, `has_demoted_external`, `eval.hpp:1671-1706`) is UNCHANGED in Phase 2.

Batched-set source for the key, at both sites (the consistency point):
- Read (Enter): `cache.batch_context()` modes (`ctx[i].first`), proto-expanded, filtered to those on `f.node`'s slots.
- Store (`place_at_this_level`): `ectx = parent_cache.batch_context()` modes (`eval.hpp:1665`), proto-expanded, filtered to those on `d`'s slots.
- Consistency argument (why the two agree for one occurrence): a batched mode `m` on the node is a residency mode, so the node's home is at-or-below `m`'s loop; therefore `m` encloses the node at ANY read-or-store scope of that occurrence, and both `batch_context`s contain it. A batched context loop NOT on the node's slots is filtered out on both sides. So `batched(node) = ctx-modes INTERSECT node-slot-modes` is site-invariant. Made an explicit test (iv).

---

### Task T1 -- the occurrence-key component + its tests

Interfaces (`occurrence_key.hpp`):
Produces:
```cpp
namespace sequant::eval {

/// Batched modes (proto-expanded) that actually appear on `node`'s own slots,
/// drawn from the ambient batch-loop modes `ctx_modes`. This is the site-
/// invariant named-index set (see the consistency argument): identical whether
/// computed at the read (cache.batch_context()) or store (ectx) site.
template <typename Node>
tensor_network::NamedIndexSet in_scope_batched_on_node(
    Node const& node, container::svector<Index> const& ctx_modes);

/// The router occurrence key: canonicalize the node's sub-expression as a
/// colored TensorNetwork with the in-scope batched indices as named_indices.
/// Symmetry is intrinsic to the colored graph; non-batched external/contracted
/// slots are renamed (ignored). Mirrors build_subnet_metadata
/// (single_term_detail.hpp:592-599).
template <typename Node>
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
  return tn.canonicalize_slots(
      TensorCanonicalizer::cardinal_tensor_labels(), &named);
}

/// Reuse of the SubNetHash/SubNetEqual pattern (single_term_detail.hpp:413-426),
/// re-declared here so eval does not depend on core/optimize.
struct RouterKeyHash {
  std::size_t operator()(
      TensorNetwork::SlotCanonicalizationMetadata const& m) const noexcept {
    return m.hash_value();
  }
};
struct RouterKeyEqual {
  bool operator()(TensorNetwork::SlotCanonicalizationMetadata const& a,
                  TensorNetwork::SlotCanonicalizationMetadata const& b) const {
    return bliss::ConstGraphCmp::cmp(*a.graph, *b.graph) == 0;
  }
};
}  // namespace sequant::eval
```
Consumes: `to_expr`/`as_tensor`/`is_tensor` (`eval_expr.hpp:637,210,160`); `TensorNetwork{ExprPtrRange}` + `canonicalize_slots` (`v2.hpp:223,314`); `TensorCanonicalizer::cardinal_tensor_labels()`; `Index::proto_indices()` for the proto-expansion in `in_scope_batched_on_node`.

`in_scope_batched_on_node` body: for each `m` in `ctx_modes`, proto-expand (`m.has_proto_indices()` -> its protos, else `m`); include an expanded mode iff it equals one of `node->canon_indices()` (proto-expanded on the slot side, same as `lifetime_mask.hpp:92-101 slot_modes_of`).

TDD steps:
1. `test_occurrence_key.cpp` case (i) -- flat non-batched genericization: build eval nodes for `A{i1,i2,i3}` and `A{i1,i2,i4}` with batched `{i1,i2}`; `RouterKeyEqual{}(occurrence_key(a,{i1,i2}), occurrence_key(b,{i1,i2}))` is TRUE (i3/i4 are renamed). Run -> fails (header absent).
2. Add `occurrence_key.hpp`. Run (i) -> passes.
3. Case (ii) -- gC proto-symmetry collapse: `gC{i1,i2,mu1,K1,a1<i1,i2>}` vs `gC{i2,i3,mu1,K1,a1<i2,i3>}` under batched `{i2}`; assert keys EQUAL (the i1<->i2 symmetry moves the batched index to a canonical slot; composite `a1<...>` coloring preserved). Run -> passes (or reveals a proto-expansion gap; fix in `in_scope_batched_on_node`).
4. Case (iii) -- distinct batched axes stay distinct: same tensor, batched `{i1}` vs batched `{i2}` produce UNEQUAL keys. Case (iv) -- read-vs-store equality: for one hoisted node, build its store-site `ctx_modes` (ectx) and a deeper read-site `ctx_modes` (consumer batch_context that adds a non-node loop) and assert `occurrence_key` is EQUAL across the two. Run -> passes.
5. clang-format; commit "eval: batching-aware occurrence key via canonicalize_slots (Phase 2 T1)".

Boundary justification: the key is the one genuinely new algorithmic component and the one correctness hazard (read==store consistency, symmetry, proto), so it is isolated and pinned by four targeted tests before any router or runtime wiring exists.

---

### Task T2 -- PlacementRouter type + CacheManager router pointer + access_at_hops + home_depth

Interfaces (`placement_router.hpp`):
Produces:
```cpp
namespace sequant::eval {

/// Override home descriptor. `residency` = enclosing batch-loop modes the value
/// should be variant to at its overridden home; resolves to a depth via the
/// rl+1 walk against the live batch_context (empty => chain root). split_index
/// is ALWAYS 0 in Phase 2 (Phase-4 peak splits; carried so no key refactor).
struct HomeTarget {
  container::svector<Index> residency;
  std::size_t split_index = 0;
};

template <typename TreeNode>
class PlacementRouter {
 public:
  using BatchContext = typename CacheManager<TreeNode, false>::BatchContext;
  using Key = TensorNetwork::SlotCanonicalizationMetadata;

  [[nodiscard]] bool empty() const noexcept { return overrides_.empty(); }
  void set_override(Key key, HomeTarget home) {
    overrides_.insert_or_assign(std::move(key), std::move(home));
  }
  /// nullptr => no override for this occurrence (the Phase-2 default) => caller
  /// derives home as today.
  [[nodiscard]] HomeTarget const* route(Key const& key) const {
    auto it = overrides_.find(key);
    return it == overrides_.end() ? nullptr : &it->second;
  }
  /// home_depth = 1 + (deepest index i in ctx whose mode is in home.residency),
  /// or 0 if none. Static mirror of place_at_this_level's rl+1 (eval.hpp:1730).
  [[nodiscard]] std::size_t home_depth(HomeTarget const& home,
                                       BatchContext const& ctx) const {
    for (int i = static_cast<int>(ctx.size()) - 1; i >= 0; --i)
      for (auto const& ix : home.residency)
        if (ctx[i].first == ix) return static_cast<std::size_t>(i) + 1;
    return 0;
  }
 private:
  std::unordered_map<Key, HomeTarget, RouterKeyHash, RouterKeyEqual> overrides_;
};
}  // namespace sequant::eval
```
Produces (`cache_manager.hpp`):
```cpp
void set_placement_router(PlacementRouter<TreeNode> const* r) noexcept {
  placement_router_ = r;
}
[[nodiscard]] PlacementRouter<TreeNode> const* placement_router() const noexcept {
  return placement_router_ ? placement_router_
         : parent_        ? parent_->placement_router()
                          : nullptr;
}
/// Fetch `key` from EXACTLY `hops` scopes up: walk `hops` parent links, then one
/// local entry::access() there (same single lifetime decay as access_at). Null
/// if that scope does not currently hold the value.
[[nodiscard]] ResultPtr access_at_hops(key_type const& key,
                                       std::size_t hops) noexcept {
  CacheManager* c = this;
  for (std::size_t i = 0; i < hops && c; ++i) c = c->parent_;
  if (!c) return nullptr;
  if (auto found = c->cache_map_.find(key); found != c->cache_map_.end())
    return found->second.access();
  return nullptr;
}
```
Consumes: `RouterKeyHash`/`RouterKeyEqual` (T1); `CacheManager::batch_context()`/`parent_`/`entry::access()`.

TDD steps:
1. `test_placement_router.cpp`: `home_depth` unit tests -- `{residency={K}}` against ctx `[{J},{K},{L}]` -> 2; `{residency={}}` -> 0; residency mode absent from ctx -> 0. Run -> fails (types absent).
2. Add `placement_router.hpp` + the `cache_manager.hpp` additions. Run -> passes.
3. `access_at_hops` unit test in `test_cache_manager.cpp` using the 3-level fixture (`test_cache_manager.cpp:405`): assert `access_at_hops(X,k)` returns the same buffer and the same single lifetime decay as `access_at(X)` for the `k` where X lives; `k` past the chain end returns null.
4. clang-format; commit "eval: PlacementRouter keyed on SlotCanonicalizationMetadata + access_at_hops (Phase 2 T2)".

Boundary justification: pure container + resolver + cache accessor, no runtime path touched, testable against hand-built chains; kept separate so the risky Enter/store wiring lands alone in T3.

---

### Task T3 -- wire the READ + STORE consult (empty => byte-identical) + shadow-assert

READ seam (`eval.hpp:717-730`), default path byte-identical:
```cpp
if (f.checked) {
  bool routed = false;
  if (auto const* router = cache.placement_router(); router && !router->empty()) {
    auto const& ctx = cache.batch_context();
    container::svector<Index> ctx_modes;
    for (auto const& e : ctx) ctx_modes.push_back(e.first);
    auto const key = occurrence_key(f.node, ctx_modes);
    if (auto const* home = router->route(key)) {
      std::size_t const use_depth = ctx.size();
      std::size_t const hd = router->home_depth(*home, ctx);
      SEQUANT_ASSERT(hd <= use_depth);
      std::size_t const hops = use_depth - hd;
      if (ResultPtr ptr = cache.access_at_hops(f.node, hops); ptr) {
        if constexpr (detail::trace(EvalTrace))
          log::cache(f.node, cache, log::label(f.node));
        finalize(slice_to_use(apply_phase(f.node, ptr), f.node, hops));
        routed = true;
      }
    }
  }
  if (routed) break;
  if (auto m = cache.access_at(f.node); m.ptr) {   // default: unchanged
    if constexpr (detail::trace(EvalTrace))
      log::cache(f.node, cache, log::label(f.node));
    finalize(slice_to_use(apply_phase(f.node, m.ptr), f.node, m.hops));
    break;
  }
  f.store_after = cache.exists(f.node);
}
```
STORE seam (`place_at_this_level`, replacing the `rl` selection at `eval.hpp:1730-1739` for a collected target `d`). SCOPE BOUND: only the LEVEL of an already-collected target changes; the collect/hoist decision is untouched.
```cpp
// Derived level (today): deepest ectx index whose mode is in d's union.
int rl = -1;
for (int i = static_cast<int>(ectx.size()) - 1; i >= 0; --i)
  if (in_union(d, ectx[i].first)) { rl = i; break; }
// Phase-2 override: if the router has an entry for this occurrence, take its
// level instead of the derived rl. Empty router => never taken => byte-identical.
if (auto const* router = parent_cache.placement_router();
    router && !router->empty()) {
  container::svector<Index> ectx_modes;
  for (auto const& e : ectx) ectx_modes.push_back(e.first);
  auto const key = occurrence_key(d, ectx_modes);
  if (auto const* home = router->route(key))
    rl = static_cast<int>(router->home_depth(*home, ectx)) - 1;  // 0 => -1 root
}
```
`cost_profile.hpp`: add `PlacementRouter<EvalNodeDryRun> const* router = nullptr` to `cost_profile()` and `cache.set_placement_router(router)` after `build_dryrun_cache` (`cost_profile.hpp:372`); internal build path passes null.

TDD steps:
1. Byte-identity guardrail (`test_eval_dryrun.cpp`): run `cost_profile()` over both witnesses with a null router and `CHECK` `peak_bytes`/`dryrun_flops`/`avoidable_flops`/`dryrun_n_ops` equal the current baseline constants. Run -> fails (cost_profile param absent).
2. Add the READ + STORE seams and the `cost_profile` param. Run the guardrail + `[eval]` + `[cache_manager]` -> byte-identical (empty router => key never computed, only default branches execute). clang-format; commit "eval: router consult on the Enter read and place_at_this_level store, seeded empty (Phase 2 T3)".
3. Shadow-assert landing (default-off, dev flag `SEQUANT_ROUTER_SHADOW`): in the READ seam, when the router routes, `SEQUANT_ASSERT(ptr.get() == cache.access_at(f.node).ptr.get())` for a NO-OP override (home resolving to the value's actual current scope). Exercise with a `test_placement_router.cpp` fixture that injects a no-op override (residency naming the value's real home, known from the fixture) and asserts routed == default pointer-for-pointer. This proves the resolver + `access_at_hops` reproduce `access_at` before any real relocation. Commit.

Boundary justification: T3 is the only task that mutates hot runtime paths; keeping the empty-default byte-identity proof (step 2) ahead of any behavior change makes the seam bisectable, and the shadow-assert validates the machinery against ground truth before T4 relocates for real.

---

### Task T4 -- the positive REAL relocation test (store participates, so it is meaningful)

TDD steps:
1. Integration relocation (`test_eval_dryrun.cpp`): over one witness forest, pick an already-hoisted target value V (one `place_at_this_level` collects, e.g. a run-scope-admitted invariant, `cache_manager.hpp:675-689`). Compute `occurrence_key(V_node, ectx_modes)` for V's store context and build a `PlacementRouter` with one `set_override(key, HomeTarget{.residency = <a shallower set>})` relocating V's home. Pass the router to `cost_profile()`. Assert: (a) V is BUILT/stored at the new home (observed via the eval trace `log::cache` store event / the `[eval]` Store record at the overridden depth); (b) the consuming read finds V there and `slice_to_use` slices it correctly for the use (result stays numerically valid -- DryRun sizing is finite/consistent); (c) `peak_bytes` moves in the expected direction vs the empty-router run (a shallower home lengthens V's live range across more inner iterations, so co-resident `chain_residency` at the inner peak point rises -- `cost_profile.hpp:196-219,428-429`). Run -> fails, then passes once T3 wiring is correct.
2. Add a focused unit mirror in `test_placement_router.cpp` (full control, deterministic): 3-level chain, store V at outer AND mid, override to outer, assert `access_at_hops` fetches the outer copy and `slice_to_use` slices 2 loops vs the default 1 -- the read-side relocation isolated from the replay.
3. clang-format; commit "eval: prove a router override relocates a value's build+read home (Phase 2 T4)".

Boundary justification: T4 is the only intentional behavior change; because the STORE side now participates (T3), the override is genuinely satisfiable (V is materialized at the new home), so this is a real end-to-end relocation, not a fixture-only demonstration.

---

### The one real risk: read-site vs store-site occurrence-key must match exactly

An override is placed against `occurrence_key(node, ectx_modes)` at the store and looked up against `occurrence_key(f.node, ctx_modes)` at the read. If those keys diverge for the same occurrence, the store puts V at the new home but the read looks up the old key, misses the override, resolves the DEFAULT home, and reads the wrong scope (or rebuilds) -- a correctness bug, not just a perf one. The divergence sources are (a) the batched named-index set differing between sites, (b) symmetry/proto canonicalization instability, (c) the frame deep-copy altering annotations.

How the tests prove it: T1 case (iv) asserts `occurrence_key` equality across a store-site and a deeper read-site context for one node (guards (a)); cases (i)-(iii) guard symmetry/proto/distinctness (guard (b)); the deep-copy preserves the `EvalExpr` annotation used to build the network (`binary_node.hpp:340-352` copies `data_`), and T4 step 1 exercises the full store-then-read round trip on a real forest (guards (c)). The empty-router guardrail (T3 step 2) guarantees none of this can perturb today's numbers.

Performance note (not a Phase-2 blocker): `occurrence_key` runs a bliss canonicalization; the Enter/store seams compute it ONLY when the router is non-empty (`router && !router->empty()` gate first), so Phase 2 (empty) pays zero hot-path cost and only T4/Phase-3 pay it. Memoizing the key per node is a Phase-3 follow-up.

---

### Genuine ambiguity flagged (not invented)

The `named_indices` set is defined as `ctx-modes INTERSECT node-slot-modes`, and I argued it is site-invariant because a batched mode on a node encloses that node at every scope. That argument assumes the collect/hoist gate keeps a value's home at-or-below all of its carried batched loops (true for `residency_all_outer` targets, `eval.hpp:1671-1677`). For a `has_demoted_external` node -- NOT hoisted, built and read locally -- the store side never overrides it (Phase 2 store override applies only to collected targets, decision C), so read/store key consistency for demoted nodes is never exercised in Phase 2, and whether a Phase-3 override could ever target a demoted node (and thus need a stable key for a value with no single home) is not settled by the current design. I did not invent a rule for it; I flag that the read/store key-equality guarantee is proven only for hoisted targets, which is exactly Phase 2's override scope, and that demoted-node keying is a Phase-3 question to resolve before overrides can target them.

---

### Deferred to later phases (state, do not plan)
- 7d `home_scope` AUTO-computation and WHICH overrides to create (incl. `sliced UNION contracted` vs `sliced UNION demoted_external`) -- Phase 3/4.
- Rational `W(cell)` reuse measure (change #2); the O2 whole-forest greedy split pass (change #3).
- `split_index > 0` same-scope stores and the path-based read-time occurrence-id (root index + L/R path) that separates two occurrences sharing one occurrence-key in one scope -- Phase 4.
- Memoization of `occurrence_key` and Phase-3 override targeting of demoted (non-hoisted) nodes.

---

## Phase 3a -- `home_scope` seed PREDICTOR (non-regressing, scope decision (A))

### Goal
Build the static `home_scope` seed predictor: a per-node UNIFIED residency = the
cross-occurrence MEET of ALL batched modes (any `BatchModeType`, not External-only)
that live on a node's own result slots. Expose it as `home_scope(node)` returning
that residency mode-set. Do this as a NEW static computation that does NOT touch
runtime `place_at_this_level` / `stamp_lifetime_masks` / `sliced_modes_`, so runtime
placement stays today's heuristic (NO regression). Reconcile the CAT-2
`occurrence_key`/`ext_modes_of` inconsistency and fix its stale comment.

### Architecture
`stamp_lifetime_masks` (`SeQuant/core/eval/lifetime_mask.hpp`) is a top-down forest
walk that accumulates ancestor+self batched modes (`ext_modes_of`, External-only),
filters each node's accumulator to the modes that live on that node's own canonical
slots (`slot_modes_of`), and stores the cross-occurrence set-intersection (meet) into
`EvalExpr::sliced_modes_`. Phase 3a parameterizes that walk by a mode-selector and a
result-setter, adds a second entry point `stamp_seed_residency` that runs the SAME
walk with an all-batched-modes selector writing a NEW `EvalExpr::seed_residency_`
field, and leaves the External-only `sliced_modes_` path byte-unchanged. The seed
predictor produces the residency mode-set only; runtime depth resolution against a
batch context is the existing rl-walk, reused later by 3b/O2 (out of scope here).

### Tech stack
C++20 header-heavy templates (`meta::eval_node` / `meta::eval_node_range` concepts),
Catch2 unit tests (`tests/unit/`), CMake + Ninja. Proto-aware `Index` (composite
indices expand to `proto_indices()`).

### Global Constraints (binding)
- No en-dashes (U+2013) anywhere in source or comments.
- Format touched files with `/opt/homebrew/opt/llvm/bin/clang-format --style=file -i`.
- No `Co-Authored-By` trailers in commits.
- RELEASE build in `cmake-build-release`; cap ninja at `-j6`.
- NON-REGRESSING: runtime `sliced_modes_` / `place_at_this_level` output stays
  BYTE-UNCHANGED. Existing `[lifetime_mask]`, `[eval]`, `[cache_manager]`, `[dryrun]`,
  and the `[.][dryrun-occ-veto]` / `[.][dryrun-extmode-avoidable]` witness figures
  stay green.
- READ-ONLY seed: do NOT delete `contracted_modes` / `has_demoted_external`; do NOT
  switch the runtime selector. Those are the Phase-4 cutover.

### File Structure
- MODIFY `SeQuant/core/eval/lifetime_mask.hpp`
  - `detail::proto_expand_into` / `detail::slot_modes_of` free helpers (hoisted from
    the two in-file lambdas), reused by `occurrence_key.hpp`.
  - `detail::stamp_residency_impl(forest, modes_of, setter)` -- the parameterized walk.
  - `stamp_lifetime_masks` becomes a thin wrapper (`ext_modes_of` selector +
    `set_sliced_modes`) -- output byte-identical.
  - NEW `stamp_seed_residency` (`all_batched_modes_of` selector + `set_seed_residency`).
  - NEW `home_scope(node)` free accessor.
- MODIFY `SeQuant/core/eval/eval_expr.hpp`
  - NEW `seed_residency()` / `set_seed_residency()` accessors (next to `sliced_modes`).
  - NEW member `container::svector<Index> seed_residency_{};` (next to `sliced_modes_`).
- MODIFY `SeQuant/core/eval/occurrence_key.hpp`
  - `in_scope_batched_on_node` reuses `detail::proto_expand_into` / `slot_modes_of`;
    stale "mirrors `ext_modes_of` exactly" comment fixed to reference the unified
    all-batched-modes selector.
- MODIFY `tests/unit/test_lifetime_mask.cpp`
  - Seed-residency against-definition tests, External-only equivalence cross-check,
    Contracted-generalization test, and the byte-unchanged-runtime guardrail.
- MODIFY `tests/unit/test_occurrence_key.cpp`
  - A pin that `in_scope_batched_on_node` is kind-agnostic (a Contracted ambient mode
    on a node's slot is included), locking the reconciliation.

---

### T1 -- Parameterize the walk; add `stamp_seed_residency` + `seed_residency_`

#### Files
- MODIFY `SeQuant/core/eval/eval_expr.hpp`
  - accessors after `set_sliced_modes` (currently ~310-312)
  - member after `sliced_modes_` (currently line 396)
- MODIFY `SeQuant/core/eval/lifetime_mask.hpp`
  - refactor `stamp_lifetime_masks` body (currently 58-134); `ext_modes_of` 73-83,
    `slot_modes_of` 92-101, walk 107-125, stamp loop 128-134.

#### Interfaces
Consumes:
- `EvalExpr::batched_here() -> container::svector<std::pair<Index,BatchModeType>> const&`
  (eval_expr.hpp:283)
- `EvalExpr::canon_indices() -> index_vector const&` (eval_expr.hpp:248)
- `enum class BatchModeType { Contracted, External }` (fwd.hpp:21)
- `detail::lifetime_mask_intersect_in_place(svector<Index>&, svector<Index> const&)`
  (lifetime_mask.hpp:20)

Produces:
- `void EvalExpr::set_seed_residency(container::svector<Index>) noexcept`
- `container::svector<Index> const& EvalExpr::seed_residency() const noexcept`
- `template <meta::eval_node_range R> void sequant::stamp_seed_residency(R const&) noexcept`
- (internal) `template <meta::eval_node_range R, typename ModeSelector, typename Setter>
  void detail::stamp_residency_impl(R const&, ModeSelector&&, Setter&&) noexcept`

#### TDD steps
1. Failing test in `test_lifetime_mask.cpp`: a forest with only External stamps, run
   `stamp_seed_residency`, assert `seed_residency()` equals what `stamp_lifetime_masks`
   produces for `sliced_modes()` on the identical forest (External-only => the two
   selectors coincide). Also assert `seed_residency()` is non-empty for a
   sliced-everywhere node. Compile fails: no `stamp_seed_residency` / `seed_residency`.
2. `eval_expr.hpp`: add member and accessors, mirroring `sliced_modes_`:
   ```cpp
   [[nodiscard]] container::svector<Index> const& seed_residency() const noexcept {
     return seed_residency_;
   }
   void set_seed_residency(container::svector<Index> m) noexcept {
     seed_residency_ = std::move(m);
   }
   ...
   /// See \c seed_residency.
   container::svector<Index> seed_residency_{};
   ```
   Doc `seed_residency` as: the cross-occurrence meet of ALL batched modes (any
   `BatchModeType`) on this node's own result slots; the perfect-CSE `home_scope`
   seed; set by `stamp_seed_residency`; empty by default (OFF path); distinct from
   the External-only runtime `sliced_modes`.
3. `lifetime_mask.hpp`: hoist `proto_expand_into` + `slot_modes_of` to `detail` free
   functions; extract `stamp_residency_impl(forest, modes_of, setter)` carrying the
   existing walk verbatim (only `ext_modes_of` -> `modes_of`, and the const_cast set
   -> `setter`). Then:
   ```cpp
   template <meta::eval_node_range R>
   void stamp_lifetime_masks(R const& forest) noexcept {
     using Node = std::ranges::range_value_t<R>;
     using Data = typename Node::value_type;
     auto ext_modes_of = [](Node const& n) {
       container::svector<Index> v;
       for (auto const& [ix, kind] : n->batched_here())
         if (kind == BatchModeType::External) detail::proto_expand_into(v, ix);
       return v;
     };
     detail::stamp_residency_impl(
         forest, ext_modes_of,
         [](Node const* n, container::svector<Index> m) {
           const_cast<Data&>(**n).set_sliced_modes(std::move(m));
         });
   }

   template <meta::eval_node_range R>
   void stamp_seed_residency(R const& forest) noexcept {
     using Node = std::ranges::range_value_t<R>;
     using Data = typename Node::value_type;
     auto all_batched_modes_of = [](Node const& n) {
       container::svector<Index> v;
       for (auto const& [ix, kind] : n->batched_here())
         detail::proto_expand_into(v, ix);  // every kind, not just External
       return v;
     };
     detail::stamp_residency_impl(
         forest, all_batched_modes_of,
         [](Node const* n, container::svector<Index> m) {
           const_cast<Data&>(**n).set_seed_residency(std::move(m));
         });
   }
   ```
4. Run `[lifetime_mask]`; the new test plus all existing meet/off-path/proto/per-slot/
   veto cases pass (proving the `stamp_lifetime_masks` refactor is behavior-preserving).
5. clang-format touched files; commit "eval: parameterized residency walk + seed_residency".

---

### T2 -- `home_scope` accessor + CAT-2 `in_scope_batched_on_node` reconciliation

#### Files
- MODIFY `SeQuant/core/eval/lifetime_mask.hpp` (add `home_scope` free fn)
- MODIFY `SeQuant/core/eval/occurrence_key.hpp` (share helpers; fix comment) --
  `in_scope_batched_on_node` 41-67, stale comment 26-28, `#include` block 4-10.
- MODIFY `tests/unit/test_occurrence_key.cpp` (kind-agnostic pin)

#### Interfaces
Consumes:
- `detail::proto_expand_into` / `detail::slot_modes_of` (from T1)
- `EvalExpr::seed_residency()` (from T1)

Produces:
- `template <meta::eval_node Node>
   container::svector<Index> const& home_scope(Node const& n)` -> `n->seed_residency()`
- reconciled `in_scope_batched_on_node` (unchanged signature/return
  `TensorNetwork::NamedIndexSet`) built on the shared helpers.

#### TDD steps
1. Failing pin in `test_occurrence_key.cpp`: a node with a Contracted-kind mode
   sitting on its own slot; call `in_scope_batched_on_node(node, ctx)` with that mode
   in `ctx`; assert it IS returned (kind-agnostic). This ALREADY passes for the
   current impl (it never inspected kind) -- it locks the reconciled semantics against
   regression and documents the fixed comment's claim. Add a `home_scope` smoke test
   in `test_lifetime_mask.cpp`: after `stamp_seed_residency`, `home_scope(node)` ==
   `node->seed_residency()`. Compile fails: no `home_scope`.
2. `lifetime_mask.hpp`: add
   ```cpp
   /// The perfect-CSE seed home residency of \p n: the unified all-batched-modes
   /// meet on its own result slots (see \c stamp_seed_residency). Phase 3a returns
   /// the residency mode-set; runtime depth resolution against a batch context is
   /// the existing rl-walk, reused by 3b/O2.
   template <meta::eval_node Node>
   container::svector<Index> const& home_scope(Node const& n) noexcept {
     return n->seed_residency();
   }
   ```
3. `occurrence_key.hpp`: `#include <SeQuant/core/eval/lifetime_mask.hpp>`; rewrite the
   two inline proto-expansion loops (46-51 ctx side, 55-60 slot side) to
   `detail::proto_expand_into` / `detail::slot_modes_of`; keep the `NamedIndexSet`
   intersection (62-66). Fix the comment (26-28): it must state this filters ALL
   in-scope batched modes (any `BatchModeType`) to `node`'s own slots, matching the
   unified `all_batched_modes_of` selector / `slot_modes_of` in `lifetime_mask.hpp`
   -- NOT the External-only `ext_modes_of` (which it never mirrored: `ctx_modes` is
   already kind-agnostic).
4. Run `[occurrence_key]` + `[lifetime_mask]`; green.
5. clang-format; commit "eval: home_scope accessor + reconcile in_scope_batched_on_node".

---

### T3 -- Validation: against-definition, cross-check, byte-unchanged guardrail

#### Files
- MODIFY `tests/unit/test_lifetime_mask.cpp` (reuses `head` / `leaf` / `inode` /
  `stamp_ext` / `stamp_ext_pair` helpers, lines 30-80; add a `stamp_con` /
  `stamp_con_pair` shim for `BatchModeType::Contracted`).

#### Interfaces
Consumes: `stamp_lifetime_masks`, `stamp_seed_residency`, `home_scope`,
`EvalExpr::{sliced_modes,seed_residency,mask_all_full}`.
Produces: `[lifetime_mask][seed]` test cases.

#### TDD steps
1. AGAINST-DEFINITION (constructed): reuse the meet / per-slot / proto forests from the
   existing `[lifetime_mask]` cases but drive `stamp_seed_residency`; assert
   `seed_residency` equals the hand-computed deepest all-batched-modes meet on result
   slots (e.g. pair-sliced-everywhere -> `{i,j}`; full-in-one-occ -> empty;
   loop-invariant sibling -> empty).
2. EXTERNAL-ONLY EQUIVALENCE CROSS-CHECK: on any External-only forest,
   `stamp_seed_residency` and `stamp_lifetime_masks` must agree node-for-node
   (`seed_residency() == sliced_modes()`), because every batched mode is External so
   the two selectors coincide. This is the "sanity cross-check vs today's sliced_modes
   where they agree" (with `contracted_modes` empty, `sliced ∪ contracted == sliced`).
3. CONTRACTED GENERALIZATION: a forest with a Contracted mode `c` on a node's own
   result slot in EVERY occurrence. Assert `seed_residency()` CONTAINS `c` while
   `sliced_modes()` (External-only) does NOT -- proving the dropped-filter
   generalization. (Use `stamp_con_pair` shim mirroring `stamp_ext_pair`.)
4. BYTE-UNCHANGED RUNTIME GUARDRAIL (the one real risk): on a mixed
   External+Contracted forest, call `stamp_lifetime_masks`, snapshot every node's
   `sliced_modes()`, THEN call `stamp_seed_residency`, and assert every node's
   `sliced_modes()` is byte-identical to its snapshot (seed writes only
   `seed_residency_`). Also assert the reverse order leaves `sliced_modes()` equal to
   a fresh `stamp_lifetime_masks`-only run.
5. Run `[lifetime_mask]`; then the full `[eval] [cache_manager] [dryrun]` suites in
   `cmake-build-release` (`ninja -j6`) plus the `[.][dryrun-occ-veto]` /
   `[.][dryrun-extmode-avoidable]` witnesses -- all figures unchanged.
6. clang-format; commit "test: seed_residency against-definition + byte-unchanged guardrail".

---

### The one real risk
The unified all-modes meet must NOT perturb runtime `sliced_modes_`. Mitigation is
structural + proven: `stamp_seed_residency` writes ONLY `seed_residency_`, and
`stamp_lifetime_masks` is refactored to a thin wrapper over the SAME
`stamp_residency_impl` with the UNCHANGED `ext_modes_of` selector and `set_sliced_modes`
setter -- so its output is behavior-preserving by construction. T3 step 4 turns that
into an executable guarantee (snapshot sliced_modes, run the seed pass, assert
byte-identical), and the existing `[cache_manager]` batch-variant veto (which reads
`mask_all_full()`, cache_manager.hpp:672) plus the witness figures stay green.

### Deferred to Phase 4 cutover (OUT OF SCOPE here)
- CAT-1 DELETE: `contracted_modes` end-to-end (field/accessors eval_expr.hpp:314-334,
  396-399; `NodeBatchAnnotation`; cost-model emission; `in_union` union) and the
  `has_demoted_external` veto (eval.hpp).
- CAT-2 unify: switch the RUNTIME selector to all-batched-modes (drop the `External`
  filter in `ext_modes_of`) -- a one-line change once T1 has parameterized the walk;
  collapses `sliced ∪ contracted` back to `sliced`.
- Wiring `home_scope` into the static peak profile (3b) and O2's `ΔPeak`, behind the
  seed+O2 coupling flag so the peak-maximal seed never ships regressed standalone.

## Phase 3b -- static peak profile (O3a; non-regressing, scope decision (A))

### Goal
Build the STATIC peak-profile analysis consuming the Phase-3a `home_scope` seed:
per-cell home-relative footprint x `[first-use, last-use]` liveness interval, swept to
`peak_bytes` plus the binding peak point and its live set (the spill candidates O2 needs
in Phase 4). A NEW additive analysis with ZERO production callers, validated by tests.
Runtime stays byte-identical: `place_at_this_level`, `stamp_lifetime_masks`, and the
Phase-1 runtime `cost_profile` path are untouched. See design doc section 9.

### Architecture
A seed cell is `(value, home_scope)` -- perfect CSE, one cell per value at its meet home,
`split_index == 0` (splits are Phase 4). Footprint reuses the EXISTING moment-aware sizer
`dryrun::CostModel::memsize(idxset, ExtentOverrides)` (`cost_model_object.hpp:120`,
`ExtentOverrides = container::map<Index,std::size_t>` at :84) with a home-relative override
table (block-extent for carried modes whose loop encloses the home, full otherwise). Home
DEPTH is the rl-walk over the enclosing-loop `BatchContext` keyed on `home_scope(node)`,
mirroring `eval.hpp:1776-1782` but using the unified seed residency instead of
`sliced ∪ contracted`. Liveness intervals come from the linearized static schedule
(one representative iteration per loop) + the Phase-2 router use-sites. An event-based
sweep gives the peak; a deliberately-different per-point live-set replay is the oracle.

### Tech stack
C++20 header (`meta::eval_node` / `meta::eval_node_range` concepts), the dry-run backend
(`backends/dryrun/`), Catch2 (`tests/unit/`), CMake + Ninja. Proto-aware `Index`.

### Global Constraints (binding)
- No en-dashes (U+2013) anywhere in source or comments.
- Format touched files with `/opt/homebrew/opt/llvm/bin/clang-format --style=file -i`.
- No `Co-Authored-By` trailers in commits.
- RELEASE build in `cmake-build-release`; cap ninja at `-j6`.
- NON-REGRESSING: ZERO production callers; runtime `place_at_this_level` /
  `stamp_lifetime_masks` / the runtime `cost_profile` path stay BYTE-UNCHANGED. Existing
  `[eval]`, `[cache_manager]`, `[dryrun]`, `[lifetime_mask]`, `[occurrence_key]`,
  `[placement_router]` suites + the `[.][dryrun-occ-veto]` / `[.][dryrun-extmode-avoidable]`
  witness figures stay green.
- READ-ONLY vs 3a/Phase-4: do NOT delete `contracted_modes` / `has_demoted_external`, do
  NOT switch the runtime selector, do NOT wire `home_scope`/`PeakProfile` into runtime.

### File Structure
- NEW `SeQuant/core/eval/peak_profile.hpp`
  - `struct PeakProfile { double peak_bytes; std::size_t binding_point;
    container::svector<std::size_t> live_at_binding; };` (cells + points identified by
    stable index into the linearization).
  - `detail::cell_footprint(...)` -- home-relative `ExtentOverrides` + `CostModel::memsize`.
  - `detail::home_depth_of(home_scope_modes, ectx)` -- the rl-walk keyed on `home_scope`.
  - `detail::linearize(forest) -> Schedule` -- static point order + per-cell
    `[first_use,last_use]` from router use-sites.
  - `peak_profile_sweep(Schedule) -> PeakProfile` -- the interval-event sweep.
  - `peak_profile_replay(Schedule) -> double` -- the per-point live-set oracle.
- NEW `tests/unit/test_peak_profile.cpp` (add to `tests/unit/CMakeLists.txt`).

---

### T1 -- footprint sizer + home-depth resolver

#### Files
- CREATE `SeQuant/core/eval/peak_profile.hpp` (the two sizing primitives + their `detail`).
- CREATE `tests/unit/test_peak_profile.cpp`; add it to `tests/unit/CMakeLists.txt`.

#### Interfaces
Consumes:
- `dryrun::CostModel::memsize(container::svector<Index> const&, ExtentOverrides const&) const`
  (`backends/dryrun/cost_model_object.hpp:120`); `using ExtentOverrides =
  container::map<Index,std::size_t>` (:84); `dryrun::SizeRegime` ctor of `CostModel` (:104).
- `home_scope(node) -> container::svector<Index> const&` (Phase 3a, `lifetime_mask.hpp`).
- `EvalExpr::batched_here() -> svector<pair<Index,BatchModeType>> const&` (eval_expr.hpp:283).
- `CacheManager::BatchContext` (`cache_manager.hpp:95`) -- the `ectx` shape
  (`svector<pair<Index, ...>>`; index is `.first`); mirror the rl-walk `eval.hpp:1776-1782`.

Produces:
- `std::size_t detail::cell_footprint(container::svector<Index> const& carried,
   container::svector<Index> const& home_modes,
   dryrun::CostModel const& cm, /* block extents */ ExtentFn const& block_of)`
  (review fix, Phase 3b: NO `ectx` param -- the footprint depends only on the
  cross-occurrence meet `home_modes` intersected with `carried`, so it is
  occurrence-independent by construction; see the T1 step-4 body below).
- `int detail::home_depth_of(container::svector<Index> const& home_modes,
   BatchCtx const& ectx)` -- deepest `ectx` level whose mode is in `home_modes`, else -1.

#### TDD steps
1. Failing test (`[peak_profile]`) for `home_depth_of`: build an `ectx` of loops
   `[a, b, c]` (innermost last) and assert:
   - `home_modes = {c}` -> depth `2`; `{a}` -> `0`; `{}` -> `-1`;
   - `{b, x}` (x not in ectx) -> `1` (deepest PRESENT residency mode);
   - proving it matches the `in_union`/rl walk keyed on `home_modes`. Compile fails: no header.
2. Write `home_depth_of` (transcribe `eval.hpp:1776-1782`, predicate = `m in home_modes`):
   ```cpp
   inline int home_depth_of(container::svector<Index> const& home_modes,
                            BatchCtx const& ectx) noexcept {
     for (int i = static_cast<int>(ectx.size()) - 1; i >= 0; --i)
       if (std::find(home_modes.begin(), home_modes.end(), ectx[i].first) !=
           home_modes.end())
         return i;
     return -1;
   }
   ```
3. Failing test for `cell_footprint`: a cell carrying modes `{p, q}` homed at scope with
   `ectx = [p]` (so `p`'s loop ENCLOSES the home -> `p` sized BLOCK, `q` unbatched -> FULL).
   With `full_extent(p)=full_extent(q)=10`, `block_extent(p)=2`, assert
   `cell_footprint == memsize({p,q}, {p->2})` == `2*10` (vs the full `10*10` when homed
   ABOVE `p`'s loop). Drive `CostModel` with a hand-built `SizeRegime` giving those
   extents (mirror the `SizeRegime` setup already used in `test_eval_dryrun.cpp` /
   `[dryrun][cost_profile]` tests -- read one for the exact constructor call).
4. Write `cell_footprint` (review fix, Phase 3b): build `ExtentOverrides` = `{ m ->
   block_of(m) : m in home_modes }` (block iff `m` is IN the meet `home_modes` -- no
   `ectx`/depth scan), then `return cm.memsize(carried, overrides);`. Modes not in
   `home_modes` get no override (full). NOTE: at the `linearize` call site, the
   `home_modes` actually passed in is the meet MINUS any mode a value ever realizes as
   its OWN loop (see T2 step 4) -- a value's own loop slices its operands, never its own
   result, so that self-contributed meet member must not be block-sized here.
5. Run `[peak_profile]`; both pass. Keep `[dryrun]` green (only a new header added).
6. clang-format; commit "eval: peak-profile footprint sizer + home-depth resolver (Phase 3b T1)".

---

### T2 -- linearization, liveness intervals, and the sweep -> `PeakProfile`

#### Files
- MODIFY `SeQuant/core/eval/peak_profile.hpp` (add `Schedule`, `linearize`,
  `PeakProfile`, `peak_profile_sweep`).
- MODIFY `tests/unit/test_peak_profile.cpp`.

#### Interfaces
Consumes (from T1): `cell_footprint`, `home_depth_of`. Use-sites are STRUCTURAL (NOT the
router -- the Phase-2 `PlacementRouter` is the seed's EMPTY override seam): a value's
consumers are its PARENT nodes in the hash-merged forest, read from the tree. Cells are
grouped by the value hash the `CacheManager` keys on (`EvalExpr::hash_value()` /
`cache_manager.hpp` key).
Produces:
- `struct Cell { std::size_t value_id; int home_depth; std::size_t footprint;
   std::size_t first_use, last_use; };`
- `struct Schedule { container::svector<Cell> cells; std::size_t num_points; };`
- `template <meta::eval_node_range R> Schedule linearize(R const& forest,
   dryrun::CostModel const&, /* block extents */ ExtentFn const&);`
- `struct PeakProfile { double peak_bytes; std::size_t binding_point;
   container::svector<std::size_t> live_at_binding; };`
- `PeakProfile peak_profile_sweep(Schedule const&);`

#### TDD steps
1. Failing test for `peak_profile_sweep` on a HAND-BUILT `Schedule` (bypass `linearize`;
   construct `Cell`s directly so the sweep is tested in isolation): three cells
   `A=[0,3] fp=100`, `B=[1,2] fp=40`, `C=[2,4] fp=10`; overlaps: point 2 has `{A,B,C}` =
   150 (max). Assert `peak_bytes==150`, `binding_point==2`, `live_at_binding=={A,B,C ids}`.
   Add a tie-break case (two points equal peak -> lowest point index wins) to pin
   determinism. Compile fails: no `peak_profile_sweep`.
2. Write `peak_profile_sweep`: emit `(+fp at first_use)` / `(-fp at last_use+1)` events,
   sort by point (deterministic), running-sum; track argmax `(peak, point)` with lowest-point
   tie-break; then a second pass (or retained live-set) collects `live_at_binding` = cells
   with `first_use <= binding_point <= last_use`.
   ```cpp
   inline PeakProfile peak_profile_sweep(Schedule const& s) {
     // events[p] = net footprint delta entering point p
     container::svector<double> delta(s.num_points + 1, 0.0);
     for (auto const& c : s.cells) {
       delta[c.first_use] += double(c.footprint);
       delta[c.last_use + 1] -= double(c.footprint);
     }
     double run = 0, peak = 0; std::size_t arg = 0;
     for (std::size_t p = 0; p < s.num_points; ++p) {
       run += delta[p];
       if (run > peak) { peak = run; arg = p; }  // strict => lowest-point tie-break
     }
     PeakProfile out{peak, arg, {}};
     for (std::size_t i = 0; i < s.cells.size(); ++i)
       if (s.cells[i].first_use <= arg && arg <= s.cells[i].last_use)
         out.live_at_binding.push_back(i);
     return out;
   }
   ```
3. Failing test for `linearize` on a SMALL constructed forest (reuse the eval-forest
   builders from `test_eval_dryrun.cpp` / `test_lifetime_mask.cpp`): a value shared by two
   consumers at different schedule points -> ONE cell whose `[first_use,last_use]` spans
   from its production to the LATER consumer; a value homed above a loop -> interval covers
   the loop subtree span. Assert the produced `Cell`s' intervals and footprints.
4. Write `linearize`: run `stamp_seed_residency(forest)`; POST-ORDER DFS the forest
   (children before parent) assigning each node a static point index (`num_points`).
   Maintain, during the DFS, the accumulated ancestor-loop stack as the per-node `ectx` --
   push each node's `batched_here()` loops on descent, pop on ascent -- since statically the
   enclosing `BatchContext` of a node IS its ordered ancestor loops (the runtime
   `parent_cache.batch_context()`). Group nodes by their VALUE hash
   (`EvalExpr::hash_value()` -- the same key `CacheManager` uses for CSE); each group is ONE
   `Cell`:
   - `first_use` = min static point over the group's nodes (its single production under
     perfect CSE);
   - `last_use` = max static point over every node that CONSUMES the value -- i.e. the
     PARENT of any node in the group (structural use-sites; NOT the empty router). A value
     homed above a loop needs no separate subtree-extension rule: its last consumer inside
     the loop already sets `last_use` past the loop body.
   - `home_depth = home_depth_of(home_scope(node), ectx-at-node)` (informational; still
     `ectx`-based and thus not itself occurrence-invariant -- only the footprint below is).
   - `footprint = cell_footprint(carried, home_modes, cm, block_of)` (review fix, Phase
     3b: NO `ectx` arg). `carried` = the node's `canon_indices()`. `home_modes` = `
     home_scope(node)` MINUS the cross-occurrence UNION of every occurrence's OWN
     realized-here modes (`batched_here()`, proto-expanded) for that same canonical
     value: `home_scope` folds a node's own loop into its own meet whenever that loop's
     mode is also one of the node's own carried slots, but a value's own loop slices its
     OPERANDS on the way down, never its own (loop-result) value, so that self-
     contribution must be excluded before block-sizing. Both the meet and the exclusion
     union are computed over ALL occurrences, so the resulting `home_modes` -- and hence
     the footprint -- is identical regardless of forest order.
5. Run `[peak_profile]`; all green. `[eval]`/`[dryrun]` green.
6. clang-format; commit "eval: peak-profile linearize + interval sweep -> PeakProfile (Phase 3b T2)".

---

### T3 -- the replay oracle + validation (oracle equality all forests; Phase-1 anchor)

#### Files
- MODIFY `SeQuant/core/eval/peak_profile.hpp` (add `peak_profile_replay`).
- MODIFY `tests/unit/test_peak_profile.cpp`.

#### Interfaces
Consumes: `Schedule`, `peak_profile_sweep` (T2); `dryrun::cost_profile(...) ->
CostProfile{ .peak_bytes }` (`backends/dryrun/cost_profile.hpp:189,221`) for the anchor.
Produces: `double peak_profile_replay(Schedule const&);` (per-point live-set sum).

#### TDD steps
1. Failing test: on the SAME hand-built `Schedule` from T2 step 1, assert
   `peak_profile_replay(s) == peak_profile_sweep(s).peak_bytes` (== 150). Compile fails.
2. Write `peak_profile_replay` -- a DELIBERATELY DIFFERENT algorithm (explicit per-point
   live-set sum, NOT interval events), so agreement is a real cross-check:
   ```cpp
   inline double peak_profile_replay(Schedule const& s) {
     double peak = 0;
     for (std::size_t p = 0; p < s.num_points; ++p) {
       double sum = 0;
       for (auto const& c : s.cells)
         if (c.first_use <= p && p <= c.last_use) sum += double(c.footprint);
       peak = std::max(peak, sum);
     }
     return peak;
   }
   ```
3. ORACLE EQUALITY (all forests): drive `linearize` on several constructed forests
   INCLUDING a demoted case (a value carried FULL above a loop it slices in another
   occurrence -> meet-empty `home_scope` -> homed at root, spanning the whole subtree) and
   assert `peak_profile_replay(sched) == peak_profile_sweep(sched).peak_bytes` EXACTLY for
   each. This is O3b; it must hold on demoted forests too.
4. PHASE-1 ANCHOR (non-demoted only): on an External-only / no-demotion forest where the
   seed placement coincides with today's heuristic, build the same `dryrun` evaluation the
   `[dryrun][cost_profile]` tests use, call `cost_profile(...)`, and assert
   `peak_profile_sweep(linearize(forest,...)).peak_bytes == profile.peak_bytes`. Document in
   the test WHY this is gated to non-demoted forests (design doc section 9.6: the seed and
   the heuristic diverge on demoted values, so the anchor is only valid where they agree).
   If no existing small non-demoted `cost_profile` fixture matches, reuse the smallest
   `[dryrun][cost_profile]` forest that has no demotion and assert equality there.
5. Run `[peak_profile]` + `[dryrun]`; then the `[.][dryrun-occ-veto]` /
   `[.][dryrun-extmode-avoidable]` witnesses -- their figures UNCHANGED (3b adds no runtime
   path; this confirms non-regression).
6. clang-format; commit "test: peak-profile replay oracle + Phase-1 anchor (Phase 3b T3)".

---

### The one real risk
The anchor (T3 step 4) equates a STATIC sweep over the SEED placement with the Phase-1
RUNTIME `peak_bytes` measured under today's HEURISTIC placement. That equality holds ONLY
where the two placements coincide (no demotion). Mitigation: gate the anchor to a
provably non-demoted forest (External-only, `home_scope` == today's `sliced_modes`), and
lean on the O3b oracle (T3 step 3) -- which shares the seed placement with the sweep -- for
the demoted cases. If the anchor forest turns out to have any demotion, the sweep and
`cost_profile` will legitimately differ and the test must NOT be "fixed" by loosening the
sizing; pick a smaller non-demoted forest instead.

### Deferred to Phase 4 (OUT OF SCOPE here)
- Wiring `PeakProfile` into runtime placement + O2's `ΔPeak` (shrink/un-hoist/split by
  `ΔPeak/ΔRecompute`), coupled with the perfect-CSE seed behind the flag.
- Incremental sweep updates (O3a incremental) -- O2 needs re-cost over one changed cell's
  span; the static full sweep here is the baseline it will make incremental.
