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
