# Lifetime-mask placement: retire the pre-slice-on-use placement apparatus

- **Date:** 2026-07-29
- **Status:** Design (approved shape; pending spec review)
- **Branch:** `evaleev/feature/suppress-heuristic-fallback` (SeQuant)
- **Depends on:** slice-on-use (phase 1), landed `af79d2707`. See
  `doc/dev/specs/2026-07-27-external-loop-scope-semantics-design.md` (esp. section 10)
  and `doc/dev/plans/2026-07-28-slice-on-use-phase1-plan.md`.
- **Scope:** SeQuant batched evaluator (`core/eval/eval.hpp`,
  `core/eval/cache_manager.hpp`) + cost-model emit (`core/eval/cost_model.hpp`).
  Runtime + emit only; no DP cost/sizing change.

## 0. Relationship to phase 1

Phase 1 (slice-on-use) made the batched evaluator *correct* regardless of where a
cached node lives: a node fetched from an ancestor scope is sliced to exactly the
batch loops the fetch crossed (`access_at` returns `{ptr, hops}`; the Enter stage
slices the `hops` innermost modes). With correctness decoupled from placement, the
C60 dry-run gate went from a 5860.877 GB replay to 443.55 GB, matching the DP model
(summand 12: 443.55 ~= 445.46).

This spec is the first phase-2 step. It does **not** chase the remaining 443 GB (that
is the DP's honest cheapest schedule with the current axes -- a separate cost-model /
new-lever question). It **retires the placement apparatus that predates slice-on-use**
-- the scalar per-occurrence `batch_scope_level`, the `consumed_upscope` bit, the
run-scope veto part (b), and `hoist_invariants` -- and replaces it with a single
per-slot lifetime mask plus per-level placement. The payoff is correctness of
*placement* (consistent across occurrences, so cross-term CSE becomes operational) and
a large reduction in accreted special-case machinery.

## 1. Problem

Placement today is three cooperating mechanisms in the batched evaluator:

1. **Run-scope membership** (`cache_manager()` builder, `cache_manager.hpp`): a node
   is registered in the run/term (root) cache iff it is an NP repeat or a P-frontier
   node, is not vetoed, and is under the footprint gate. The **veto** excludes
   (a) a node whose `batched_here()` carries a *contracted* mode free on its result
   (`sliced_batch_axis`), and (b) any node with `batch_scope_level() >= 0`
   (batch-variant: placed inside an enclosing batch loop).
2. **`hoist_invariants`** (`eval.hpp`, fired per batch loop): collects descendants
   whose scalar `batch_scope_level` is *strictly outer* to the firing loop, builds
   each once on a fresh cache (so it is sliced only to *its own* enclosing modes, not
   the deeper loop's), and stores it at its level by walking up the scope chain.
   `consumed_upscope` is a patch bit that stops an inner External carrier (which is
   loop-local, not invariant) from being mis-collected as an invariant.
3. **Normal per-batch store** on the current scratch (dropped at reset).

Two problems.

**(P1) The scalar scope annotation is per-occurrence, so placement is inconsistent
across occurrences of the same canonical node -- and the runtime cache is keyed by
canonical identity.** The IR survey
(`.superpowers/sdd/lifetime-mask-ir-survey.md`) shows the `s*C` overlap-times-PNO
intermediate (canonical hash `1989507463377952644`) appearing 52 times, at
`scope=-1` in ~28 occurrences and `scope>=0` (re-sliced) in ~24. The cache takes
whichever occurrence stores first; the other 24 terms re-slice or rebuild the shared
node. Cross-term reuse (CSE) of a genuinely shared node is therefore *not* operational
under external batching -- the very recompute the batched path is meant to avoid.

**(P2) The apparatus is an accretion of patches a per-slot mask subsumes.** The
scalar `scope_level` cannot distinguish "invariant outer to this loop" from
"loop-local at this loop" for a node carrying an External mode, which is exactly why
`consumed_upscope` had to be added (Tasks 1-3 of the phase-1 effort). A per-slot mask
encodes that distinction structurally (section 2), retiring the bit.

## 2. The per-slot lifetime mask

Define, for each **canonical** eval node (identity = `EvalExpr::hash_value()`, the
cache key), a mask over its canonical slots (`canon_indices()`), each slot classified
**SLICED** or **FULL**, by a **meet over all occurrences** of that node in the forest:

- **Occurrence-local classification.** In a given occurrence, a slot is *sliced* if
  its index is batched-External in the realized batch nest at/above that node (the
  accumulated `batched_here()` External stamps of ancestors + self). Contracted (aux)
  modes are summed away and never appear in `canon_indices()`, so only External (occ)
  loops classify slots. **Proto-aware:** a composite/PNO slot `a<i,j>` is sliced iff
  its proto pair `(i,j)` is batched -- a PNO index is domain-tied to its occ pair, so
  slicing the pair slices it.
- **Cross-occurrence meet.** A canonical slot is SLICED iff *every* occurrence slices
  it to a **proto-compatible** block; it is FULL if any occurrence leaves it unsliced,
  or occurrences slice it to incompatible blocks. "Proto-compatible" means the block
  is determined by the same (possibly proto-tied) batched index -- not raw canonical
  label equality. Using label equality here mis-fires on PNO slots (`a_3<i1,i2>` vs
  `a_4<i1,i2>` are the same sliced pair under different free-index labels); the
  physical rule is proto-aware block agreement.

**Derived quantity -- mask-scope.** The node's placement level is the deepest
(innermost) batch-loop level among its SLICED slots. All-full => mask-scope `-1`
(block-agnostic; the root/run cache). This subsumes the scalar `batch_scope_level`.

**Why the meet is the fix for P1.** The mask is one value per canonical node, so a
node has one placement regardless of how many terms use it. `s*C` has occurrences
that leave it fully unsliced => meet = all-full => one top-cache entry, reused across
all 52 terms. A mixed PNO node whose occ-pair slot is sliced in *every* occurrence
(to the same proto pair) stays mixed => per-pair scratch, correctly distinct per pair.
Proto-aware block compatibility is what separates these two -- resolving survey
caveat #1 (block-agnostic reuse vs per-occ-pair instantiation) and caveat #2
(PNO-slot label mismatch) together.

**Why the mask subsumes `consumed_upscope`.** A node whose *own result* carries a
sliced-External slot has that mode in its own SLICED set, so its mask-scope is *at*
that loop -- it is loop-local, placed at that level, never hoisted above it. The
scalar `scope_level` could not see this (hence the bit); the per-slot mask encodes it.

**Caveat -- the demoted-external carrier still needs a per-occurrence check.** The
subsumption above holds for the mask read *per occurrence*, but the canonical mask
`S(n)` is the cross-occurrence *meet* (Eq. below): a mode is demoted out of `S(n)`
when occurrences bind it to different (proto-incompatible) blocks, even though the
node carries it free in every occurrence. Such a *demoted-external carrier* has an
empty (or under-reporting) `S(n)` yet must NOT be hoisted to the run scope -- doing so
would materialize its full external extent (the cross-pair giant). The runtime
therefore keeps a per-occurrence cross-check: a member with empty combined residency
is hoisted only if it carries no External `batched_here()` stamp at all (the
`has_demoted_external` exclusion; see section 4). So the mask retires
`consumed_upscope` as an *emitted annotation*, but the loop-local/invariant distinction
it drew is recovered at runtime from `batched_here()` vs `S(n)`, not from `S(n)` alone.

## 3. Placement rule

A canonical node's home cache follows directly from its mask:

- **all-full (mask-scope `-1`)** -> the **root / run-scope cache**, block-agnostic,
  one entry shared across all terms.
- **mask-scope `L >= 0`** -> the **level-`L` scratch**, stored FULL over its full
  slots and its deeper (invariant) modes, SLICED over its slots up to and including
  level `L`.

Slice-on-use (phase 1) serves a node from its home scope to any deeper use, slicing
the modes the fetch crosses. The root builder's footprint gate (`max_footprint`) is
unchanged: an all-full node over the gate is still recomputed rather than held whole.

## 4. Mechanism: per-level placement (option A)

Retire `hoist_invariants`. Realize placement at each batch-loop level directly:

> At each realized loop level `L` (its block scratch created, before the deeper nest
> replays within a block), materialize -- exactly once per block of `L` -- every
> member-subtree node whose **mask-scope == `L`** on level `L`'s scratch, built via the
> sliced-leaf path (so it is sliced to this block and its own outer modes via
> slice-on-use, full over its full/deeper modes). Then descend. Deeper per-block
> replays find these via the scope chain (fall-through + slice-on-use) and do **not**
> rebuild them.

Nodes with mask-scope `< L` were materialized by an outer level; mask-scope `> L`
are left for their own (deeper) level; mask-scope `-1` are placed at the root by the
builder. Each level thus owns exactly its mask-scope-`==L` nodes -- no strictly-outer
walk-up, no fresh-cache pre-build triggered from an inner loop, no `consumed_upscope`
special case. A mask-scope-`==L` node is `L`-sliced and invariant across `L+1..inner`,
so building it once per `L`-block and reusing it across inner blocks is exactly
correct.

This is the difference from `hoist_invariants`: hoist looked *down* from a firing
inner loop for outer invariants and walked *up* to store them; per-level placement has
each level build *its own* nodes and store locally. The scatter and group/contracted
branches keep their per-block replay structure; only the `hoist_invariants(...)` call
inside each is replaced by the per-level materialize step.

## 5. What retires, what stays

**Retire:**

- `hoist_invariants` (the function, its strictly-outer collect predicate, and the
  walk-up target-location), in `eval.hpp`. The register-then-store primitive
  (`ensure_hoist_slot` + `store` on a scratch whose key set is otherwise seed-only) is
  **reused** by per-level placement -- local to each level, no walk-up.
- Veto **part (b)** only (`n->batch_scope_level() >= 0`) in `cache_manager.hpp`;
  replaced by mask-based run-scope membership (all-full => eligible).
- `consumed_upscope`: the `NodeBatchAnnotation` field, the `EvalExpr` accessor
  (`batch_consumed_upscope`), the cost-model emit, and the `!ext_loop_local`
  hoist-exclusion (all landed in Tasks 1-3 of the phase-1 effort).
- The scalar `batch_scope_level` emit and its runtime uses; superseded by the mask and
  its derived mask-scope.

**Keep:**

- Slice-on-use (phase 1): `access_at` `{ptr, hops}`, the Enter-stage `slice_to_use`,
  `sliced_leaf`, the `batch_context`.
- Veto **part (a)** (`sliced_batch_axis`, contracted mode free on result) -- a distinct
  mechanism (a contracted-accumulation free-large-index intermediate), not a placement
  scope.
- The `max_footprint` gate, the scope chain (`set_parent`/`parent`/`access_at`),
  `make_batched_scratch`, and the scatter / group-contracted per-block replay branches.

## 6. Payoff: CSE operational across terms

With one consistent mask per canonical node, a genuinely shared intermediate lands in
one home cache and is reused everywhere. `s*C` (52 occurrences) becomes one top-cache
entry instead of being re-sliced in 24 terms; shared mixed intermediates land once at
their occ-pair scratch and are reused across the terms that share that pair. This is
the "CSE if we batch over externals" the phase-1 discussion flagged as in-scope: the
per-canonical mask is the mechanism that makes cross-term reuse under external batching
correct, because placement no longer depends on which occurrence the runtime met first.

## 7. Scope boundary

This spec changes **emit** (the mask, replacing the scalar scope annotation and the
`consumed_upscope` bit) and **runtime** (mask-based membership + per-level placement,
replacing veto-b + hoist). It does **not** change the DP's cost/sizing (`szcell`,
cell-restricted sizes). "Cost-model lifetime-awareness" -- having the DP price
intermediates by the mask's residency -- is a separate later spec; the C60 gate shows
the DP already sizes the chosen schedule accurately post-slice-on-use, so a placement
change is not expected to need a simultaneous sizing change. If it does, the gate
(section 8) catches it as a modeled-vs-replay divergence.

## 8. Acceptance gate and risks

**Gate (measurement-driven, same instruments as phase 1):**

- C60 `tests/unit/test_eval_dryrun.cpp` `[.][dryrun-occ-veto]` (`SEQUANT_UT_DRYRUN_NTERMS=55`,
  `cmake-build-release`): replay peak **<= ~443 GB** (same-or-better; the mask-placement
  is meant to preserve hoisting, so `avoidable_time` should be **no worse, ideally
  better** as cross-term CSE removes re-slicing of `s*C` and shared mixed nodes).
- MPQC he10 CSV-CCSD with forced aux+occ batching: energy correct
  (Δ < 1e-7 vs `-0.33231474200227867`; occ-batching's accumulation-order Δ ~= 6.6e-8 is
  expected and unrelated).
- **OFF-path byte-identity:** on the order-blind path every node is all-full /
  mask-scope `-1`, the veto and per-level placement are inert, and the per-batch replay
  runs exactly as before. This must be byte-identical (a unit assertion, as in phase 1).

**Risks:**

1. **Eval-order change.** Per-level placement builds mask-scope-`==L` nodes at level
   `L` before descending; this reorders when some intermediates are formed relative to
   hoist. Numerics (accumulation order) and peak can shift. The C60 gate + energy check
   bound this; expect only accumulation-order-level deltas.
2. **Per-block scratch interaction.** A mask-scope-`==L` node is per-`L`-block; the
   materialize step must build it inside the per-block loop (once per block, on the
   reused scratch) and must not leak a wrong-block value into `L+1`. This is the
   subtle implementation point; the plan must pin the exact hook site in each replay
   branch and test wrong-block isolation.
3. **Cross-occurrence meet correctness.** The meet is a new forest-wide, per-hash
   analysis. Proto-aware block compatibility must be right or `s*C`-class nodes get
   mis-placed (over-sliced full-extent storage, or a shared node re-sliced). The
   survey's fully-characterized nodes (`1989507463377952644` all-full;
   `15545560759149115397` mixed at occ-pair scope) are the ground-truth fixtures.

## 9. Open questions / deferred

- **Cost-model lifetime-awareness** (DP prices intermediates by mask residency) --
  separate later spec, as above.
- **Reducing the 443 GB** (summand-12-class peak) -- needs a new batching lever
  (contracted-occ, or cross-term forest merge), tracked separately.
- **The mask as the single source for both placement and cost** -- once the cost-model
  spec lands, the per-slot meet should be computed once and consumed by both the DP and
  the runtime (the "shared analysis" from the phase-1 discussion). This spec computes
  it for placement only; the interface should be written so the DP can consume the same
  result without recomputation.
