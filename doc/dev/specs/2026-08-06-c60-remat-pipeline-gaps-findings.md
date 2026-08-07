# C60 batched-eval peak non-monotonicity: the remat pipeline is disconnected

_2026-08-06. Branch `evaleev/feature/multimode-batched-eval`._

## Summary

Pricing C60 PNO-CCSD (pVDZ-F12 doubles residual) through the dry-run cost
model shows the modelled peak is **non-monotonic in `peak_threshold`**: a
tighter budget can yield a *higher* realized peak.

| `peak_threshold` | forest peak | fits 1 TB |
|---|---|---|
| 100 GB | 333 GB | yes |
| 500 GB | **3152 GB** | **no** |
| 1000 GB | 813 GB | yes |

Adding batching should never raise the realized peak under full optimization,
so a heuristic is failing. Root-causing it (on the single peak-driving term,
#4) shows the fault is **not** in the cost model's sizing, the factorizer, the
proto/ToT index handling, or the batched evaluator's slicing. It is that the
**rematerialization pass that is supposed to enforce the peak budget is not
connected end-to-end.** Three independent gaps, all measured.

## Methodology note (important)

The per-node "realized size" numbers must come from the authoritative source —
`cost_profile()`'s DryRun replay trace, which sizes each op via
`EvalExpr::canon_indices()` (the real array-mode slots) and the actual runtime
`ExtentOverrides`. A hand-rolled "crude walk" over `node_free_indices()` (which
reads result modes off the *expression*, `as_tensor().const_braketaux_indices()`,
and hides a ToT composite's proto/outer occ slots) produces wrong per-node sizes
and led to two incorrect diagnoses before being discarded. Do not size IR nodes
off expression reads; use `canon_indices()`.

## The peak driver (term 4, thr=500, authoritative)

The peak is driven by a **g·C intermediate** `I(i_1,i_2,i_3,Κ_2;a_3<i_1,i_2>)`
built **once, full, 2509 GB** (`builds=1`) and held resident through the whole
225-iteration (15×15) `i_3,i_4` loop. At thr=1000 the analogous g·C is built
**per-Κ_2-block (60× at 41.8 GB)** and the peak is 757 GB.

Why full at thr=500: the g·C is **CSE'd** — its two occurrences (the `a_3`
copy binding `i_3` and the `a_4` copy binding `i_4`) share one canonical
identity (`hash=15545560759149115397`). `stamp_lifetime_masks` computes each
node's residency as the cross-occurrence *meet* of its batched modes; here that
is `{i_3} ∩ {i_4} = ∅`, so its `home_scope` (`sliced_modes()`) is empty and it
is homed at **run scope (full)**. This is *correct* CSE-safe behavior — one
materialization must serve both occurrences, and it cannot be sliced by `i_3`
for one and `i_4` for the other.

The DP prices the same factorization at **368 GB**: its `subtree_peak`
propagates the enclosing `{i_3,i_4}` sliced-set into `U` for the descendant g·C
and shrinks it, i.e. it prices per-occurrence slicing that CSE forbids. This
per-occurrence-sliced peak *is* realizable — by recomputing the `i_4` slice of
the second occurrence inside the `i_3` loop (trading time for peak). **That is
exactly rematerialization's job**, and it is where the pipeline breaks.

## The three gaps (measured)

| thr | no-remat realized | remat pass output | remat router passed → realized |
|---|---|---|---|
| 500 | 3079 GB | status=Feasible, modelled **191 GB**, router non-empty | **3079 GB (unchanged)** |
| 100 | 240 GB | status=Feasible, modelled **66 GB**, router empty | 240 GB |

**Gap 1 — remat is not invoked; the budget is decorative.**
`dryrun::cost_profile(...)`'s `router` parameter defaults to `nullptr`, and
`remat_to_router` is called only from `test_placement_remat.cpp`. No production
or measurement caller builds a router, so `peak_threshold` is enforced by
nothing. It should not be possible to default to "no remat" — the threshold is
a budget.

**Gap 2 — a populated router is not applied by the evaluator. [FIXED 2026-08-07]**
Building a remat router at the 500 GB budget yields `status=Feasible`, a
non-empty router, and a remat-modelled peak of **191 GB** (it demotes the
run-scope g·C). But passing that router to `cost_profile` initially left the
realized peak at **3079 GB, unchanged**. Root cause: the hoist consults the
CONTEXT-DEPENDENT `occurrence_key` via `route()`. At the OUTERMOST hoist
(`ectx={}`, before the i_3 loop opens) the query misses (the override is keyed
at `{i_3}`), so the g·C is built FULL at the chain root first; inside the i_3
loop `route()` hits (rl -1→0) but the value is already alive at the broader
scope, so `if (target->alive(d)) continue;` skips the per-block rebuild.

Fix: a context-invariant `PlacementRouter::moved(value_hash)` flag (populated by
`remat_to_router` alongside the per-occurrence overrides). In the hoist, when
`route()` misses at the current context but the value is `moved()`, DEFER it
(do not build it full here) so it is built per-block at its deeper home where
`route()` hits. Byte-identical for an empty router. Result: term-4 thr=500
realized peak drops **3079 → 569 GB** and **569 < 757** (thr=1000), restoring
monotonicity. (The residual 569 vs remat's modelled 191 is Gap 3.)

**Gap 3 — remat's peak model disagrees with the realized peak.**
At thr=100 remat reports the *seed* placement already Feasible at **66 GB** and
moves nothing (empty router) — but the realized seed peak is **240 GB**. So
`peak_profile_sweep` (remat's internal oracle) under-estimates what the DryRun
replay realizes, meaning remat's own feasibility test is against the wrong
number.

## Fix direction

Complete the pipeline so budget → remat → router → runtime → realized peak are
one consistent loop:

1. **Invoke remat by default** in `cost_profile` (and the real eval path) at
   `peak_threshold`, not via an opt-in router argument.
2. **Apply the router in `evaluate()`** — consume `route()`/`home_depth` to home
   a moved value at its overridden (deeper) scope, so the plan is realized.
   (Gap 2, DONE 2026-08-07: the hoist now defers `moved()` values on a
   context-miss; 3079 → 569 GB.)
3. **Reconcile the two peak oracles** — `peak_profile_sweep` (remat) and the
   `cost_profile` DryRun peak must be the same measure, so remat optimizes and
   checks against the peak that is actually realized. (Gap 3.)

## What it is NOT (ruled out by measurement)

- Not the cost model's sizing (`inner_aware_volume` sizes proto as outer modes
  correctly).
- Not a proto-vs-free slicing bug — `canon_indices()` exposes `i_3,i_4` and
  `index_position` finds them; the evaluator slices free occ fine (thr=1000
  builds the g·C per-block).
- Not a degenerate loop — `eff_count=1` is the residual's accumulate-once count;
  the loop iterates 225× (inner nodes `builds=225`).
- Not the DP being CSE-blind in a way that is unrealizable — the DP's
  per-occurrence-sliced plan is realizable via remat recompute; remat just never
  runs / applies.

## Reproduce

Hidden test `[.][c60-term-dump]` (uncommitted spike harness), env:
`SEQUANT_UT_REGIME=c60 SEQUANT_UT_DRYRUN_TERM=4 SEQUANT_UT_PEAK_THR_GB=<gb>`
plus `SEQUANT_UT_C60_CPTRACE=1` (authoritative per-op sizes),
`SEQUANT_UT_C60_FULLTREE=1` (schedule + residency + hashes),
`SEQUANT_UT_C60_REMAT=1` (build+pass a remat router).
