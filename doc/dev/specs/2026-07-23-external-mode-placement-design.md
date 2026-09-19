# External-mode placement in the order-aware batching cost model (stage 3)

**Status:** designed 2026-07-23 (brainstorming). Stage 3 of the order-aware
multilevel batching cost model (`2026-07-20-order-aware-multilevel-batching-design.md`,
whose contracted peak+flops halves are implemented and validated). Supersedes the
root-level external seeding of `2026-07-20-external-mode-batching-design.md` with
node-level placement decided INSIDE the cost model.

**Repos:** SeQuant (`core/optimize/cost_model.hpp`). Consumed by MPQC CCk; any
behavior change needs an MPQC repin. Gated behind `order_aware_recompute` +
`batch_spectator_indices` (default off => byte-identical).

---

## 1. Problem

The C60 PNO-CCSD residual's memory blocker is a single ~2231 GB operand,
`g(μ̃,μ̃,Κ_2) * C(i_2,i_1,μ̃; a_1<i_1,i_2>)`, full in the free PAO index `μ̃` and the
occ pair `(i_1,i_2)`. These are **external** modes -- open on the term result,
contracted at no node. The contracted cost model (resident-scan peak + ordered-key
flops) fixes recompute (a TIME problem) but cannot move this peak (a MEMORY problem):
no contracted reordering touches an externally-indexed operand. Bounding C60 memory
REQUIRES batching external modes on the peak-setting operand.

The existing external machinery (`seeded_forest_peak`, `reconstruct_batched_modes`)
cannot do this. It is:

- **root-level only** -- it seeds an external mode into the batch context at the
  ROOT (the loop wraps the whole term), never at a subtree;
- **work-neutral-or-decline** -- it vetoes an external mode unless EVERY node carries
  it (no enclosed non-carrier), the opposite of how a contracted mode is treated (a
  contracted mode with a non-carrier just pays: hoist or recompute);
- **outside the DP** -- placement is decided in the emit walk AFTER the DP has priced
  the schedule, so the factorization is chosen blind to external batching.

On C60 it emitted 24 `External` stamps across 55 terms and left the 2302 GB peak
untouched (D5/A0 forensics).

**Array-level framing (removes a false blocker).** Factorization operates on ARRAYS,
not tensors: no bra/ket/aux, no protoindices, no symmetry -- only outer modes and
inner (ToT) modes. An outer mode from a protoindex slot is indistinguishable from one
from a plain slot; both are plain sliceable outer modes. So the occ pair `(i_1,i_2)`
appearing as protoindices of the PNO composite `a_1<i_1,i_2>` is NOT a special case --
they are ordinary sliceable external modes.

---

## 2. The unified model: external is a contracted mode placed at a chosen node

**`B` at a node is every mode whose batch loop encloses that node -- contracted OR
external, no distinction.** Once a mode is in `B`, the identical machinery prices it:
`escaped_outer` charges recompute for any enclosed node that does not carry the mode,
and the resident-scan prices hoisting that node out instead. `i` hoisting `I2` (below)
is the same event as a contracted mode hoisting a non-carrier.

The only asymmetry is **identification and placement**, not treatment:

- a **contracted** mode is summed at exactly one node; its loop wraps that node's
  subtree (placement fixed by the tree);
- an **external** mode is summed at NO node; the DP must CHOOSE which node its loop
  is placed at (an external mode is a contracted mode whose "contraction subtree" is
  a placement decision).

**External loops are LOCAL, not root-reaching.** A node can be evaluated one slice at
a time; placing an external loop around a subtree does NOT propagate to the top of the
tree. So a per-node external placement's peak benefit is LOCAL to that node's subtree.
This is what makes "place it where the peak is" clean -- the win lands where you
triggered it.

### Worked example (the R{i,m} note that grounds this)

`R{i,m} = (A{i,j} B{j,k})(C{k,l} D{l,m})`, tree `I1{i,k}=Σ_j AB`, `I2{k,m}=Σ_l CD`,
`R=Σ_k I1 I2`. Batching contracted `k` (outer) and external `i`:

```
for K:                        # k (contracted) OUTERMOST (carried everywhere)
  for kk,l,m:  I2[K] += C[K] D          # I2 hoisted ABOVE the i-loop, inside k
  for I:                      # i (external) INNER to k, placed LOW
    ... I1[I,K] ...  R[I] += I1 I2      # I2 (no i) reused per i-block
```

`i` is external (on `R`) yet its loop sits **inner to `k`**, placed low, with the
non-carrier `I2` hoisted to the `k`-level. `[k outer, i inner]` vs `[i outer, k inner]`
is a real cost difference the DP must price -- exactly like it now prices contracted
order. And "external batch loops do not just go to the top" means: don't naively
enclose everything; hoist the non-carriers -- the same hoist-vs-recompute trade.

---

## 3. Two-phase design (external decided AFTER contracted)

Deciding external placement jointly with the contracted DP is circular (a node's peak
depends on whether you inject external there). **Stage the decision:**

### Phase 1 -- contracted ordered DP (already implemented)

Run the contracted order-aware DP unchanged. It fixes the contraction tree, the
contracted batching, and therefore every node's PEAK. `build_cells` enumerates
**contracted** batchable modes only (external modes are excluded from the phase-1
enumeration -- today they sit in `ctx.batchable_modes` but are never injected, so
their cells are dead; phase 1 should not enumerate them).

### Phase 2 -- external placement on a fixed-contracted schedule

Because external batching is applied to a schedule whose CONTRACTED batching is
already fixed, that schedule's node peaks are fixed, so external placement is **unique
and deterministic** (no fixpoint). Given such a schedule, the pass is:

- **Trigger (per-node `peak_threshold`).** At each node whose modeled peak exceeds
  `peak_threshold`, consider batching. (This is the NEW per-node use of the budget;
  see Section 4.)
- **Try.** For each external mode the node carries, produce an injected variant:
  append that mode to the node's cell as INNER (via the same `descend` contracted
  modes use -- so its nest position, inner to every contracted loop enclosing the
  node, is determined by the injection node, not searched). Re-price the subtree's
  peak with the external mode in `B` -- `escaped_outer` charges the enclosed
  non-carriers, the resident-scan prices their hoist.
- **Keep.** `pareto_insert` the variant; if it does not cut the node's peak it is
  dominated and dropped.
- **Iterate for cascade.** Slicing `μ̃` at the peak node may drop it below a sibling
  that is now the peak and wants `i` sliced. Re-scan for over-threshold nodes after
  each placement; repeat until nothing exceeds budget or no external mode helps. A
  loop, not a combinatorial search.

**Where the pass runs -- two variants (start cheap):**

- **CHEAP (recommended first): process the ROOT frontier before `select_root`.** After
  phase 1 produces the root frontier (candidate whole-term schedules), apply the pass
  to each frontier point (walk its tree, place externals), producing externally-
  augmented `(peak, flops)` points; then `select_root` picks among them. Simple, cheap
  (few schedules). Blind spot: a schedule dominated on contracted `(peak, flops)` and
  pruned before phase 2 -- even one that would have been best WITH external -- is
  already gone.
- **COMPLETE: apply the pass to every candidate inside `relax` before `pareto_insert`.**
  External then influences frontier membership at every level, removing the "tree
  chosen blind to external" residual (Section 8). Costs external pricing at every DP
  candidate; more complex.

The choice IS the "does tree-blindness matter" question: cheap accepts it, complete
removes it. **Start cheap and measure** -- for C60 the peak-setting operand (`g·C`)
sits in essentially every reasonable factorization, so the schedule carrying it is on
the root frontier and the cheap pass reaches it. Only if the C60 peak stays stuck
because the winning schedule was pruned pre-external do we pay for the complete variant.

**Faithful peak credit is an annotation CONTRACT, not a DP guarantee.** A placed
external mode's credited peak drop is realized only if the runtime evaluates that
node's subtree one external-slice at a time (streamed assembly), never held whole. The
DP has no way to ENFORCE this -- and it never could for contracted batching either.
The DP emits the placement annotation; a well-behaved runtime is REQUIRED to respect
it, exactly the annotate-then-obey contract the resident-scan and contracted batching
already rely on.

---

## 4. Model parameters (precise definitions)

- **`peak_threshold`** -- the memory budget, in bytes. TWO well-defined uses:
  - **per-schedule** (phase 1, unchanged): `select_root` selects, among a term's
    candidate schedules, the min-flops one whose ROOT peak fits the ceiling (perf-first;
    fall back to global min-flops if none fit). A term that fits unbatched keeps the
    unbatched schedule -- a finite threshold gates selection, it does not force batching.
  - **per-node** (phase 2, new): a node whose modeled peak exceeds it triggers external
    injection at that node.
- **`order_aware_recompute`** -- master gate for the order-aware model (resident-scan
  peak + ordered-key flops + this external placement). Default off => byte-identical.
- **`batch_spectator_indices`** -- gate for external batching specifically (external
  modes are identified by `is_external_mode`: open on the root, contracted nowhere).

---

## 5. Relationship to existing code

- **Replaces `seeded_forest_peak` + its emit-walk trial loop** (root-level,
  work-neutral-or-decline, outside the DP) with the phase-2 node-level placement
  (priced inside the cost model, allowing enclosed non-carriers via hoist-vs-recompute).
- **`build_cells`** enumerates contracted modes only (phase 1); external modes are
  appended on demand in phase 2 via `descend`.
- **`escaped_outer` / resident-scan / `descend`** are reused verbatim -- an external
  mode in `B` is priced identically to a contracted one.
- **`BatchModeType::External`** emission stays; the phase-2 pass stamps it for the
  modes it places.

---

## 6. Constraints

- **Byte-identical when off.** `order_aware_recompute` and/or `batch_spectator_indices`
  off => phase 2 does not run => byte-identical to today.
- **Annotation contract (not a DP guarantee).** The DP emits placement annotations; a
  well-behaved runtime is REQUIRED to respect them (evaluate a placed node's subtree
  one external-slice at a time, streamed assembly). The DP does not and cannot enforce
  this -- the same annotate-then-obey contract as contracted batching and the
  resident-scan.
- **Cross-repo.** Validate on the C60 dry-run (DP-side modeled peak, the trustworthy
  metric -- NOT the replay `avoidable_time`) before any MPQC repin / Owl run.

---

## 7. Success criterion

On the C60 giant (fixture summand 38), with `order_aware_recompute` +
`batch_spectator_indices` on, the phase-2 pass PLACES external modes (`μ̃`, the occ
pair) on the peak-setting operand and its modeled peak DROPS (by roughly
block/extent per placed mode) -- the DP now REACHES the operand the root-level seeding
missed. Stretch: the giant's modeled peak moves under (or materially toward) the
budget. Measured DP-side via `cost_profile`, not the replay.

---

## 8. Known limitations (carry forward)

- **Tree chosen blind to external** -- the CHEAP phase-2 variant (Section 3) processes
  only schedules that survived contracted-only frontier pruning, so a tree that is
  min-flops-under-contracted-peak but not the best admitter of external placement wins.
  The COMPLETE variant removes this. Deferred BY MEASUREMENT: start cheap; upgrade only
  if the C60 peak stays stuck because the winning schedule was pruned pre-external.
- **Greedy per-node placement** with cascade iteration (Section 3) -- not a joint
  optimum over interacting externals; the frontier + iteration cover the common case.
- **Depth cap** (from the contracted work) also bounds how many modes co-nest at a
  node including placed externals.
