# CSE-aware rematerialization: splitting cells along relabeled modes

_2026-08-07. Branch context: `evaleev/feature/multimode-batched-eval` (SeQuant), paired MPQC batched-CCk work._

## Summary

The batched-eval peak model is non-monotonic in `peak_threshold`: a tighter
budget can raise the realized peak. Root-causing it end-to-end reveals a
**latent correctness hole**, not merely a mispricing: the rematerialization
placement router can direct the two occurrences of a **CSE-folded, divergently
relabeled** intermediate to share **one sliced materialization**, feeding one
occurrence the wrong slice. The runtime's tree-recursive hoist path
(`place_at_this_level`) has no guard against this, even though a *different*
runtime path (`make_batched_scratch`) already contains exactly the criterion
that would prevent it.

The fix is to make one criterion govern both paths: a value may be materialized
once and shared across occurrences only when those occurrences are
**slicing-signature consistent** at the shared scope; otherwise it must be
**split** — materialized per occurrence (recompute). remat must obey the same
criterion when it emits router overlays and when it prices peak.

This document records the model (canonical vs physical indices under CSE), the
DAG-scope framing that makes placement well-defined across trees, the verified
location of the hole, the unified fix, and an explicit scope boundary
(intra-tree; a scope-oriented executor is future work).

## Glossary

Reuses the batched-DAG vocabulary of
`2026-08-02-home-scope-placement-design.md` (§3-4); terms introduced by *this*
document are marked _(new)_.

- **Value** — what SeQuant hashes (`TreeNode` / `EvalExpr::hash_value()`): a
  canonical (sub)expression. Perfect CSE is "one value = one node," so the
  reused evaluation is a DAG of values.
- **Occurrence / instance** — a single *use-site* of a value: a position in a
  specific term's tree. One value has many occurrences (what CSE folds together).
- **Cell** — a *chosen physical copy* of a value serving a **subset of its
  occurrences**: one build, one live range, one cache entry — a partition-cell of
  the value's occurrences. `cell = (value, home-scope)`: the two identify a cell
  on their own, because split cells of one value land at *distinct* home-scopes
  (the DAG-scope of each occurrence subset), so no separate split-index is needed.
- **Split** — un-folding a value into more than one cell (a *partial un-CSE*),
  each serving a disjoint occurrence subset, when one shared copy cannot serve
  all. The central mechanism of this document.
- **Home scope** — which level of the batch-loop nest a cell is built at and
  lives at. A cell **carries** (slices) the home loops it depends on and is
  **invariant to** (rebuilt/replicated within) the rest — the latter is its
  recompute cost.
- **Canonical indices** (`EvalExpr::canon_indices()`) — a node's *physical* index
  labels in canonical *order* (labels preserved). Distinct from the structural
  hash, which is label-agnostic.
- **Occurrence key** — a **DAG-global** identity of an occurrence's *batched-slot
  structure* (which slots batched, of which space), independent of tree-local
  labels; what router overlays are keyed by. Detailed in Background.
- **DAG-scope** _(new)_ — a schedule-global batch-loop **nest position**: an
  **ordered sequence of spaces** (a nest prefix), e.g. `[occ, occ, aux]`,
  independent of any value or tree. A level is named by its **position** in the
  sequence, which disambiguates a space repeated *along the path* (`occ→occ`
  nesting) — so no separate instance index is needed **given our construction has
  no same-space siblings** (see "DAG-scopes form a tree" for that assumption and
  when it would break). The label-free replacement for a tree-local home.
- **Shared-label vs relabeled mode** _(new)_ — a mode whose physical label is
  identical across *all* a value's occurrences (sliceable *while folded*) vs one
  whose label differs across occurrences (sliceable *only by splitting*).
- **Slicing signature** — per occurrence, the position of the batch mode in
  `canon_indices()`, or its absence; the runtime shares a materialization across
  occurrences only when their signatures are consistent (`make_batched_scratch`).
- **GCP (greatest common prefix)** _(new here as a placement rule)_ — the deepest
  DAG-scope shared by all a value's consumers; a shared cell's home is bounded by
  it, and a too-shallow GCP forces either a full (top) home or a split.

## Background: two identity notions

An eval-IR node has two distinct identity/label notions, and conflating them is
the source of the bug:

- **`EvalExpr::hash_value()`** — the **structural DAG-node hash**: it identifies
  a node by tensor-network topology up to α-equivalence of index labels.
  `g·C(i_1,i_2,i_3,Κ; a_3<i_1,i_2>)` and `g·C(i_1,i_2,i_4,Κ; a_4<i_1,i_2>)` —
  same connectivity, different labels — **collide**. This is what CSE folds by.
- **the occurrence key** (`TensorNetwork::SlotCanonicalizationMetadata`, keyed
  in `PlacementRouter`) — a **DAG-global** identity of an occurrence's
  *batched-slot structure*: which slots are batched, of which space, up to the
  value's symmetries. It is what the router keys overlays by, so it **must be
  independent of tree-local index labels** — two occurrences of one value in
  different trees, or two legs in one tree bound to different externals, must map
  to the **same** key. It therefore colors a batched slot by its **space color
  combined with a single batched salt** — the salt is **the same for every
  space** (it only marks a slot as batched); the per-space *color* is what keeps
  batched-occ distinct from batched-aux. So all batched-occ slots share one
  color (distinct from non-batched-occ, renamable *within* the class), all
  batched-aux slots share another, and a non-batched slot is colored by space
  alone. The *specific* index bound to a batched slot (`i_3` vs `i_4`) is **not**
  encoded here — it is resolved later, in home-scope computation.

`EvalExpr::canon_indices()` is a third thing: the node's **physical** index
labels in canonical *order* (labels preserved, not relabeled). Proof, from a
dry-run residency dump:

```
hash=15545…  canon=[i_2 i_1 i_3 Κ_2 a_3<i_1,i_2>]
hash=15545…  canon=[i_2 i_1 i_4 Κ_2 a_4<i_1,i_2>]
```

Same `hash_value()`, different `canon_indices()`. So CSE folds by structural
hash while each occurrence retains distinct physical labels.

## The model: how CSE and slicing interact

CSE folds the two occurrences into **one value**, materialized **once**. This is
sound *because a full materialization's free indices are just names*:
`g·C[·,·,x,·,·]` is one array over the whole third-occ dimension `x`; the two
occurrences are two reads that *name* `x` differently (`i_3` vs `i_4`). Renaming
a free output index is a **no-op on data**, so one full array serves both reads.

The no-op-relabel property holds **only for a full materialization**. Slicing
the third-occ (home over that slot) binds the slice to a **physical loop
variable** — the `i_3`-block. The occurrence that relabels the slot to `i_4`
needs the `i_4`-block; relabeling the `i_3`-block's *data* to "`i_4`" yields
wrong numbers. Therefore:

> **A CSE-folded value cannot be sliced along a mode its occurrences relabel
> differently. Slicing such a mode requires *un-folding* — a separate
> materialization per occurrence (recompute).**

`stamp_lifetime_masks` already encodes this: the residency meet for such a value
is `∅` precisely because the relabeled mode does not survive the
cross-occurrence intersection. The meet is the correctness guard. The bug is
that **`shrink_candidates` bypasses it** — it offers the relabeled mode from
`(enclosing_modes ∩ carried) − home_modes`, read off the *first* occurrence's
physical labels (`placement_remat.hpp`).

### Shared-label vs relabeled modes

- **Shared-label modes** — identical `canon_indices` across *all* occurrences
  (`i_1`, `i_2`, `Κ`): sliceable *while folded* (survive the meet).
- **Relabeled modes** — differing `canon_indices` across occurrences (third-occ
  `i_3`/`i_4`, PNO `a_3`/`a_4`): sliceable *only by un-folding*.

## Scope must be a DAG-level coordinate, not a tree-local label

A home written as `{i_3}` is meaningful only inside the one tree whose dummy is
`i_3`. CSE merges trees into a DAG; a shared value's other uses live in trees
with their own labels. So a label-based home is not a well-formed property of a
DAG value. Scope must be **schedule-global**:

- A **DAG-scope** is an ordered sequence of **spaces** — e.g. `[occ, occ, aux]`
  — describing a batch-loop *nest position*, independent of any value. Sequence
  order carries the nesting (interleaving like `[occ, aux, occ]` is
  representable); a level is named by its **position**, which disambiguates
  repeated spaces without a separate instance index.
- **Per-occurrence translation** comes from the occurrence key: it encodes that
  occurrence's tree-scope (physical batch context, tree-local labels); matching
  each tree loop-index to the space at its **nest position** yields the map
  *(tree-index → DAG-scope position)*.
- A value's **home is a *prefix* of the nest**, and its own carried indices tell
  you, per level, whether it is *sliced by* that loop or merely *homed within
  it* (invariant → replicated → the recompute cost). Example: leg A homes at
  `[occ]` (position 0) and carries it (sliced on `i_3`); leg B homes at
  `[occ, occ]`, carries the loop at **position 1** (sliced on `i_4`) but **not**
  position 0 — homed within the `i_3` loop it doesn't carry, hence rebuilt per
  `i_3`-block.

Value-relative coordinates (enumerating the *value's own* slots) are rejected:
they cannot name a level a value is homed *within but invariant to* (B's `i_3`),
which is exactly the level that carries the cost.

### DAG-scopes form a tree

The set of all DAG-scopes is a **tree** — the loop-nest tree — rooted at the
empty (index-free) scope, with each edge opening one more loop. A cell's home is
a **node** in this tree; a value's consumers are cells at **descendant** scopes
(they slice the value on use). This is the structure the placement rule below,
the peak/co-residency model, and (eventually) execution all operate on: not the
contraction tree, but the nest tree the batching induces.

> **Coordinate caveat: the space-sequence names a node only when no level has
> same-space siblings.** A node of this tree is a *path* from the root, so the
> space-sequence (the edge labels along that path) is a unique coordinate **only
> if no node has two children of the same space**. Position disambiguates a space
> repeated *along one path* (`occ→occ` nesting — `i_1` at position 0, `i_2` at
> position 1) but **not** same-space *siblings*: if `i_1` opened two sibling occ
> loops `i_2` and `i_3`, both children key to `{occ, occ}` and the prefix names two
> distinct scopes. **Our construction cannot produce that**, and that is why the
> space-sequence suffices (and why no instance index is needed): the DAG is a
> CSE-fusion of independent trees, and the router resolves each `(space, position)`
> scope to an occurrence's *own* physical index per use, so tree-A's `{i_1, i_2}`
> and tree-B's `{i_1, i_3}` deliberately **collapse to one** `{occ, occ}` scope
> (that unification is the DAG-globality). We only ever get same-space *nesting*
> (the g·C `i_3`-outer/`i_4`-inner chain), never un-collapsed same-space siblings.
> A **future** executor that kept same-space-same-position loops distinct (loop
> fission, or a nest that genuinely branches) would violate this and require a
> path/instance disambiguator on the coordinate — the very "instance index" this
> design otherwise argues away. Note it before relying on space-sequence
> uniqueness in any new nest-forming pass.

### Share vs split = greatest-common-prefix + fusibility

A shared node's home is the **greatest common prefix (GCP)** of all its
consumers' DAG-scopes. Incompatible nests (`[occ, occ]` vs `[aux, occ]`) share
only `[]` → the node is homed full at top and each consumer slices on use —
correct, if unshareable. When the GCP is too shallow (forces full, too big),
**partition consumers into clusters that share a deeper DAG-scope** and home one
cell per cluster at *that cluster's* GCP, paying recompute per extra
materialization. "Slice at top / can't slice" is the degenerate one-cluster-at-
`[]` outcome. This consumer-clustering *is* the split, and it reduces to the
signature criterion below.

## The verified correctness hole

The runtime has **two** paths that materialize a value for reuse, and only one
is guarded.

**Guarded (safe): `make_batched_scratch`** (`eval.hpp` ~1443-1460). It registers
a subnode for sharing only if its **slicing signature** — *the position of the
batch mode in that occurrence's `canon_indices()`, or its absence* — is
consistent across all occurrences. Inconsistently-sliced subnodes are
"evaluated per occurrence, unshared." For the g·C over `i_3`,`i_4`: leg A carries
`i_3` (pos 2), absent `i_4`; leg B carries `i_4` (pos 2), absent `i_3` →
inconsistent for both modes → **not shared → per-occurrence → correct**. This is
why base runs and the reference energy are correct.

**Unguarded (the hole): `place_at_this_level`** (`eval.hpp` ~1801-1895), the
tree-recursive hoist path the router drives. It:
1. **dedups occurrences into one target** by canonical identity —
   `std::none_of(targets, eq(*p,n))` (1801-1802);
2. **registers one slot** — `target->ensure_hoist_slot(d)` (1883), keyed by hash,
   *not* through the signature check;
3. **builds once and shares** — `if (target->alive(d)) continue;` (1884),
   `evaluate_impl(d, sliced_leaf)` builds one materialization sliced to `d`'s
   home; every occurrence fetches it.

And the router-directed fetch `access_at_hops` (`cache_manager.hpp` 557-565) is a
`cache_map_.find(key)` by canonical identity at a fixed scope — **no signature
check**.

Consequence:
- **No router:** the divergent value's home is the meet `∅` → root → **full →
  shared safely**.
- **Router homes it at a sub-scope** (`{i_3}`): one **sliced** materialization is
  shared by both divergent legs → the `i_4` leg reads the `i_3`-slice →
  **corruption.**

Router-specific, matching observation exactly: the w8 full residual is correct
because remat did not move the g·C there; the isolated c60 term-4 (where remat
moves the g·C to `{i_3}`) is where the hole fires. `make_batched_scratch`'s guard
does not cover this path.

## Design

One criterion governs both paths and remat.

**0. Make the occurrence key DAG-global.** Change `occurrence_key`
(`occurrence_key.hpp`) to color batched slots by **space color + a batched
salt** instead of the current **space + label** (73-80). The salt is a **single
value shared across all spaces**; combined with the per-space color it makes
batched-occ distinct from non-batched-occ (and from batched-aux, via the space
color), while all batched-occ slots share one color and are renamable within the
class.
Then a value's occurrences map to **one** occurrence key regardless of which
external each is bound to (`i_3`/`i_4` in one tree, `i_7` in another), so the
router carries **one DAG-global overlay** per occurrence pattern, and the
specific binding is resolved per use (below). Without this, overlays are
keyed by tree-local labels and cannot be shared across occurrences/trees — the
whole DAG-scope approach reduces to "get lucky when labels line up."

**1. Unify the hoist path with the signature criterion.** In
`place_at_this_level`, an occurrence is only deduped into an existing target and
allowed to share its materialization when it is **slicing-signature consistent**
with that target at the home scope. When it is inconsistent (a relabeled mode
sliced at the home), register a **per-occurrence** hoist slot and build
**per-occurrence** (the split), keyed so the two materializations do not collide.
The signature test is the *same* one `make_batched_scratch` uses (position of the
home's batch mode in `canon_indices`, or absence); factor it into a shared helper
so both paths call one implementation.

**2. remat emits one DAG-global overlay; the home resolves per use.**
`remat_to_router` emits a **single** overlay per moved value, keyed by its
DAG-global occurrence key (§0), whose `HomeTarget` carries a **DAG-scope**
(ordered sequence of spaces), not a tree-local label set. At each use site,
`home_depth` resolves that generic scope through the occurrence's **physical
binding** of the batched slot — A's batched-third-occ is `i_3` → depth(`i_3`);
B's is `i_4` → depth(`i_4`) — so the two occurrences land at **different scopes**
(distinct `cache_map_`s → no collision). The split is thus realized by *one*
overlay resolving to different physical depths per use; no per-occurrence overlay
is needed. This requires `home_depth` (and the router read) to take the
occurrence key so it can translate a DAG-scope position → physical loop, replacing today's
raw `Index`-identity match of a canonical residency against the physical `ctx`
(`placement_router.hpp`).

**3. remat prices the split — as two real cells with a *replication factor*, not
a 2× fudge.** `shrink_candidates`/`apply_shrink` become CSE-aware: a candidate
mode is classified shared-label (slice in place, the existing behavior) vs
relabeled (offered only as a *split*). Applying a split **un-folds the one
`ValueCell` into two cells of the same hash**, each serving one occurrence subset,
each with its own home-scope, carried set, and liveness interval (subset-local
`first_use`/`last_use`); each split cell is non-divergent. This requires retaining
the per-occurrence records (`carried`/`ectx`/`point`/`consumer_point`/`home`,
computed per `NodeRec` in `compute_dag_boulevard`, discarded at grouping today) so
the split can partition occurrences by physical binding.

The cost of a split cell has **two** components, and the current 2× in
`cell_footprint` captures neither correctly — it is a peak-only fudge that
silently drops the dominant term:

- **Peak (memory).** The two split cells' footprints are co-resident only where
  their liveness intervals overlap; the interval sweep prices this structurally
  (double-counted on overlap, not on disjoint occurrences). Deleting the flat 2×
  and letting the sweep decide is correct *for peak*.
- **Recompute (the dominant, mis-priced term).** A split cell homed at a nest
  prefix is **rebuilt once for every iteration of each enclosing loop it is homed
  *within* but does not *carry*** (invariant → replicated → recompute, per "The
  model"). So its recompute cost is
  `build_cost × ∏ (block count of each enclosing-but-not-carried loop level)`,
  **not** a constant 2×. Canonical case (`i_3` outer, `i_4` inner; V bound to
  `i_3` in leg A, `i_4` in leg B): cell A homes at `i_3`, carries it → **1× full
  V**; cell B homes at `i_4`, carries `i_4` but is **invariant to `i_3`** → it
  materializes a full V **every `i_3`-block** → **N₃ × full V**. The split costs
  `(1 + N₃) × full V`, which equals the 2× fudge only at N₃ = 1 and diverges
  above it.

**Does this recompute term change the schedule? Not by itself — it is a
prediction.** remat operates over a **fixed DAG** and only chooses placement; no
pricing can alter what is computed or the numerics. The split itself is forced by
*correctness* (a divergent value homed at a sub-scope must un-fold, §1), not by its
price, so a mis-priced replication factor never yields an incorrect schedule. And
remat's placement is gated by **peak** (`cell_footprint`, the hard budget
constraint), which is a *memory* quantity — the replication factor is *flops* and
does not enter peak. `rematerialize_to_budget`'s selection objective today is
**pure `DeltaPeak`** (its doc comment: "all shrinks are zero-recompute this phase,
so the metric is pure `DeltaPeak`"), so the recompute term is **report-only**: it
sharpens the dry-run's cost *forecast* (does C60 fit, at what recompute?) without
moving a single home. The one thing it *breaks* is that invariant — a split is no
longer zero-recompute. **Whether to then make recompute a *secondary* objective**
(prefer, among peak-feasible placements, the one that recomputes least — so remat
never picks an expensive split when a free shared-label shrink would fit) **is a
separate, deliberate decision, not implied by pricing the term.** Until it is made,
price the replication factor for the forecast and keep the objective peak-driven.

**Hard dependency on §2.** This pricing is **impossible on a value-relative home
coordinate** (`slot_positions`) and *requires* the DAG-scope home of §2. The
replication factor is a product over the levels a cell is *homed-within but does
not carry* — exactly the levels (leg B's `i_3`) that a value's own slot positions
**cannot name** (`i_3 ∉ B's canon_indices`; see "Scope must be a DAG-level
coordinate"). A DAG-scope home (an ordered nest prefix) does name them, so §3 is a
downstream consumer of §2: **§2 must land before §3 can be priced correctly.**

Consumer-clustering (GCP per cluster) is the general form; the minimal realization
is per-occurrence split along relabeled modes. `remat_to_router` still emits **one
DAG-global overlay per value hash** (§2): both split cells share the hash and
resolve to the *same* DAG-scope, so moved-detection must become hash-aware (split
cells carry fresh `value_id`s absent from the seed map) while still emitting one
overlay.

**4. Defense-in-depth on the router read.** Have the router-directed read verify
the fetched entry's slicing signature matches the occurrence's before using it,
refusing a mismatched share (recompute instead). This is the runtime analogue of
`SEQUANT_ROUTER_SHADOW`, but a live correctness check rather than a no-op
assertion, so no future overlay can silently reintroduce the hole.

## Scope boundary

**In scope (correctness, this design):** the split is realized on the **current
tree-recursive + scope-chain + router** runtime. The divergence this fixes is
**intra-tree** — the g·C's two legs are operands of one contraction in one term —
and the split cells land at **distinct scope-chain depths** within that tree, so
the existing scope chain holds both. Cross-tree divergence is already correct:
scratch caches are per-pass (non-persistent), so the only entries surviving
across trees are persistent-at-root (full, safe); recursion recomputes per tree
rather than sharing a sliced value.

**What recursive one-tree-at-a-time descent sacrifices.** The current executor
descends each term's tree in turn. A sub-top scratch scope opened while
descending one tree is closed before the next tree runs, so a *sliced*
intermediate materialized there cannot be seen by any other tree. Genuine CSE at
a sub-top DAG-scope means "materialize once, sliced, and share across *all*
trees contributing to that scope" — which requires those trees to be
**co-evaluated inside that scope's loop**. Recursion cannot do that. So it keeps
CSE only at **two** places: the **top** (persistent/root) scope, where the value
is full, and **within a single tree**. **All cross-tree CSE below the top level
is sacrificed** — such values are either homed full at top or **recomputed per
tree**. This is not a correctness loss (recompute is exact), but it is real:
most of the potential sub-top CSE across the forest is given up.

**Out of scope (future work):** a **scope-oriented executor** that walks the
DAG-scope (loop-nest) tree and, at each scope, opens the loop **once** and
computes all cells homed there **across every contributing tree** — recovering
exactly the cross-tree sub-scope CSE recursion sacrifices. Since schedules vary
greatly with problem size, we cannot assume this is never needed: at scales
where fitting the budget *requires* sharing a sliced intermediate across trees,
recursion's recompute is correct but simply **cannot fit**, and no split rescues
it — that is the case that makes the scope-oriented executor mandatory. It is a
separate, larger project and is **not** required to close the correctness hole
this document fixes (which is intra-tree).

## Testing

1. **Occurrence-key DAG-globality**: two occurrences of one value bound to
   different externals (both intra-tree, e.g. `i_3` vs `i_4`, and across trees,
   e.g. `i_7`) produce the **same** occurrence key under space+salt coloring; a
   batched-occ slot and a non-batched-occ slot produce **different** colors.
2. **Signature-criterion unit test**: a minimal forest with a CSE-folded value
   whose two occurrences slice a relabeled mode; assert the hoist path builds it
   **twice** (per-occurrence) and never shares a divergent slice — mirroring the
   `make_batched_scratch` behavior.
3. **Router no-op invariance**: empty/no-op router remains byte-identical
   (existing guarantee); the new per-occurrence path is inert when nothing is
   split.
4. **Wet correctness under a forced split**: drive a small system (w8 PNO-CCSD)
   through a schedule that *does* move the divergent value to a sub-scope (via
   the schedule-size override that makes remat fire), and assert the CCk energy
   matches the reference to lossless tolerance (~1e-8). Without the fix this run
   must diverge (the corruption); with it, it must track.
5. **Peak monotonicity**: the modeled and realized peak are non-increasing as
   `peak_threshold` tightens, on the C60 term that exposed the non-monotonicity.
6. **Defense-in-depth check**: with the router read's signature verification
   enabled, a deliberately mis-emitted shared overlay is caught (recompute), not
   silently corrupting.

## Non-goals (YAGNI)

- No `(hash, home-scope)` cache re-keying: divergent occurrences land at distinct
  scope-chain depths, so the hash-keyed cache already holds both.
- No scope-oriented executor (see boundary).
- No change to the structural hash or the occurrence-key coloring; the fix is to
  make placement/sharing *obey* the criterion those already express.

## Reproduce / evidence

- Divergent g·C structure and per-occurrence `canon_indices`: SeQuant dry-run
  residency dump (`[c60-term-dump]`), term 4.
- Hole location: `eval.hpp` `place_at_this_level` (dedup 1801-1802, register
  1883, build/share 1884-1895); `cache_manager.hpp` `access_at_hops` (557-565).
- Guard that the hoist path lacks: `eval.hpp` `make_batched_scratch` (~1443-1460).
- Router miss/defer seam already present: `eval.hpp` Enter-stage router consult
  (718-776).
- Wet confirmation that base runs are correct (guarded path) and the hole is
  router-specific: w8 PNO-CCSD reference `CSV-CCk Energy = −2.0005382446100`,
  matched by the c60-size-override run whose remat moves did not include the g·C.
