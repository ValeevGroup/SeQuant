# External batch loops in the multilevel-cache scope model: design spec

**Status:** design. Revised 2026-07-28 (2nd re-pin, after verbatim trace analysis
of summand 46). The C60 peak root cause is **NOT** a cache-hoist. That earlier
framing -- and the `consumed_upscope` emit + hoist-exclusion built on it (SeQuant
Tasks 1-3, commits `e66fdb6a4`/`0711a1586`/`e3174a468`) -- are a **correct
latent-bug fix but ORTHOGONAL to C60**: a controlled probe showed the fix
byte-inert on this witness (the C60 peak is the batched-*scratch*, not the cache,
which is only ~94 GB < budget). The real cause (Sec. 1, verbatim-grounded): a
**shared array used in more than one batch context** is folded by the runtime
cache onto one entry, stored full at its outermost use scope (correct, and
peak-acceptable), and then served **whole** to a *nested external* use -- because
the fallback `access` does not slice the fetched array to the accessing scope's
block. `le_g` slices only leaves, never a cached intermediate. The missing model
element is **use-requires-slicing (slice-on-use)** (Sec. 2.5). `consumed_upscope`
(Sec. 2.3) remains a correct sub-model but is not the C60 lever. The general
CSE x multimode-batching interaction (how a shared array's single lifetime scope
is chosen across all its use contexts, including cross-term) is an OPEN model
question under active discussion (Sec. 9), not settled here.

Amends the multilevel-cache / nested-batching lifetime model of
`mpqc4/doc/dev/specs/2026-07-20-order-aware-multilevel-batching-design.md`
(the `scope_level`/scope-chain caching contract, implemented as SeQuant Tasks
A4/B1-B4, commits `0c36ebe24`..`235309f01`) and the external-mode batching of
`SeQuant/doc/dev/specs/2026-07-20-external-mode-batching-design.md`. Those two
were designed separately; this spec makes external batch loops first-class
members of the cache-scope hierarchy, which neither currently does.

**Repos:** SeQuant (`core/optimize/cost_model.hpp`, `core/eval/eval.hpp`,
`core/eval/cache_manager.hpp`, `core/eval/eval_expr.{hpp,cpp}`). Consumed by MPQC
CCk (CSV/PNO-CCSD).

**Reproducer / evidence:** the C60 pVDZ-F12 PNO-CCSD residual replayed through the
real batched runtime (`tests/unit/test_eval_dryrun.cpp` `[.][dryrun-occ-veto]`),
role-separated external-occ config (occ in the EXTERNAL role only:
`is_batchable_external_index = {occ i}`, `is_batchable_contracted_index = {Κ}`,
`μ̃`/pao contracted or unbatched, occ block 8, `peak_threshold = 100 GB`,
`order_aware_recompute = true`). Post role-split (SeQuant `79540c831`..`7f5240014`,
MPQC `1685c7ac1e`..`5a3d94dc63`) the honest replay of `[.][dryrun-occ-veto]`
(nterms=55) measures peak = **5860.877 GB**, avoidable_time = **44.6%**,
`contracted_occ_stamps == 0`, `external_occ_stamps == 244`. The peak-replay term
(summand 46) is modeled at **35.08 GB** -- a **167x** replay-vs-model gap FOR THAT
TERM. **NB: 35.08 GB is summand 46's per-term model, NOT the whole-forest target**
(an early framing error corrected below). The forest modeled peak is the MAX over
summands, dominated by summand 12 at ~445 GB.

> **RESULT (2026-07-28): slice-on-use LANDED (SeQuant `ac4f62377`, phase 1) and
> did its job.** Re-measured `[.][dryrun-occ-veto]` (nterms=55): peak
> **5860.877 -> 443.55 GB (13.2x drop)**. The peak is now summand 12, whose
> **replay (443.55 GB) MATCHES its factorizer model (445.46 GB)** -- the 167x
> replay-vs-model gap is ELIMINATED; the runtime now realizes what the DP predicts.
> Summand 46's giants are sliced (it is no longer the peak). Avoidable_time
> 44.6 -> 42.67% (barely changed -- slice-on-use fixes the PEAK/residency, not the
> recompute/CSE). The remaining 443 GB is the DP's honest schedule cost for summand
> 12; reducing it (toward the aspirational `< 100 GB`) is a factorizer/cost-model
> job = phase 2, not slice-on-use. Real-MPQC batched-energy correctness is not yet
> validated (OFF-path is byte-identical, synthetic slice-exactness passes).

---

## 1. Problem

Batching a **contracted** index and batching an **external** index are not
symmetric in the cache model, and the current model was written for the
contracted case only.

- A **contracted** index is summed away: the batched node has exactly ONE
  result (the accumulator), so it needs exactly ONE cache entry, at its
  lifetime-scope. The multilevel-cache model assigns each node one scope and
  one cache. This works.
- An **external** index is free -- carried onto the result. Batching it is
  inherently **dual**: for a node `R{i}` batched over external `i`, there is the
  per-block slice `R[I]{i}` (computed inside the loop) AND the full assembled
  result `R{i}` (all `i`, the scatter destination). Because the cache keys nodes
  by canonical identity with NO slice index, `R[I]` and `R{i}` collide on the
  same key -- the single-entry model cannot hold both.

### 1.1 Root cause (RE-PINNED, verbatim): a shared array used in two batch contexts, served whole to the nested one

The C60 peak is the batched-**scratch** high-watermark (two ~2930 GB giants
co-resident), not the cache (~94 GB, under budget). The giants are full because a
**single physical array is used in two different batch contexts within one term**,
folded onto one runtime-cache entry, and served **whole** to the context that
needs it sliced.

**Verbatim evidence (`[.][dryrun-occ-veto]`, nterms=55, summand 46 -- the peak).**
The full term (stock `SEQUANT_DUMP_PEAKTERM`, `to_latex`):

```
2 · g^{i3}_{i1}[K2] · g^{i4}_{μ̃1369}[K2] · C^{μ̃1369}_{a3<i2>}
  · s^{μ̃1371}_{μ̃1370} · C^{μ̃1370}_{a2<i1,i2>} · C^{a5<i3,i4>}_{μ̃1371}
  · s^{μ̃1373}_{μ̃1372} · C^{μ̃1372}_{a4<i3,i4>} · C^{a1<i1,i2>}_{μ̃1373}
  · t^{i2}_{a3<i2>} · t^{i3 i4}_{a4<i3,i4>,a5<i3,i4>}
```

The residual is `R = I(a1<i1,i2>; a2<i1,i2>)` -- so `i1,i2` are **external** (the
residual pair) and `i3,i4` are **contracted** (summed, absent from `R`). The
intermediate that drives the peak is the overlap x PNO-coefficient product
`s · C = I(μ̃; a<occ-pair>)` (NOT `g · C`). It appears at **two distinct tree
nodes** of the same canonical hash `1989507463377952644`:

| node | verbatim | occ pair | role | `batched_here` | `scope_level` | consumer |
|---|---|---|---|---|---|---|
| #8  | `I(μ̃1372; a1<i1,i2>)` | (i1,i2) | **external** | `{i1:EXT i2:EXT}` | 1  | giant #7 `I(a1;a4)` |
| #10 | `I(μ̃1370; a5<i3,i4>)` | (i3,i4) | **contracted** | `{}` | -1 | giant #9 `I(a5;a2)` |

Both are the **same physical array** (canonicalization relabels the bound occ-pair
domain; the `s` and `C` children are hash-identical). CSE *in factorization* is
off, so the optimizer emits them as two independent subtrees, but the **runtime
cache** dedups them by canonical identity onto one entry. #10's occ pair is
contracted, so it is (correctly) unbatched and stored **full** at the outermost
use scope (`scope_level = -1`, 8.7 GB -- peak-acceptable, because it *is* used
unbatched there). #8's occ pair is external, so its use must be **sliced** to the
current `i1,i2` block -- but #8 reads the array through the scope-chain fallback
`access`, which returns it **whole** (`le_g` slices only leaf results, never a
cached intermediate). So the external consumer contracts the full 8.7 GB `s · C`
into the **full 2930 GB giant**; two such half-giants (`#2` and `#9`, contracted
over `a5<i3,i4>` to form `R`) are co-resident = **5860.877 GB**. The factorizer
models the external use sliced -> 35.08 GB; the 167x gap is exactly the unrealized
slice-on-use.

### 1.2 What is and is not the cause

- **NOT a cache-hoist / not `hoist_invariants`.** A `SEQUANT_UT_PEAK_COMPONENTS`
  probe split the peak: scratch = 5860.877 GB, cache.hwmark = 94.018 GB, and
  neutralizing Tasks 1-3's `!ext_loop_local` conjunct left **both byte-identical**.
  The `consumed_upscope` hoist-exclusion cannot touch #10 anyway -- #10 carries the
  *contracted* pair (i3,i4), which is not in the external mask, so it never has an
  External stamp to key on.
- **NOT a canonicalization / correctness bug.** #8 and #10 are the same abstract
  tensor; the shared hash is correct. The folding serves a correct value -- the
  defect is purely that the external use is not sliced.
- **NOT `scope_level == -1` per se.** The array's outermost-use lifetime scope
  (full at `{}`) is *correct* -- that is where it is peak-acceptable. The bug is
  the absent slice on the *nested* use.

**Root cause, one line:** the model deduces a shared array's lifetime scope from
its outermost use (via the scope-chain fallback) but does **not** slice the array
when it is *used* in a nested batch scope that batches some of its modes -- so a
full-at-outer-scope array is served whole into an inner external block that
required its slice.

---

## 2. Design

### 2.1 External loops are members of the scope hierarchy

The lifetime-scope of a node is the innermost enclosing **batch loop** -- whether
**external or contracted** -- among the batched indices the node's result carries.
Concretely, `scope_level` for node `n` must be computed over the full realized
loop nest enclosing `n` (external loops included), not only over `cell_seq_`'s
contracted placement.

- A node carrying a batched external index `i` is scoped to the loop of `i`
  (or an inner loop it also carries), NEVER run/term scope. `scope_level == -1`
  is reserved for a node that carries NO batched index of any kind.
- `cell_seq_` stays contracted-only (its job is contracted *nesting order*, which
  is meaningless for an external mode that is never a contracted-here set). The
  external contribution to `scope_level` comes from the node's carried external
  modes and the external loop nest the emit already placed on its ancestors --
  a separate input to the scope computation, unioned with the `cell_seq_`
  placement.

This is supporting machinery: it makes the cache chain place `R[I]` (loop scope)
and `R{i}` (enclosing scope) at distinct levels. It does NOT by itself stop the
full hoist of a loop-local intermediate -- that is Sec. 2.3.

### 2.2 Slice-vs-full is a CACHE-LEVEL property, not a key property

This is the load-bearing idea for disambiguating the dual entries. `R[I]{i}`
(slice) and `R{i}` (full) share one canonical cache key; they are disambiguated
purely by WHICH scope's cache holds them:

- `R{i}` (full, all `i`) lives in the cache of the scope **enclosing** the
  external `i`-loop.
- `R[I]{i}` (the current block) lives in the cache of the **loop** scope.
- `access(R)` from INSIDE the loop hits the loop-scope cache first and returns
  `R[I]`; from OUTSIDE it falls through the scope chain to the enclosing cache
  and returns `R{i}`. Same key, correct value, no slice index in the identity.

The scope chain (`CacheManager::parent_` / `set_parent` / fall-through `access`,
`cache_manager.hpp:166,205,257`) is exactly this mechanism; the model just has to
place external loops into it.

### 2.2b Lifetime = outermost use scope; a nested use must slice (slice-on-use) -- THE C60 LEVER

Sec. 2.2 assumed the loop-scope cache already HOLDS the slice `R[I]`. That is true
for a scatter target (Sec. 2.4 writes `R[I]` per block). It is NOT true for a
**shared array used in more than one batch context** -- the C60 case (Sec. 1.1),
where `s*C` is used with its occ pair external at one consumer and contracted at
another. The complete model has two parts, and only the first was in the spec:

- **Lifetime = outermost use scope.** A physical array (one canonical identity,
  possibly reached from several tree positions and several terms) is assigned ONE
  lifetime scope: the OUTERMOST (least-nested) of all its use scopes. If it is used
  anywhere unbatched (scope `{}`), it lives full there, and that full size is
  peak-acceptable *by construction* -- an unbatched use means the placer did not
  need to batch it (its full footprint was under budget on that path). The
  scope-chain fallback already effects this: the outermost-scope copy is the one
  that survives resets and is served to inner scopes.
- **A use requires slicing (slice-on-use) -- MISSING.** When that array `A` (living
  at scope `S`) is USED inside a nested evaluation scope `S + {I}` that batches a
  mode `I` which `A` carries, the use must fetch `A` and take the slice
  `A[I-block]` for the current block, NOT consume `A` whole. Today the fallback
  `access` returns `A` whole and `le_g` slices only LEAF results, so a cached
  intermediate carrying a batched mode escapes the slice entirely. `A` living full
  at the outer scope is correct and cheap (8.7 GB for `s*C`); the fix is to slice
  it on the inner external use so the consumer gets `A[i1,i2-block]` (a ~13 GB
  giant) instead of contracting full `A` (the 2930 GB giant). Slice-on-use is the
  cached-intermediate analogue of `le_g`.

This is the load-bearing fix for C60. It subsumes the Sec. 2.3 dichotomy: a
consumed-upscope array (used outside its loop) lives full at the outer scope and is
sliced on inner uses; a strictly loop-local array (used only in one batched scope)
lives sliced at that loop scope with no outer full. The general rule -- **lifetime
= outermost use, slice on every nested batched use** -- covers both, and covers the
mixed external/contracted-context array that neither Sec. 2.3 nor the scatter
handled.

### 2.2c The access mechanism: `materialize` + `slice_to_use`

Slice-on-use lives at the point a node's value enters the current compute context.
Exactly two code paths return a directly-available value today -- the cache hit
(`eval.hpp:619`) and the leaf eval (`:659`) -- and they unify. Terminology:
**scope** = a scope ID (the axis set `{I1,I2,...}`); **scope_level** = a scope
depth (an int); a node's **use scope** is the axis set active at its depth, its
**lifetime scope** is the (prefix) scope where its value actually lives.

- **`leaf_evaluator`** (the `le` rename) -- user-provided RAW leaf compute:
  evaluate a leaf to its FULL value, NO slicing.
- **`materialize(node) -> {value, lifetime_scope}`** -- the unified accessor: the
  cache hit if cached, else the `leaf_evaluator` value if a leaf, else null (an
  internal not-yet-computed node -- produced by the iterative eval loop, not here).
  It surfaces the **observed** lifetime scope: where the value resolved (top/full
  for a freshly materialized leaf; the hit level for a cached intermediate).
- **`batch_context`** -- the ordered stack of enclosing `{axis, block}`, one entry
  per realized loop, outermost-first (what `le_g` closed over, made explicit).
- **`slice_to_use(value, node, lifetime_scope)`** -- for each `{axis, block}` in
  `batch_context` BELOW the lifetime scope that `node` carries, apply the existing
  `slice_mode`. The slice set is exactly (use scope MINUS lifetime scope) intersect
  the node's carried modes.

The Enter stage collapses the hit and leaf branches into one, sliced once:

```
if (auto m = materialize(node); m)
  return finalize(slice_to_use(apply_phase(node, m.value), node, m.lifetime_scope));
// null: internal-not-cached -> custom-evaluator interception / push children
```

Because `slice_to_use` reads the **observed** lifetime scope -- where the value
actually resolved -- not an emitted per-node bit, it is robust to the
inconsistent-stamp problem of Sec. 1.1: the two positions of a shared array fold
onto one cached value, and every use slices that one value from wherever it truly
lives. `le_g` is deleted (its eval half is `leaf_evaluator`, its slice half is
`slice_to_use`). `materialize` is the synchronous accessor only -- computing a
missing internal node stays the iterative eval loop's job (stack-safe), and a
caller that wants to force a full subtree build still uses `evaluate(node, ...)`.

**Fate of the emitted `scope_level` metadata (investigated).** Its ONLY readers
are the run-scope veto (`cache_manager.hpp:622`, `batch_scope_level() >= 0` => not
run-cacheable) and `hoist_invariants` (`eval.hpp:1492,1515`). Both are the
order-aware *hoisting* optimization (build a loop-invariant once at an outer scope
so a per-iteration `reset()` does not force its recompute) -- NOT correctness.
Under `materialize` + `slice_to_use`, correctness derives from the observed
lifetime scope, so the emitted `scope_level` (and `consumed_upscope`, and the
veto) degrade to hoisting perf hints. Whether hoisting is still worth keeping --
or whether store-at-computation-scope + slice-on-use subsumes the whole Tasks
1-3 / A4 `scope_level`/`consumed_upscope`/`hoist_invariants` apparatus -- is an
open question (Sec. 9). Do NOT assume that apparatus survives unchanged.

### 2.2a The scope hierarchy IS the loop nest: one cache per loop level

"Scope = innermost carried loop" means there is one cache per loop level, i.e.
per scope (each loop is a scope with its own per-iteration reset). For the mixed
`[e outer, c inner]` nest of 2.4a that is THREE cache levels: outermost, the
e-scope, and the (e-then-c)-scope. This is not new machinery -- it is the
existing per-loop scratch structure:

- The batched evaluator already creates one scratch `CacheManager` per batch loop
  (`make_batched_scratch`, reset once per iteration); a nested loop reinstalls a
  fresh inner scratch. The committed `parent_`/`set_parent`/fall-through `access`
  (`cache_manager.hpp:166,205,257`) is exactly the link between these per-loop
  caches. This spec wires that chain and assigns each node to the correct level;
  it adds no caches.
- **The count is nest DEPTH, not combinatorial.** At any instant the runtime is
  on one path (e-block `E`, c-block `C`), so exactly the caches on that path are
  live: outermost, the e-scratch for `E`, the c-scratch for `C`. Advancing `e`
  scatters `R[E]` out and recreates the inner scratches. You never hold
  `nBatch(e) * nBatch(c)` caches -- it is a chain of length = loop-nest depth
  (bounded by the `depth < 8` runtime backstop). Where things live: full `R{e}`
  at the outermost, `R[E]` accumulator at the e-scope, transient partials
  `A{e,c}` at the (e,c)-scope.

### 2.3 A full entry exists ONLY for a node consumed up-scope (landed; NOT the C60 lever)

> **Status note (2nd re-pin).** This subsection drove SeQuant Tasks 1-3
> (`consumed_upscope` emit + hoist-exclusion). Those LANDED and are a correct
> latent-bug fix, but they are **orthogonal to the C60 peak** (Sec. 1.2): the C60
> lever is slice-on-use (Sec. 2.2b), which this subsection did not contain. Keep
> this as the sub-model for the loop-local single-context case; do not expect it to
> move the C60 peak.

**Terminology.** A node is *consumed up-scope* when its result is read from an
enclosing scope -- one level UP the scope chain (Sec. 2.2), i.e. OUTSIDE its own
external loop. "Up-scope" and "outside its loop" are the same condition; the name
ties it to where the full entry lives (the enclosing-scope cache).

Peak correctness hinges on this. The full `R{i}` entry is created **only** for a
node whose result is read up-scope (outside its external loop) -- the genuine
outputs:

- the scatter TARGET (a term result / residual contribution carrying `i`), and
- any shared final also consumed by a different loop nest (handled by the
  external forest-merge of Sec. 8; for a single term in isolation this reduces to
  the scatter target).

A **pure intermediate** consumed entirely WITHIN the loop (e.g. the summand-46
operand, contracted into the next factor in the same block) gets ONLY the
loop-scope block slice -- **no full entry at all.** Its block contribution flows
up the contraction chain into the target's block `R[I]`, which is what the scatter
writes into `R{i}`.

Giving every external-carrying intermediate a full entry re-materializes the full
intermediate by construction -- the exact 5860 GB failure. So the rule is: full
entry iff consumed-upscope; loop-scope-only otherwise.

**Classification location (DECIDED): the DP emits a `consumed_upscope` bit.**
During the reconstruct/emit walk the DP already builds `NodeBatchAnnotation` per
node (carrying `scope_level`, `effective_count`, batch stamps); it adds a boolean
`consumed_upscope`. Per term the rule is local and exact:

> a node is `consumed_upscope` iff it is the **outermost carrier** of a batched
> external mode -- i.e. its parent does not carry that mode.

Because an external mode is free and propagates up the contraction chain to the
term result, the outermost carrier is the term's external output (the residual
contribution `R{i}`, the scatter target). Every external-carrying **intermediate**
has a parent that also carries the mode, so it is `consumed_upscope = false`
(loop-local). No use-count reasoning is needed for the single-term case; the
tree structure decides it. (Cross-term shared finals are Sec. 8.)

### 2.4a Mixed external+contracted nests: scope is the innermost CARRIED loop

A (sub)tree may be wrapped in a nest of several loops mixing external and
contracted modes -- e.g. `R{e} = sum_c A{e,c} B{c}` batched over external `e`
(outer) AND contracted `c` (inner). This composes without any `R[E,C]` entry,
because a node's scope is its innermost **carried** loop (the innermost loop
whose index appears on the node's result), NOT the innermost loop of the nest:

- `R` carries `e` but NOT `c` (`c` is summed away). So `R`'s scope is the
  **e-loop**, above the inner c-loop. `R` lives in two caches: full `R{e}` at the
  enclosing scope and `R[E]` at the e-loop scope. There is no `R[E,C]` -- `R`
  never carries `c`.
- The inner contracted c-loop, which `R` is invariant to, is realized by
  **accumulation** into `R[E]`: `R[E]` lives at the e-loop scope, ABOVE the
  c-loop per-block reset line, so the c resets do not clear it (the accumulator
  pattern of the contracted model). Each c-block adds its contribution into
  `R[E]`.
- The outer external e-loop is realized by **scatter**: `R[E]` is written into
  `R{e}`'s e-slice in the enclosing cache.

So `R` is simultaneously a contracted accumulator (over the inner `c` it does not
carry -> accumulate into `R[E]`@e-scope) and an external scatter source (over the
`e` it carries -> scatter `R[E]` into `R{e}`@enclosing). A node carrying BOTH `e`
and `c` (an un-summed partial like `A{e,c}`) scopes to the innermost carried loop
(`c`) and is loop-local/transient; a full-over-`c` entry can never arise, because
a contracted index has no full extent on any result -- so the full/dual entry of
2.3 is only ever over the EXTERNAL modes a node carries.

**Consequence:** external and contracted CAN be co-batched in one nest. The load-
bearing constraint is that the scope rule is innermost-CARRIED (never
innermost-of-nest); placing a node at the innermost loop of the nest would send
it looking for a nonexistent `R[E,C]` (the failure mode this rule avoids).

### 2.4 The scatter IS the enclosing-cache write

`pre_sized_zeros_over_mode` + `write_into_slice` (`eval.hpp` external branch,
~`:1564-1636`) is the realization of "assemble `R{i}` in the enclosing scope by
writing each block's `R[I]` into its `i`-slice." It fires ONLY for the dual
(consumed-upscope) nodes of 2.3; it must NOT be applied to a loop-local
intermediate (that is how the full intermediate gets built today).

---

## 3. Correctness invariants

- **Exactness.** `R{i} = concat_I R[I]{i}` over disjoint external blocks; the
  external partition is work-neutral (no cross-block coupling, flops unchanged).
  Serving `R[I]` inside the loop and `R{i}` outside is slice-exact by
  construction (distinct caches, distinct scopes).
- **No full intermediate.** No node whose result carries a batched external index
  and is consumed only within its loop may have a full-extent entry anywhere. The
  batched-scratch peak must reflect only loop-local block slices plus the genuine
  outputs' full results.
- **Contracted unchanged.** Contracted batching keeps its single-entry model;
  this spec adds external as a parallel case, and reduces to the current behavior
  when no external mode is batched.
- **OFF-path byte-identical.** With `order_aware_recompute = false` /
  `batch_spectator_indices = false`, no external loop is emitted, every node keeps
  the sentinel scope, `consumed_upscope` is never consulted, and nothing here
  engages -- byte-identical to today.

---

## 4. Interaction with existing machinery

- **`scope_level` emit (A4, `cost_model.hpp:1794-1809`).** Extend to union the
  external enclosing loops the node carries into the placement, so an
  external-only-enclosed node scopes to its external loop, not `-1`. `cell_seq_`
  and the contracted placement are unchanged. (Supporting; see 1.2 -- necessary
  but not sufficient.)
- **`consumed_upscope` emit (new).** Set on `NodeBatchAnnotation` during the same
  emit walk (2.3): outermost carrier of an external mode -> true, else false.
- **The hoist store (B2, `eval.hpp` `hoist_invariants`).** Today it stores a
  non-volatile node at its `scope_level` cache. Add the gate: a node carrying a
  batched external index and `consumed_upscope == false` is **excluded from the
  real-cache hoist** -- it falls through to the per-block loop scratch (reset each
  block) and is rebuilt sliced per block. This is the fix that the A/B in 1.2
  showed `scope_level` correction alone does not deliver.
- **The scatter (`eval.hpp` external branch).** Restrict `pre_sized_zeros_over_mode`
  + `write_into_slice` to `consumed_upscope` nodes (2.3/2.4); a loop-local
  intermediate must not get a full destination.
- **The veto (B4, `cache_manager.hpp` batch-variant veto, commit `235309f01`).**
  `batch_variant = sliced_batch_axis || scope_level >= 0` already refuses
  run-scope residence for a node with `scope_level >= 0`. With external in the
  scope (2.1) it keeps external-carrying intermediates out of the outermost cache;
  the hoist-exclusion above is what additionally keeps a loop-local one out of its
  own `scope_level` real cache.

---

## 5. Implementation shape (for the follow-on plan; not this spec)

**0. Slice-on-use (Sec. 2.2c) -- the C60-critical piece.** Designed (Sec. 2.2c):

- **Rename `le -> leaf_evaluator`** across the `evaluate` overloads and
  `make_batched_custom_evaluator` (own mechanical commit; behavior-neutral).
- **Make `access` surface the resolved level.** The scope-chain `access`
  (`cache_manager.hpp`) returns not just the value but the chain level where it
  hit, so the caller knows the value's lifetime scope.
- **Thread `batch_context`** (the ordered `{axis, block}` stack) through the
  batched custom evaluator (it currently lives implicitly in the composed `le_g`
  closures; make it explicit) so the value path can compute (use MINUS lifetime).
- **Add `materialize` + `slice_to_use`** and collapse the `evaluate` Enter stage's
  hit (`:619`) and leaf (`:659`) branches into the single sliced access of
  Sec. 2.2c. Delete `le_g`.
- **Validate on the C60 gate** (Sec. 6.3): the summand-46 `s*C` is served sliced,
  the giants drop to block extent, peak 5860.877 -> 443.55 GB (the replay now
  MATCHES the DP model; summand 46 sliced, the forest peak becomes summand 12 at
  ~445 GB -- the DP's honest cost), slice-exact. (Achieved, SeQuant `ac4f62377`.)

Items 1-4 below (`consumed_upscope` / `scope_level` / hoist-gate, landed as Tasks
1-3) do NOT deliver this and do not move the C60 peak; under slice-on-use they
degrade to hoisting perf hints whose necessity is open (Sec. 9). Do not extend
them expecting a C60 win.

1. **DP emit `consumed_upscope`.** In the reconstruct/emit walk, set the bit on
   each node: outermost carrier of a batched external mode -> true, else false
   (2.3). Thread it through `NodeBatchAnnotation`.
2. **Scope-level for external (supporting).** Union the node's carried external
   batched modes + the external loop nest into the `scope_level` computation so
   external-enclosed nodes get `scope_level >= 0` (2.1).
3. **Runtime hoist gate.** In `hoist_invariants`, exclude from the real-cache
   hoist any node carrying a batched external mode with `consumed_upscope == false`
   -> it lands in the per-block loop scratch and is rebuilt sliced per block (4).
4. **Runtime scatter gate.** Restrict `pre_sized_zeros_over_mode` +
   `write_into_slice` (full destination) to `consumed_upscope` nodes (2.4).
5. Rely on the existing veto + scope chain for run-scope exclusion and full/slice
   disambiguation.

---

## 6. Validation

1. **Unit (emit).** A small network with an external index `i` free on the result
   and an intermediate `M` carrying `i` under the term: assert `M`'s emitted
   `consumed_upscope == false` and the output carrying `i` is `consumed_upscope ==
   true`; assert `M`'s emitted `scope_level` is the `i`-loop level (>= 0).
2. **Unit (runtime).** Replay that network: `M` is realized at block-`i` extent
   (not full) -- assert no full-extent cache entry for `M`; the output `R{i}` is
   assembled full via the scatter; result equals the unbatched reference
   (slice-exact).
3. **C60 gate.** On `[.][dryrun-occ-veto]` with external-occ on, the summand-46
   operand/giant is realized at block occ, and the batched-scratch peak drops from
   5860.877 -> 443.55 GB, matching the per-term DP model (summand 12 replay
   443.55 ~= model 445.46); NOT the ~35 GB (that was summand 46's per-term model,
   not the forest peak); converged energy unchanged.
4. **OFF-path parity** -- byte-identical with external/order-aware off.

---

## 7. Open questions

- **Multiple external modes (2-nest `i₁,i₂`).** Each adds a scope level; a
  loop-local intermediate scopes to the innermost external loop it carries. The
  full/slice dual composes (full only at the outputs); confirm the scatter's
  nested `pre_sized`/`write_into_slice` writes into the right enclosing level for
  each of the two loops. (The summand-46 term is exactly this 2-nest.)
- **`consumed_upscope` and the emit walk order.** The bit is "parent does not
  carry this external mode." Confirm the emit walk has the parent's carried-mode
  set available when it stamps the child (top-down), or compute it in a short
  post-pass over the annotated tree.

---

## 8. In-scope follow-on (b): cross-term CSE via external forest-merge

**Why this is in scope, not an aside.** Fix (a) confines each term's external
intermediates to loop-scope block slices, rebuilt per block. But the batched
evaluator's cross-term CSE machinery -- the `BatchGroup` that co-evaluates
persistent finals shared across terms so a shared sub-intermediate is built once
per batch instead of once per consumer -- is **contracted-only**. The external
scatter is explicitly single-node ("an external mode is not a persistent-final
sharing mode, so the group/replay machinery ... does not apply -- scatter just
this node", `eval.hpp:1585-1589`, `solo{{&node, K}}`). So once we batch over an
external mode, a final shared across terms is **rebuilt per term** -- the CSE the
contracted `BatchGroup` (`eval.hpp:1654-1673`, joins persistent finals realizing
the identical batch partition) would have shared is lost.

For the C60 residual this is real: shared external-carrying finals (e.g. the
half-transformed DF factors reused across summands) recompute per term under (a).
So (a) makes the peak correct but leaves per-term recompute on the table.

**(b) = external forest-merge across terms.** Give external batching the analogue
of the contracted `BatchGroup`: collect, across terms, the `consumed_upscope`
external finals that stream over the **same external partition**, and co-evaluate
them in one shared external loop nest (their shared sub-intermediates built once
per block). This is the external counterpart of the group-collection at
`eval.hpp:1654-1673` and the layered co-evaluation that follows it.

**Sequencing.** (a) first (this spec) -- it fixes the peak and is a precondition
(the loop-scope block-slice discipline is what external group members co-evaluate
over). (b) second, its own spec/plan -- it recovers the cross-term CSE that (a)
alone drops. Both are in scope for the overall external-batching effort.

---

## 9. OPEN: how CSE and multimode batching interoperate (needs a settled model)

**What the code actually does today** (verified, `[.][dryrun-occ-veto]`):

- **Two independent CSE axes.** `OptimizeOptions::CSE.subnet` controls *factorization*
  CSE (merging identical subnetworks in the binarized tree); it is OFF in this run.
  Separately, the **runtime cache** (`CacheManager`, keyed by `hash::value` +
  connectivity-graph equality, `eval_node_compare.hpp`) dedups by canonical
  identity regardless of that flag. So identical subexpressions ARE shared at
  runtime even with factorization CSE off.
- **Cross-term sharing is ON.** `cost_profile` (`cost_profile.hpp:263-277`) replays
  every summand through ONE shared cache, calling `reset()` between summands, which
  drops non-persistent scratch but **keeps persistent (cross-term) entries**
  (`:73-75`). A `min_repeats`-persistent intermediate is therefore shared across
  summands. (The C60 blow-up is intra-summand -- #8/#10 in summand 46 -- but the
  same array also lives cross-term.)

**The model tension this exposes.** One physical array (one canonical identity) can
be *used* in several batch contexts (loop nests) -- e.g. `s*C` with its occ pair
external at one use and contracted at another, and potentially in different terms.
The current lifetime model is "one lifetime scope (one cache level) per array,"
deduced from uses via the scope-chain fallback (Sec. 2.2). Sec. 2.2b adds the
missing half (slice-on-use). But several questions must be settled before the
model is trustworthy:

1. **Lifetime-scope choice across ALL uses.** Sec. 2.2b says "outermost use scope."
   Confirm this is well-defined and correct when uses span terms and mix
   external/contracted roles: is it always the least-nested enclosing loop common
   to all uses, and is "unbatched somewhere => full is peak-acceptable" always
   true (or can a cross-term use force a full entry that a single term's budget
   never sanctioned)?
2. **Slice-on-use vs. store-the-slice.** Two realizations give the same result:
   (i) store `A` full at the outer scope and slice on each inner use (Sec. 2.2b),
   or (ii) store the per-block slice `A[I]` at the loop scope (Sec. 2.2/2.4). When
   is each preferred -- e.g. does a cross-term-shared `A` want the full outer entry
   (shared once) plus slice-on-use, while a single-term-loop-local `A` wants only
   the block slice? Does mixing them per-array stay coherent under the cache's
   single-entry-per-identity constraint?
3. **Recompute vs. residency.** Slice-on-use keeps `A` resident full at the outer
   scope (cheap here, 8.7 GB) and re-slices per block (cheap). The (b) forest-merge
   instead co-evaluates shared finals per block. These are different points on the
   recompute/residency trade; the cost model (`n_replay`, persistence weighting)
   must price whichever the runtime actually does, or the DP's schedule and the
   realized peak diverge (as they did here, 35 vs 5860 GB).
4. **Does the DP even know?** The emit annotates per tree position
   (`node_batch_axes`), but the runtime shares per canonical identity. A single
   identity reached from positions with different batch roles (external #8 vs
   contracted #10) has conflicting per-position annotations. The model must define
   which annotation governs the shared entry -- or make the runtime slice-on-use
   independent of the per-position stamp (preferred: slice from the *use site's*
   active block stack, not from an emitted bit).

**Resolution so far (design discussion, 2026-07-28).** The `materialize` +
`slice_to_use` model (Sec. 2.2c) settles Q4 and most of Q2:
- **Q4 (which annotation governs a shared entry):** none -- the runtime slices
  from the **observed** lifetime scope (`materialize` surfaces the resolved chain
  level) and the current `batch_context`, independent of any per-position emitted
  stamp. The two positions of a shared array fold onto one cached value; every use
  slices it correctly regardless of how it was stamped.
- **Q2 (slice-on-use vs. store-the-slice):** not a free choice. Cache lifetimes
  are: the top cache is forest-persistent, every nested scratch is
  per-iteration-transient (verified, Sec. 2.2a). So **cross-term-shared => lives
  full at the top => slice-on-use is the only coherent realization**;
  store-the-slice at a loop scope is available only for a strictly loop-local,
  non-shared array. Persistence forces the choice.
- **Q1 (outermost-use lifetime cross-term):** still to confirm -- is the full
  entry always <= some single term's sanctioned budget, or can cross-term sharing
  force a full residency no single term's `place` approved? (Likely fine because a
  persistent full entry is 8.7 GB-class, but state the bound.)
- **Q3 (recompute vs. residency) and the `scope_level`/hoist apparatus:** OPEN.
  `scope_level`'s only readers are the veto + `hoist_invariants` (Sec. 2.2c),
  both hoisting perf, not correctness. The clean question for the follow-on: once
  slice-on-use lands, is store-at-computation-scope + slice-on-use enough, letting
  us **retire** the `scope_level` emit, `consumed_upscope`, the run-scope veto, and
  `hoist_invariants` -- or is per-iteration invariant hoisting still a needed
  optimization the cost model must price? Do not treat the Tasks 1-3 / A4 apparatus
  as load-bearing; slice-on-use may subsume it.

## 10. Lifetime = per-slot mask (model, validated on the real IR)

The lifetime scope of a canonical node is a **per-canonical-slot bitmask**
(each mode SLICED or FULL), computed by a MEET over the node's occurrences: a slot
is SLICED iff every occurrence batches it to the same block (proto-aware: a
composite PNO index `a<i,j>` is sliced when its proto occ pair is batched), FULL if
any occurrence is unbatched over it or occurrences disagree on the block. Placement:

- **all-full** -> the top cache (canonical, block-agnostic, full). E.g. `s*C`.
- **all-sliced** -> the deepest loop scratch (fully block-local, transient).
- **mixed** -> the scratch of the **deepest sliced mode's loop**, stored full over
  the full modes -- the *hoisting* case (built once at that scope, reused across
  inner loops it is invariant to).

CSE operates at the placement scope: every consumer at/below it hits the shared
entry; the scope's reset boundary ends the sharing. This replaces the veto + hoist
+ `consumed_upscope` machinery with placement-by-rule.

**Cost model must be lifetime-aware.** The per-slot mask sets each intermediate's
residency (full vs sliced), hence the peak, hence the cost -- so the DP must run
the SAME meet to size intermediates (this is exactly why the DP modeled summand 46
at 35 GB while replay hit 5860: the DP assumed `s*C` sliced, the runtime kept it
full). A transposed use (`M(y,x)` binding the external axis to a different slot) or
an unlucky loop order forces a slot FULL; the cost model must charge that, not
optimistically assume sliced. The per-slot meet is one **shared analysis** feeding
both the cost model (sizing/peak/order) and the runtime (placement/slicing).

**Validated (IR survey, `[.][dryrun-occ-veto]`, 224 canonical nodes,
`lifetime-mask-ir-survey.md`).** Distribution: 172 all-full, 15 all-sliced, 37
mixed. The hoisting shape occurs (occ-pair sliced, aux `K` full -> occ scratch,
hoisted out of the K loop); `s*C` is all-full and currently re-sliced/recomputed in
24/52 terms (the recompute slice-on-use removes). Findings folded in above:
(i) the meet must be **proto-aware** -- and this is a *correctness* matter (a
mis-marked node lands in the wrong tier), not just sizing; (ii) the mask **tightens
top-cache membership** (only all-full is block-agnostic-cacheable), which resolves
the pair-specific-data hazard the current persistence heuristic risks; (iii) the
masks are dominated by **within-occurrence** structure (the cross-occurrence
demotion never fires here), so the analysis is per-occurrence-cheap.

**Sequencing.** Phase 1 = slice-on-use (Sec. 2.2c) -- independent of the mask
(slices from the OBSERVED lifetime), fixes the C60 gate. Phase 2 = mask-based
placement + cost-model lifetime-awareness + retire veto/hoist/`consumed_upscope` --
its own spec+plan after phase 1 lands.
