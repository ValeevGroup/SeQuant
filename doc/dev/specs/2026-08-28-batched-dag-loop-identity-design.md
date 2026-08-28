# Batched forest → DAG → schedule → execution: loop identity and the per-occurrence mode↔loop atlas

**Status:** design draft — 2026-08-28
**Supersedes:** `2026-08-27-frame-correct-use-induced-slicing-design.md` (which
conflated `ordinal`; see §11).
**Area:** `SeQuant/core/eval` (`compute_dag_boulevard`, `analyze_legality`,
`build_ordered_schedule`, `dag_scope.hpp`, `value_id.hpp`, the runtime slice
seam) and `placement_remat.hpp`.
**Motivates:** the ordered (DAG) evaluator crashes when a value carries more than
one batched mode of one space; the concrete failure is w20 aux+occ CCSD
(`TA_ASSERT is_range_set_congruent`). Simpler cases (one batched mode per space
per value) work.

This spec is **domain-neutral**: spaces are compared only as abstract colors,
never by concrete identity. The running example uses an abstract space `o`.

---

## 1. What breaks, and why (measured)

The ordered evaluator realizes a nest of batch loops and, per value, a map from
each enclosing loop to the physical mode it slices. When a value carries **two
batched modes of one space** (e.g. a rank-2 amplitude `T[i,j]` with both indices
batched), the current loop **identity cannot tell the two loops apart**, so both
modes are driven by one loop — the *collapse*. Measured:

- **w20 aux+occ** (ordered, two batched external modes of one space) dies with
  `TA_ASSERT is_range_set_congruent` (TiledArray `tensor/kernels.h:552`) on a
  multi-member product `C(i4,i3,μ) * I(i2,i1,μ) -> I(i2,i1,i4,i3)`: four modes of
  the batched space under **one** loop, each operand sliced on a single mode to
  the current block (`via=exact`) while the others stay full → incongruent
  result. (Traced from the recorded crash logs.)
- **w8** with the same class of terms runs **lossless** (batched energy Δ≈4.6e-13
  vs unbatched) because its ranges never form the incongruent block — the defect
  is latent at small sizes and **scales in**.
- The per-value role source (`classify_axis` → `per_axis`) was measured to emit
  **correct per-instance roles** on w8 (each batched external → LoopCarried, each
  contracted → Reduction, both members of a two-member group present as separate
  entries). So the roles are usable; the defect is downstream, in **loop
  identity and the mode↔loop map**, not in role classification.
- A *separate* w8 crash (`ordered_home_reads` premature eviction of a
  reduction-space intermediate) is **not** this bug and is out of scope here.

The root cause is that **loop identity is under-specified**. This spec defines it
completely and specifies how the map from loops to a value's modes is produced
and tracked end-to-end.

---

## 2. The pipeline

```
factorizer (per tree)      -> batched FOREST: each tree's nodes carry batched
                              modes in that tree's own index-frame; same-space
                              batched modes form an UNORDERED group.
compute_dag_boulevard      -> batched DAG: trees fused by structure; each value
   (fusion)                   has occurrences, each keeping its own frame.
rematerialize_to_budget    -> peak-driven placement (shrink/split); may break a
   (remat)                    group into subgroups/singletons. Peak only.
build_ordered_schedule     -> the realized loop nest + the per-occurrence
   (nest realization)         mode↔loop map. Assigns absolute loop numbering.
evaluate_ordered_multiroot -> execution; consumes the map via the slice seam.
   (runtime)
```

---

## 3. Loop identity: three orthogonal axes

A realized loop is identified by **`(depth, altitude_ordinal, latitude_ordinal)`**.
`space` is carried for diagnostics and abstract-color matching **only**; it is
never part of identity.

- **`depth`** — cross-space-group nesting position: which space-group encloses
  which. Fixed by the schedule's type order; one `depth` per space-group. (Two
  loops of *different* spaces are ordered by `depth`.)
- **`altitude_ordinal`** — which member *within* one same-space group of
  interchangeable loops. This is the vertical nesting position inside the group.
  The **order among members is free** (interchangeable); the ordinal fixes a
  chosen order. **This axis does not exist in the code today** — it is the
  missing piece.
- **`latitude_ordinal`** — which sequential **pass** from the legality
  producer→consumer (PROCON) split: one loop-nest re-traversed in sequence
  (producer pass, then consumer pass). This is the current code's `ordinal`
  (`emit_pass(0,…)`/`emit_pass(1,…)`, `ordered_schedule.hpp:1402-1405`), to be
  **renamed** `latitude_ordinal`.

`depth` and `altitude_ordinal` are both *spatial* (loops live on the stack
simultaneously); `latitude_ordinal` is *temporal* (one loop-nest, two passes in
sequence). Keeping them as three named fields removes the ambiguity that sank the
previous spec (which used one `ordinal` for two different jobs).

### Proposed types

```cpp
// dag_scope.hpp
struct DagScopeLevel {
  std::size_t depth;            // cross-space-group nest position
  int         altitude_ordinal; // member within a same-space group
  int         latitude_ordinal; // PROCON pass index (was: ordinal)
  std::wstring space;           // COLOR ONLY -- never identity
  // identity/hash/== over (depth, altitude_ordinal, latitude_ordinal) only.
};
using LoopId = std::size_t;     // index into the schedule's canonical level list
```

---

## 4. The per-occurrence mode↔loop atlas

For each **occurrence** (a use of a value in the DAG), the runtime needs a map
**`loop-id → position`**: given an enclosing loop, which physical slot of this
occurrence's result it slices. `position` is the slot index in the occurrence's
own `canon_indices()` frame (what `LoopColoredSliceSeam::by_hash` already stores,
`dag_scope.hpp:90-97`) — **frame-local, always available**.

The only quantity that must be *decided* is, per batched slot, **which group
member** (which altitude) it binds to. Two facts make this per-**occurrence**, not
per-value, and not a within-frame rank:

**(a) The map is not a rank.** Different occurrences of one loop assign it to
different positions. In the worked example below, `C[1,2]` maps position 0→member
A, but `E[2,1]` maps position 0→member **B**.

**(b) The correspondence is a tracked *choice*, because connectivity can
conflict.** Two occurrences that share a value can imply *contradictory*
member identifications, and fusion must pick one and record the transposition the
loser suffers.

### Worked example (why (a) and (b))

Forest (two same-space groups, `o`):
```
for {i1,i2}
  C[i1,i2] = A[i1,p] * A[p,i2]
  D[i1,i2] = A[i1,p] * A[i2,p]
for {i3,i4}
  E[i4,i3] = C[i3,i4] * D[i4,i3]
```
Fuse the two `{o,o}` groups. Trace the shared arrays:
- through **C** (`C[i1,i2]` produced, `C[i3,i4]` consumed): slot 0 ties `i1≡i3`,
  slot 1 ties `i2≡i4`;
- through **D** (`D[i1,i2]` produced, `D[i4,i3]` consumed): slot 0 ties `i1≡i4`,
  slot 1 ties `i2≡i3`.

**C and D disagree** (`i1≡i3` vs `i1≡i4`). There is no clean partition; fusion
**chooses** one — say follow C: member A = `{i1,i3}`, member B = `{i2,i4}` — and
the loser **D is read transposed**:
```
for {A@o, B@o}            // one nest, members A,B of space o
  C[A,B] = A[A,p] * A[p,B]
  D[A,B] = A[A,p] * A[B,p]
  E[B,A] = C[A,B] * D[B,A]   // E, D consumed in the chosen member frame
```
The other choice (member A = `{i1,i4}`) is equally valid and transposes C
instead. **Both compute the same result**; the choice is free (§10 non-goal:
optimal choice is future work).

Legality (PROCON) then finds `E` reads `C,D` cross-iteration within `{A,B}` and
**splits** the nest into two sequential passes — this is the `latitude` axis:
```
for {A@o, B@o}   // producer pass, latitude 0
  C[A,B] = ...
  D[A,B] = ...
for {A@o, B@o}   // consumer pass, latitude 1
  E[B,A] = C[A,B] * D[B,A]
```

The atlas the runtime consumes, per occurrence (position → member; absolute
altitude numbering applied at §5.4):
- `C` build: `0→A, 1→B`
- `D` build: `0→A, 1→B`; `D` **as consumed by E**: `0→B, 1→A` (transposed)
- `E` build: `0→B, 1→A`

### The tracked datum

Per occurrence, per batched slot, a **group-member id** (a small integer local to
the `(scope, space)` group), shared across occurrences that bind the same member.
This is the load-bearing data. Fusion assigns it by unifying occurrences over the
DAG's producer→consumer slot connectivity (a union-find), resolving conflicts by
a recorded choice. The **absolute `altitude_ordinal` numbering** (member-id →
0,1,…) is applied later (§5.4), so fusion commits to *no* order.

---

## 5. Who produces what, and when

| Stage | Produces | On the loop group |
|---|---|---|
| **Factorizer** (per tree) | per-node batched modes in the tree's own frame; same-space modes as an **unordered group** | order undecided |
| **Fusion** (`compute_dag_boulevard`) | per-occurrence `position → group-member-id`, unified over slot connectivity, conflicts resolved by a tracked choice | members identified; **order still free** |
| **Remat** (`rematerialize_to_budget`) | peak-driven shrink/split; **may break a group into subgroups/singletons**; peak only, prices no time | order still free; group *membership* may change |
| **Nest realization** (`build_ordered_schedule`) | the realized nest; assigns **absolute `altitude_ordinal`** + `depth`; realizes **`latitude`** PROCON passes → full `loop-id`; the per-occurrence `loop-id → position` map | order fixed here |
| **Runtime** | consumes `loop-id → position` per occurrence via the seam | — |

**Timing rule.** The *correspondence* (member identity, incl. transpositions)
is established at **fusion** and tracked through. The *absolute numbering*
(`altitude_ordinal`) is deferred to **nest realization** — because it is free
until then and remat may still change group membership. Remat does **not** need
the absolute numbering: it distinguishes values homed in different members by
their own home-mode positions (frame-local), which suffices for the peak sweep's
distinctness requirement.

### 5.4 Value-id coloring

A value is cache-/remat-keyed by **sliced-mode-id = space + loop-id** for each of
its home-sliced modes (`CachedValue{node, coloring}`, `value_id.hpp`). The
coloring must key by the **full `loop-id`** (`depth, altitude_ordinal,
latitude_ordinal`). Today it keys by **`depth` only**
(`value_id.hpp:42`, `colors.emplace(m, depth)`), and `home_mode_depth` stores
`pair<Index,int depth>` — dropping altitude *and* latitude. So two values homed
in different members of one group are colored identically today. Keying by the
full loop-id (available once §5's numbering is done, on the post-remat cells the
schedule is built from) makes storage (producer frame) and lookup (consumer loop)
agree on one identity.

---

## 6. What is broken today (code cites)

1. **Loop identity is incomplete.** `DagScopeLevel = {depth, space, ordinal}`
   (`dag_scope.hpp:34-47`); the level enumeration keys `(depth, space, ordinal)`
   (`ordered_schedule.hpp:1745`); `ordinal` is **latitude** only
   (`emit_pass 0/1`, `1402-1405`). **No `altitude`** → same-space group members
   collide on `(depth, space, latitude=0)`.
2. **The nest collapses the group.** `types` dedups by space
   (`ordered_schedule.hpp:992-999`); `depth_of_type` is space→depth
   (`1017-1023`); the escape loop folds a value's multiple same-space non-local
   modes into **one** escape per space-depth (`1049-1068`).
3. **Coloring is depth-only** (`value_id.hpp:42`; `home_mode_depth` drops the
   ordinal).
4. **Fusion tracks the correspondence by label set-arithmetic**
   (`compute_dag_boulevard` `home/enclosing/divergent_modes` folds,
   `peak_profile.hpp` ~490-513), which cannot represent the per-occurrence
   member map or the transposition of §4(b) robustly.
5. **The schedule is built from the seed, not post-remat.** The runtime schedule
   uses `rich = compute_dag_boulevard(roots)` (pre-remat) (`cck.ipp:2222-2244`);
   remat's result (`res.cells`) only builds a router overlay
   (`sequant.h:269-276`). So the schedule never sees remat's split/placement.
6. **The consumer-disambiguation seam is a workaround for missing altitude.**
   `LoopColoredSliceSeam::by_hash_consumer` / `mode_of`'s consumer arm
   (`dag_scope.hpp:102-186`) disambiguates "two modes under one loop" by
   consumer. With a real `altitude_ordinal`, each member is its own `LoopId` and
   each mode binds directly; this machinery can retire (or shrink to the
   genuinely-divergent case).

---

## 7. Design

1. **Complete the loop identity.** Add `altitude_ordinal`; rename `ordinal` →
   `latitude_ordinal` (§3). Identity/hash over the three fields.
2. **Track the per-occurrence member map through fusion** (§4). Fusion assigns
   per-occurrence `position → group-member-id` by unifying occurrences over the
   producer→consumer slot connectivity, resolving conflicts by a recorded choice
   (any valid choice for now, §10). This is a **new first-class datum**, not
   label set-arithmetic.
3. **Realize per-member loops in `build_ordered_schedule`.** Stop collapsing:
   emit one loop per group member at its space-depth, assign `altitude_ordinal`
   (a canonical numbering of the members), and realize `latitude` PROCON passes
   as today. Build from **post-remat cells** so the numbering and map reflect the
   real placement.
4. **Color value-ids by the full `loop-id`** (§5.4).
5. **The seam consumes per-occurrence `loop-id → position` maps.** With altitude
   present, `mode_of` resolves each member directly; retire the consumer-arm
   workaround where altitude now suffices.

Roles (`per_axis`, `classify_axis`) are **kept as-is** — measured correct
per-instance (§1). This design changes *loop identity and the mode↔loop map*, not
role classification.

---

## 8. Invariants

- **Frame purity.** No cross-frame **label** match anywhere. Cross-frame
  **space** comparison is allowed only as an abstract color (matching group
  spaces across trees during fusion); concrete space names never appear in
  scheduling/identity logic.
- **Per-occurrence.** The mode↔loop map is per occurrence; a value with divergent
  occurrences carries per-occurrence transpositions (§4b), never a single
  per-value rank.
- **Identity completeness.** Every realized loop has a distinct
  `(depth, altitude_ordinal, latitude_ordinal)`; no two same-space group members
  share a loop-id.
- **Numbering freedom.** `altitude_ordinal` is a free canonical choice; the
  schedule must be correct for *any* consistent numbering.

---

## 9. Validation

- **w8 forced-multi-member-batching**: crash-free, zero unresolved facts, batched
  energy matches unbatched to full precision.
- **w20 aux+occ** (the `is_range_set_congruent` crash): completes and matches the
  unbatched/aux-only reference.
- Existing `[ordered-schedule]` unit tests stay green.
- A schedule-level dump shows, for a two-member-group value, **two distinct
  loop-ids** with distinct `altitude_ordinal`, and the per-occurrence map
  reproduces the §4 transposition.

---

## 10. Non-goals

- **Optimal fusion choices.** Where connectivity conflicts (§4b), fusion picks
  *a* valid member assignment; choosing the *best* one (fewer transpositions /
  less split recompute) is a **future** task. Correctness of the basic machinery
  (batched forest → DAG → schedule → execution) comes first.
- **Optimal altitude order / group reordering across spaces.** The canonical
  numbering is fixed; a future optimization may reorder.
- No change to the factorizer, the role classifier, or the reduction/contracted
  batch path beyond giving contracted modes the same per-member identity.
- The separate `ordered_home_reads` premature-eviction bug is tracked elsewhere.

---

## 11. Why the previous spec was withdrawn

`2026-08-27-frame-correct-use-induced-slicing-design.md` reasoned about `ordinal`
as the within-group selector (what is actually **altitude**), while the code's
`ordinal` is **latitude** (PROCON passes). It therefore proposed "distinct depths
per instance" and mixed the two axes. Its **measured** content (the w20 trace, the
w8 losslessness, the per-instance role verification) is preserved in §1 here; its
model sections are replaced by §3-§7.
