# Subsuming `demoted_external` into a corrected occurrence-key

Status: design, agreed via brainstorming 2026-08-03. Amends
`doc/dev/specs/2026-08-02-home-scope-placement-design.md` (Sections 7d, 9, and the
roadmap). Depends on Phase 1 (peak oracle) and Phase 2 (router as read+store override
seam), both complete on `evaleev/feature/multimode-batched-eval` (`...0bae1273e`).

## 1. Problem and finding

Phase 2 landed the occurrence-key as `canonicalize_slots(node sub-network,
named_indices = the in-scope batched modes on the node's slots)` ->
`SlotCanonicalizationMetadata`. The open Phase-3 question was whether the
`has_demoted_external` compensation in `place_at_this_level` (eval.hpp:1678-1701) is a
genuine model component or a removable hack.

Two probes settled it:

- **Demotion is a RUNTIME-block incompatibility, not a structural one.** A demoted
  external carrier is the SAME canonical value in occurrences that bind the external occ
  slot to different, proto-incompatible PNO pairs (`a<i1,i2>` vs `a<i3,i4>` -- the
  cross-pair two-PNO giants). The cross-occurrence meet (`stamp_lifetime_masks`)
  intersects the occ modes by Index identity (`{i1,i2} INTERSECT {i3,i4} = {}`), empties
  `sliced_modes`, and fires the guard. The occurrences differ only by concrete labels of
  one space -- runtime, not structure.

- **The occurrence-key as SHIPPED collapses them, but a corrected key SEPARATES them
  losslessly.** `TensorNetworkV3::canonicalize_slots` hardcodes
  `distinct_named_indices=false` (v3.cpp:697), so a named index is colored by
  `Index::color()` = space only; `i1..i4` share space `i` -> identical bliss graphs ->
  collapse. Coloring the BATCHED externals by their position in a shared FOREST-GLOBAL
  named set instead separates them, and -- verified -- does so WITHOUT breaking the
  symmetric-domain collapse the key must preserve (see Section 3).

**Conclusion.** `demoted_external` is a symptom of two fixable defects in the current
occurrence-key and residency derivation, not a genuine model primitive:
1. the key colors batched externals by space, discarding the block-identity that
   distinguishes cross-pair occurrences (Section 3);
2. residency comes from the cross-occurrence meet, which is exactly what empties on the
   block-incompatible collision (Section 4).

Fix both and demoted carriers become distinct key-classes, each homed inside its own
block loop, each with per-block residency -- so the scattered giant is never formed and
the `has_demoted_external` veto is retired (Section 5).

## 2. Why demotion is genuinely a runtime concept the key can and must encode

Two occurrences of one canonical value that bind an external batched mode to different
blocks are the same STRUCTURE but different DATA (different PNO domains). For placement
they must not share a hoisted full copy. The distinguishing signal is the concrete
batched index each carries (the block), which is:
- not the label per se (labels are arbitrary), but
- the index's IDENTITY within the forest -- which, given a consistent global numbering,
  distinguishes `i1` from `i3` while treating same-space non-batched externals uniformly.

Complete invariance to external relabeling is neither achievable nor required (you cannot
bootstrap a canonical external ordering -- topo-sorting external slots needs a coloring,
which needs the ordering; and bliss canonical form is heuristic and choice-dependent
anyway). The bar that IS met and sufficient: **forest-consistent determinism** -- the
same occurrence keyed at the store site and the read site produces the same graph. Bliss
gives that (same colored graph in -> same canonical form out), and store/read key the
same node with the same labels within one forest evaluation.

## 3. Linchpin 1 -- the mixed-tier occurrence-key coloring

**Requirement (three tiers):**
- BATCHED externals: a distinct color per index, sourced from a shared FOREST-GLOBAL
  batched-index numbering (so `i1` and `i3` differ). Separates demotion; pins each
  batched slot by color uniqueness.
- NON-BATCHED externals: kept NAMED (so their external slots stay pinned -- this is what
  preserves `A[i1,_]` vs `A[_,i1]` distinctness and the slot structure generally) but
  colored by SPACE (so they do not inject label dependence into the canonical slot
  ordering).
- INTERNALS: anonymous / space, as today.

The single `distinct_named_indices` bool cannot express this -- it colors ALL named
indices uniformly. A per-index color override is required.

**Mechanism (validated against the exact router path).** `VertexPainterImpl::operator()
(const Index&)` (vertex_painter.cpp) currently branches named/anonymous and, for named,
`distinct` vs `idx.color()+0xabcd`. Add an optional, non-owning per-index color hook,
consulted first; a returned value overrides, `nullopt` falls through to the existing
default (preserving the `+0xabcd` named-vs-internal shift and `ensure_uniqueness`):

```cpp
// create_graph options + canonicalize_slots parameter (nullable, BY POINTER so the
// stateful closure is never copied per-call -- mirrors `const NamedIndexSet* named_indices`):
const std::function<std::optional<Color>(const Index&)>* index_color = nullptr;
```
`std::optional` makes "no value => use default" manifest; the pointer means the caller
owns the closure (built once per forest, capturing the global batched set by reference)
and `canonicalize_slots`/`create_graph`/`VertexPainter` only dereference it.

**The occurrence-key call becomes:** `named_indices = all externals` (default -- every
external slot pinned), `distinct_named_indices = false` (unchanged), plus
`index_color = [&global_batched](const Index& i) -> std::optional<Color> {
  return in(i, global_batched) ? std::optional{ordinal_of(i, global_batched)}
                               : std::nullopt; }`.
The forest-global batched set (a label-sorted `NamedIndexSet` of every batched external
in the forest) is built once and threaded to every per-node key call, so the ordinal is
GLOBAL, not per-node (per-node positions collapse demotion; global positions separate
it).

**Why it does not over-separate (validated, "needle threaded").** Distinct coloring pins
index IDENTITY; the tensor's declared symmetry pins SLOT interchangeability -- orthogonal.
A symmetric PNO domain in opposite order (`a<i1,i2>` vs `a<i2,i1>`) is a slot permutation
over the SAME two colors (`{0,1}` present in both graphs) -> the automorphism fires ->
COLLAPSE (correct). A different pair (`a<i1,i2>` vs `a<i3,i4>`) is a different color SET
(`{0,1}` vs `{2,3}`) -> no permutation reconciles -> SEPARATE (correct). Measured through
`create_graph -> bliss canonical_form -> ConstGraphCmp::cmp` on all of: demotion
(SEPARATE), symmetric-order (COLLAPSE), gC-style symmetry with both occ batched
(COLLAPSE), gC different-pair (SEPARATE), and identical/diff-label controls -- all as
desired.

## 4. Linchpin 2 -- residency from the key-class, not the meet

Today residency (`sliced_modes`) = the cross-occurrence meet (`stamp_lifetime_masks`),
which empties on the block-incompatible collision -- the root of the demotion symptom.
Once the corrected key SEPARATES demoted occurrences into distinct key-classes, the meet
is unnecessary: within a key-class every occurrence has the same batched-slot structure
BY CONSTRUCTION, so **residency = the batched (`sliced UNION contracted`) modes that
key-class carries** -- a per-key-class quantity, no intersection. `home_scope` = deepest
enclosing loop over that residency.

- NON-DEMOTED values: their occurrences already agree; same key-class; residency = what
  the meet gave -> placement UNCHANGED.
- DEMOTED values: now distinct key-classes; each carries its own block's residency ->
  homed INSIDE its own block loop -> no full-extent scattered giant; `has_demoted_external`
  falls out.

This retires the `has_demoted_external` veto in `place_at_this_level` and replaces the
meet-based `sliced_modes` derivation. (Design-level; the exact per-key-class residency
computation is an implementation task -- Section 6 -- and MUST be validated, not assumed:
non-demoted homes unchanged, demoted homes per-block.)

## 5. Consequences

- **Not byte-identical.** This CORRECTS demoted-value placement (per-block instead of
  vetoed-to-local) while leaving non-demoted placement identical. The two hidden witnesses
  re-baseline; the occ-veto aux+occ leg should IMPROVE (the cross-pair giant it carries is
  exactly what per-block homing removes). Non-demoted figures must be unchanged.
- **Validation shifts** from "byte-identical" to: (a) non-demoted placement/peak unchanged;
  (b) demoted values homed per-block (no full-extent giant), confirmed against the eval
  trace; (c) the static/replay peak oracle agrees; (d) an integrated store->read relocation
  test where the relocated value CARRIES a batched mode (the non-empty-named-set round trip
  the Phase-2 review flagged as a prerequisite -- Section 6).
- **It amends Phase 2's occurrence-key** (the coloring was incomplete) and **moots the
  "demotion veto in `home_scope`"** design that the byte-identical framing implied.

## 6. Reshaped roadmap

The parent spec's roadmap (P3 = static peak sweep) is revised to insert the
demotion-subsumption work first, since it corrects the seed placement the sweep profiles.

- **Phase 3a -- corrected occurrence-key coloring.** Add the `index_color` hook to
  `create_graph`/`VertexPainter` and expose it on `canonicalize_slots` (nullable, by
  pointer, `nullopt` => default -- all existing callers byte-unchanged). Build the
  forest-global batched-index numbering and thread it. Occurrence-key uses
  `named_indices = all externals` + the `index_color` closure. Tests: the Section-3 table
  (demotion separates, symmetry collapses, controls), plus read-site == store-site key
  equality for a batched-mode-carrying value.
- **Phase 3b -- residency from key-class + retire `demoted_external`.** Replace the
  meet-based `sliced_modes` with per-key-class residency; remove the `has_demoted_external`
  veto in `place_at_this_level`. Validate: non-demoted homes unchanged; demoted homes
  per-block; witnesses re-baseline (non-demoted figures identical, demoted leg improved);
  the integrated batched-mode relocation round-trip.
- **Phase 3c -- static peak profile (the former Phase 3, spec 7c/O3a)** on the corrected
  seed, validated == the Phase-1 replay oracle.
- **Phase 4 -- O2 greedy** (spec 7b); **Phase 5 -- feedback** (spec 7b/O6). §7d's
  demotion-FOLD (home a carrier inside its external loop but above invariant inner loops)
  is a Phase-4 O2 quality move, not needed here.

## 7. Open items to validate in implementation (do not assume)

- **O-a: forest-global batched-index numbering.** Confirm the label-sorted global
  `NamedIndexSet` is computed identically at every store and read site of one forest (the
  determinism the key rests on), and cheap enough to build once per forest.
- **O-b: residency-from-key correctness.** Linchpin 2 is reasoned but NOT yet probed.
  Verify per-key-class residency reproduces the current non-demoted homes exactly and
  yields per-block demoted homes -- this is the gate for retiring the meet and the veto.
- **O-c: the `index_color` hook must not perturb other `canonicalize_slots` callers.**
  Default `nullptr` keeps `build_subnet_metadata` etc. byte-identical; confirm.
