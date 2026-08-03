# Meet-based `home_scope` for Phase 3, and why the occurrence-key cannot subsume demotion

Status: design, agreed via brainstorming 2026-08-03. Supersedes an earlier
same-day draft that proposed subsuming `demoted_external` into the occurrence-key
(retracted -- see Section 1). Amends `doc/dev/specs/2026-08-02-home-scope-placement-design.md`
(Section 7d and the roadmap). Depends on Phase 1 (peak oracle) and Phase 2 (router as a
read+store override seam), both complete on `evaleev/feature/multimode-batched-eval`
(`...0bae1273e`).

## 1. Negative result: the occurrence-key cannot replace the meet

We explored eliminating the `has_demoted_external` compensation in `place_at_this_level`
(eval.hpp:1678-1701) by (i) a mixed-tier occurrence-key coloring that separates demoted
occurrences, and (ii) deriving residency from the occurrence-key-class instead of the
cross-occurrence meet (`stamp_lifetime_masks`). Result: **the idea is refuted; the meet is
load-bearing and `has_demoted_external` stays.**

What the probes established (all reverted, nothing committed):

- **Demotion is a RUNTIME-block incompatibility.** A demoted external carrier is the same
  canonical value in occurrences binding the external occ slot to different proto-
  incompatible PNO pairs (`a<i1,i2>` vs `a<i3,i4>`, the cross-pair giants); the meet
  intersects the occ modes by Index identity to empty, firing the guard. The occurrences
  differ only by concrete labels of one space.
- **A mixed-tier coloring CAN separate demoted occurrences** (batched externals colored by
  a shared forest-global ordinal via a `distinct`/`index_color` hook), and does so without
  breaking the symmetric-domain collapse. Recorded for reference -- it may matter for
  Phase-4 router-override granularity -- but it is NOT used in Phase 3.
- **Residency-from-key is refuted.** The occurrence-key partitions occurrences by their
  EXACT batched-slot pattern: it can SEPARATE but has no operation that INTERSECTS.
  Partial-overlap occurrences of one value -- `A[i1,i2]` (slices `{i1,i2}`) and `A[i1,_]`
  (slices `{i1}`) -- must share ONE cell at their intersection home `A[i1,_]`, each slicing
  further on use; that home is a cross-occurrence set-INTERSECTION only the meet computes.
  Probe (faithful: one canonical node, `TreeNodeEqualityComparator` YES so one meet
  bucket): the key gives `{P=A[i1,i2], Q=A[i1,_]}` DIFFERENT keys under BOTH the shipping
  and the mixed-tier coloring -- the split is structural (named-set cardinality 2 vs 1,
  non-isomorphic graphs), not a coloring artifact -- while the meet gives them `{i1}`
  (non-empty). Deriving residency from the key-class would drop that shared cell,
  regressing legitimate CSE into avoidable recompute.

**Why the meet is exactly right.** Demotion is the meet-EMPTY case (`A[i1,i2]` vs
`A[i3,_]`, meet `{}`). The meet handles BOTH regimes with one operation: partial overlap
-> non-empty intersection -> a real shared home; disjoint -> empty -> full home ->
`has_demoted_external` declines to hoist the scattered giant. The occurrence-key can do
neither, so it cannot subsume the meet or the demotion veto.

**Corollary.** The occurrence-key coloring fix is NOT needed for Phase 3 -- `home_scope`
comes from the meet, not the key. It is at most a Phase-4 question (whether the O2 pass
ever needs the router to override demoted occurrences individually). Deferred.

## 2. Phase 3 = meet-based `home_scope` (the direction)

`home_scope(occurrence)` is a STATIC replica of `place_at_this_level`'s runtime decision,
computed without running eval:
- **is it hoisted?** `residency_all_outer(n) && !has_demoted_external(n) &&
  batch_order_aware(n) && !subtree_any(n, is_volatile)` (eval.hpp:1671-1706);
- **at what level?** `rl + 1`, where `rl` is the deepest enclosing `batch_context` index
  whose mode is in the residency union `sliced_modes() UNION contracted_modes()` (from the
  meet; `in_union`/`rl` walk, eval.hpp:1724-1739); `rl == -1` => chain root;
- **else** (not hoisted, incl. `has_demoted_external`): local / at-use.

Section-7d reconciliation (settled): this is the LIVE rule -- residency `sliced UNION
contracted`, `has_demoted_external -> local`. Spec 7d's demotion-FOLD (`home = deepest over
sliced UNION demoted_external`, i.e. home the carrier INSIDE its external loop but above
invariant inner loops) is a strictly BETTER placement and a **Phase-4 O2** quality move
(reachable as a relocation override), NOT the Phase-3 seed.

`home_scope` is a static PREDICTOR consumed by the peak sweep (Phase 3b) and Phase-4
`ΔPeak`; it does NOT populate the runtime router (the seed stays derived / empty-router,
per the Phase-2 override-seam model). So Phase 3a adds a function plus its validation and
changes no runtime placement -- byte-identical.

**Validation (per-value -- the reason we split B, not A).** The static
`home_scope(occurrence)` must resolve to the SAME scope the runtime actually stores/finds
it at: emit `place_at_this_level`'s chosen store scope in the eval trace and assert
equality per hoisted value (the T3 shadow-assert idea generalized). Per-cell home equality
is tighter than an aggregate peak match, and O2 (Phase 4) consumes per-cell homes. Plus the
integrated store->read relocation test where the relocated value CARRIES a batched mode --
the non-empty-named-set round trip the Phase-2 whole-branch review flagged as a Phase-3
prerequisite.

## 3. Roadmap (revised)

- **Phase 3a -- static meet-based `home_scope` predictor** + per-value validation against
  the runtime store scope + the integrated batched-mode relocation test.
- **Phase 3b -- static peak profile** (spec 7c / O3a): the weighted-interval sweep over the
  `home_scope` predictor + schedule, validated to equal the Phase-1 replay oracle.
- **Phase 4 -- O2 greedy** (spec 7b), including 7d's demotion-fold as an O2 relocation move.
- **Phase 5 -- feedback** (spec 7b / O6).

## 4. Open items to validate in implementation

- **O-a: the per-value home oracle.** Emit `place_at_this_level`'s chosen store scope per
  hoisted value so the static `home_scope` predictor can be asserted equal to it. Confirm
  the static computation reproduces every non-demoted home exactly and the demoted/local
  homes too.
- Retired by the negative result: the residency-from-key derivation (refuted) and the
  `index_color` occurrence-key change (not needed for Phase 3; possible Phase-4 item).
