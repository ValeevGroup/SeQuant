# Perfect-CSE `home_scope` seed, O2-owned demotion, and why the meet (not the occurrence-key) computes the seed home

Status: design, agreed via brainstorming 2026-08-03. Supersedes two earlier same-day
drafts (one proposing to subsume `demoted_external` into the occurrence-key -- retracted;
one proposing to *retain* the `has_demoted_external` veto in a byte-identical seed --
also corrected here). Amends `doc/dev/specs/2026-08-02-home-scope-placement-design.md`
(Sections 7b, 7d, and the roadmap). Depends on Phase 1 (peak oracle) and Phase 2 (router
as a read+store override seam), both complete on
`evaleev/feature/multimode-batched-eval` (`...0bae1273e`).

## 1. The model, and the correction

Placement is **perfect-CSE seed + O2 spill pass** (register allocation): the seed hoists
every value to its maximal shared home (max sharing => peak-MAXIMAL), and O2 relaxes that
just enough to fit the peak budget by evicting/splitting the cells alive at the binding
peak point.

The correction this note records: today's `place_at_this_level` bakes an **ad-hoc spill
decision into the seed** -- `has_demoted_external` (eval.hpp:1678-1706) un-hoists any value
whose meet-dropped external mode would otherwise materialize a full "scattered giant." That
is precisely an O2 decision (an over-budget cell at a binding peak point), pre-empted by a
hard-coded "don't hoist." Since we are building the O2 pass to make exactly this
peak-vs-recompute trade cost-based, the veto is redundant with it and strictly worse (it
always un-hoists to local; O2 would weigh un-hoist vs split vs shrink).

There is **no external/contracted asymmetry** to preserve. A contracted batched mode
hoisted above its loop is a big full copy too (`A[c,_]` and `A[_,c]` -- two occurrences
carrying a contracted `c` in different slots -- perfect-CSE hoists `A` above the `c`-loop,
full on `c`; O2 weighs that copy exactly like an external giant). "External vs contracted"
is a property at the loop node, irrelevant to how a hoisted copy below it is sized/spilled.
So `has_demoted_external` is a hack, and it is **removed**; the seed has zero
external/contracted special-casing.

## 2. The seed: pure perfect-CSE `home_scope` from the meet

`home_scope(value)` = deepest enclosing scope over its residency

```
residency = sliced_modes (cross-occurrence meet, stamp_lifetime_masks)  UNION  contracted_modes
home_scope = deepest enclosing batch loop whose mode is in residency, else chain root
```

with **NO veto and NO demotion fold**. For a demoted value (an external mode dropped by the
meet) this homes the value ABOVE that loop -- full on that mode, i.e. the giant -- BY
DESIGN. The seed is peak-maximal; O2 (Section 4) makes it feasible.

The **meet is load-bearing** here and cannot be replaced by the occurrence-key (Section 3):
it is the cross-occurrence set-INTERSECTION that produces the shared home; the key can
separate but not intersect.

`home_scope` is a STATIC predictor (computed without running eval) consumed by the peak
profile (Phase 3b) and O2's `ΔPeak` (Phase 4). Spec 7d's demotion-FOLD ("home the carrier
inside its external loop, above invariant inner loops") is NOT a seed rule -- it is one
candidate O2 move (Section 4).

## 3. Negative result: the occurrence-key cannot replace the meet

We explored deriving residency from the occurrence-key-class instead of the meet (to also
drop `stamp_lifetime_masks`). Refuted. The occurrence-key partitions occurrences by their
EXACT batched-slot pattern -- it can SEPARATE but has no operation that INTERSECTS.
Partial-overlap occurrences of one value -- `A[i1,i2]` (slices `{i1,i2}`) and `A[i1,_]`
(slices `{i1}`) -- must share ONE cell at their intersection home `A[i1,_]`, each slicing
further on use; that home is a set-intersection only the meet computes. Probe (faithful:
one canonical node, single meet bucket): the key gives `{A[i1,i2], A[i1,_]}` DIFFERENT
keys under both the shipping and a mixed-tier coloring -- the split is structural
(named-set cardinality 2 vs 1), not a coloring artifact -- while the meet gives `{i1}`.
Deriving residency from the key-class would drop that shared cell, regressing legitimate
CSE into avoidable recompute. So the meet stays.

Corollary: the mixed-tier occurrence-key coloring (which does structurally separate the
demoted cross-pair occurrences) is NOT needed for the seed. It is at most a Phase-4
question -- whether O2 ever needs the router to override such occurrences individually --
and is deferred.

## 4. O2 owns demotion (and all spill)

A demoted giant is just an over-budget cell alive at a binding peak point. O2's moves --
shrink (re-home into an existing carried loop), un-hoist (lower the home toward a loop it
is invariant to), split (partition use-sites into shorter-lived cells) -- handle it
cost-based. The current veto ("home local") and 7d's fold ("home inside the external
loop") are both candidate O2 moves, chosen by `ΔPeak/ΔRecompute`, not hard-coded.

## 5. Consequences

- **Not byte-identical.** The seed differs from today's heuristic on demoted values (the
  seed hoists them; the heuristic vetoes to local). So "validate `home_scope` == runtime
  store scope" is wrong for demoted values. Validate the seed against its DEFINITION (the
  meet home); validate the whole placement (seed + O2) against the peak oracle + result
  correctness.
- **The seed alone is peak-maximal -- worse than today until O2 lands.** So the seed
  (Phase 3) and O2 (Phase 4) are a COUPLED unit: keep today's heuristic placement as the
  default/fallback until O2 is ready (gate the perfect-CSE-seed + O2 path behind a flag),
  rather than shipping a regressed seed standalone.
- **Removes `has_demoted_external`** (and the meet-veto coupling in `place_at_this_level`).

## 6. Roadmap (revised)

- **Phase 3a -- pure perfect-CSE `home_scope` seed** from the meet (`sliced ∪ contracted`,
  no veto, no fold). Static predictor; validated against its definition (the meet home).
- **Phase 3b -- static peak profile** (spec 7c / O3a): weighted-interval sweep over the
  seed + schedule, validated to equal the Phase-1 replay oracle (with co-resident summing).
- **Phase 4 -- O2 greedy** (spec 7b): owns demotion/spill (shrink/un-hoist/split by
  `ΔPeak/ΔRecompute`); lands COUPLED with the seed behind a flag, replacing the current
  heuristic (incl. the removed veto) when enabled.
- **Phase 5 -- feedback** (spec 7b / O6).

## 7. Open items

- **O-a: static meet-home computation** validated against its definition (the deepest
  `sliced ∪ contracted` enclosing loop) across the witness forests.
- **O-b: the seed+O2 coupling / flag** so the perfect-CSE seed never ships as a regressed
  standalone; today's heuristic stays the default until O2 is feasible.
- Retired by the analysis above: residency-from-key (refuted); the `index_color`
  occurrence-key change (not needed for the seed); `has_demoted_external` (an O2 move, not
  a seed rule).
