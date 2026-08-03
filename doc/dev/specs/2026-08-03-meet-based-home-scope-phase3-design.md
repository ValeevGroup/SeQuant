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
residency  = sliced_modes
sliced_modes = cross-occurrence MEET (stamp_lifetime_masks) of ALL batched modes
               (any BatchModeType) that live on the node's own result slots
home_scope = deepest enclosing batch loop whose mode is in residency, else chain root
```

with **NO veto and NO demotion fold**. For a demoted value (a mode dropped by the meet)
this homes the value ABOVE that loop -- full on that mode, i.e. the giant -- BY DESIGN. The
seed is peak-maximal; O2 (Section 4) makes it feasible.

**Unify `sliced_modes`; remove `contracted_modes`.** Today residency is `sliced_modes UNION
contracted_modes`, an artifact: `sliced_modes` is built External-only (`stamp_lifetime_masks`
walks `ext_modes_of`, which filters `batched_here()` to `BatchModeType::External`), so it
"structurally cannot express a node variant to an aux loop it carries free," and
`contracted_modes` was bolted on to patch that gap -- emitted *per-occurrence* by the cost
model (`cost_model.hpp:1898-1911`), never meet'd. There is no principled reason a node inside
an aux loop differs from one inside an occ loop; both make it variant to that loop's index.
The fix: drop the `External` filter so `stamp_lifetime_masks` meets ALL batched modes on the
result slots, and delete `contracted_modes` (field, accessors, cost-model emission,
`NodeBatchAnnotation::contracted_modes`, and the `in_union` union). This also *corrects* a
latent bug -- the per-occurrence `contracted_modes` under-meets aux modes that differ across
occurrences (`A[c,_]` vs `A[_,c]`); routing them through the same meet fixes it. A node that
*sums* `m` (m not on its result slots) still homes outside the `m`-loop correctly: the walk
passes `m` DOWN to the operands (which carry it and get it in their `sliced_modes`) while the
node's own contribution is `acc INTERSECT its slots`, excluding `m`.

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

- **Phase 3a -- the `home_scope` seed PREDICTOR (non-regressing; scope decision (A)).** Build
  the unified-meet `home_scope(value)` as a NEW static computation -- the cross-occurrence meet
  of ALL batched modes on the node's result slots (Section 2), then `home_scope = deepest loop
  over that residency`, no veto, no fold -- consumed by 3b and O2. It does NOT touch runtime
  `place_at_this_level` / `stamp_lifetime_masks`, so runtime placement stays today's heuristic
  (NO regression) and `contracted_modes` / the veto stay in place for now. Deliverables: the
  predictor, its validation against its definition (the deepest all-batched-modes-meet loop),
  and reconciling the CAT-2 `occurrence_key`/`ext_modes_of` inconsistency. The CAT-1 deletions +
  the `stamp_lifetime_masks` unification are the **CUTOVER**, deferred to Phase 4.
- **Phase 3b -- static peak profile** (spec 7c / O3a): weighted-interval sweep over the
  seed predictor + schedule, validated to equal the Phase-1 replay oracle (co-resident summing).
- **Phase 4 -- O2 greedy + CUTOVER** (spec 7b): owns demotion/spill (shrink/un-hoist/split by
  `ΔPeak/ΔRecompute`); lands COUPLED behind a flag. The cutover to the perfect-CSE seed --
  unify `stamp_lifetime_masks` (drop the `External` filter), DELETE `contracted_modes` +
  `has_demoted_external` (CAT-1) -- happens HERE, replacing today's heuristic when the flag is on.
- **Phase 5 -- feedback** (spec 7b / O6).

## 7. Open items

- **O-a: static meet-home computation** validated against its definition (the deepest
  `sliced_modes` enclosing loop, with `sliced_modes` = the unified all-batched-modes meet)
  across the witness forests.
- **O-b: the seed+O2 coupling / flag** so the perfect-CSE seed never ships as a regressed
  standalone; today's heuristic stays the default until O2 is feasible.
- Retired by the analysis above: residency-from-key (refuted); the `index_color`
  occurrence-key change (not needed for the seed); `has_demoted_external` (an O2 move, not
  a seed rule).

## 8. Audit: placement/residency hacks and where they go (2026-08-04)

A read-only sweep of the placement/residency/batching code (`eval.hpp` `place_at_this_level`,
`lifetime_mask.hpp`, `eval_expr.{hpp,cpp}`, `cost_model.hpp`, `node_batch_annotation.hpp`,
`cache_manager.hpp`, `schedule_dump.hpp`) inventoried the ad-hoc heuristics the clean model
retires. **Key structural finding: the entire `contracted_modes` mechanism exists solely
because `ext_modes_of` (`lifetime_mask.hpp:73-83`) is hardwired to `BatchModeType::External`.
Dropping that one filter clause collapses `sliced ∪ contracted` in `place_at_this_level` back
to `sliced` alone, making the CAT-1b removal a wholesale DELETE, not a rewrite.**

### CAT-1 -- confirmed DELETE (Phase 3a)
- `has_demoted_external` veto: `eval.hpp:1720-1743` (lambda+comment), `:1748` (collect-gate
  conjunct), `:1784` (comment); tests `test_eval_ta.cpp:3019,3033,3050,3053,3118,3124,3175`
  (expectations invert under the hoisting seed).
- `contracted_modes` end-to-end: `eval_expr.hpp:314-334` (get/set/doc), `:398-399` (member);
  `node_batch_annotation.hpp:25,33,44-46`; `cost_model.hpp:1885-1911` (doc+emission);
  `eval_expr.cpp:604` (threading); `eval.hpp:1663-1666,1716-1718,1766-1770` (`residency_all_outer`
  contracted loop + `in_union` union half + doc); `schedule_dump.hpp:133-134`; tests
  `test_eval_ta.cpp:2743-2760,2837-2845,2883,4772-4773`, `test_eval_dryrun.cpp:4320`.

### CAT-2 -- External-only assumptions to GENERALIZE (Phase 3a)
- `lifetime_mask.hpp:73-83` -- `ext_modes_of`'s `BatchModeType::External` filter. THE
  load-bearing line; drop it so `sliced_modes` = meet of ALL batched modes. HIGH.
- `lifetime_mask.hpp:31-52,86-101` -- surrounding docs assume External (`slot_modes_of` code
  already generalizes). HIGH.
- `occurrence_key.hpp:17-67` -- `in_scope_batched_on_node` is already kind-agnostic but claims
  to "mirror `ext_modes_of` exactly"; today they DISAGREE -- a latent Phase-2-key
  inconsistency to reconcile. MED.
- `cost_model.hpp:1841-1882` -- External-before-Contracted emission order has a real
  scatter-widening effect; flag, do not blind-delete. LOW-that-it's-pure-hack.

### CAT-3 -- genuinely O2-adjacent (code-verified 2026-08-04; the audit's first pass over-reached)
Re-reading the code shrank this list to two -- most CAT-3 candidates turned out to be seed cache
semantics (moved to CAT-5) or compile-time factorizer knobs (CAT-3b):
- `cost_model.hpp:1912-1920`, `eval_expr.hpp:344-353` -- `effective_count`/`batch_effective_count`,
  the integer reuse count = O2's rational `W` recompute measure. Genuinely -> `W`.
- `cache_manager.hpp:633-644` -- `max_footprint` "refuse to cache anything bigger than X."
  Default `0` => INERT (perfect CSE). Only its NON-default use is a seed-baked peak trade (the
  crudest spill, at cache-construction) that O2's evict-where-peak-binding supersedes.

### CAT-3b -- compile-time factorizer/cost-model knobs O2 RETIRES by taking over placement (not "moved")
- `cost_model.hpp:688,1868-1876` + `optimize.cpp:110` -- `node_level_placement`: the DP deciding
  per-node placement INSIDE `optimize()`. O2 (a post-factorizer placement pass) subsumes that
  FUNCTION, so the knob becomes moot -- not code relocated into O2.
- `cost_model.hpp:636,848-1040,1892-1897` + `optimize.cpp:109` -- `order_aware_recompute`: the DP's
  ordered-cell recompute-accounting mode; its placement-recompute part is superseded by `W`.
  Compile-time; not consulted at runtime.

### CAT-4 -- other smells
- `eval.hpp:1747` (+`eval_expr.hpp:355-368`, `node_batch_annotation.hpp:49-56`,
  `schedule_dump.hpp:118`) -- `batch_order_aware` gate; its own comment admits it's an OFF-path
  discriminator workaround.
- `cache_manager.hpp:641-649,713-744` -- the batch-variant veto (negative-half placement); to be
  subsumed by the router, not deleted in isolation. MED.
- `eval.hpp:1782-1795` -- router override replaces the LEVEL only; today two placement
  authorities coexist (router seam + `place_at_this_level` heuristics). MED.
- `eval.hpp:1536,1598` -- `depth < 8` magic recursion backstop, unexplained.
- `eval.hpp:1639-1652` -- "absent from `batched_here()` == Contracted, byte-identical" rationale;
  evaporates under uniform-mode residency.

### CAT-5 -- genuine, KEEP (do not over-delete)
- `cache_manager.hpp:129-164` (+ `is_volatile` gate `eval.hpp:1616-1626`) -- `is_volatile` is the
  P/NP cache-PERSISTENCE classifier (a volatile amplitude intermediate is rebuilt every CC
  iteration => must be NP; persisting it across iterations is a stale-data CORRECTNESS bug). Not a
  peak trade. The `persistent_only` batching gate is opt-in, default OFF, about not-wasting-batching
  effort. KEEP.
- `cache_manager.hpp:555-589` -- `min_repeats` is the CSE-cacheability threshold at cache
  construction; default `2` = perfect CSE ("cache everything reused"). Part of the seed's
  *what-is-CSE-able*, NOT a spill knob and NOT replaced by `W`. KEEP (only a user `>2` throttle is a
  peak knob O2 supersedes).
- `eval.hpp:1857+` -- External SCATTER (`write_into_slice`, disjoint) vs Contracted ACCUMULATE
  (`add_inplace`): the real runtime meaning of `BatchModeType` (`fwd.hpp:21`) -- irrelevant only
  to *residency*, not removed.
- `eval.hpp:1713-1719` -- `residency_all_outer` as the "home = deepest residency loop" walk
  (structural, not a peak trade); only its `contracted_modes`/veto conjuncts go, input
  generalizes via CAT-2.
- `eval_expr.hpp:340-342` + `lifetime_mask.hpp` -- `sliced_modes`/`mask_all_full`/
  `stamp_lifetime_masks`: the meet the whole model is built on.
- `single_term_detail.hpp:296-316`, `cost_model.hpp:643-667` -- `is_batchable_{contracted,
  external}_index` (batching candidacy; spec 9.4 untouched).
- `cost_model.hpp:581-594` -- `accumulation_factor` (arithmetic cost), `peak_threshold` (O2's
  constraint).

**Scope for the plan:** Phase 3a = CAT-1 delete + CAT-2 generalize (+ reconcile the CAT-2
occurrence-key inconsistency). CAT-3 (small: `effective_count`->`W`, `max_footprint`'s non-default
use) + CAT-3b (compile-time factorizer knobs `node_level_placement`/`order_aware`) are Phase-4/O2
scope -- O2 supersedes their function; only `effective_count` is literally replaced. CAT-4
router-seam duplication resolves as the seam matures. CAT-5 (incl. `is_volatile`, `min_repeats`)
stays -- these are cache semantics / CSE-cacheability, not spill heuristics.
