# Slot-symmetry deduction via TN canonicalization (design note)

**Date:** 2026-07-01
**Status:** design record; supersedes the deduction approach in
`2026-06-28-general-symmetry-sector-storage-design.md` (Phase 0). The Phase-0
`SlotSymmetry` carrier and its out-of-band (1a) wiring stand; the *deduction
algorithm* is what this note replaces.
**Repos:** SeQuant (`core/eval/slot_symmetry.*`, `core/tensor_network/v3.*`),
consumed by mpqc4 CCk/CSV.

## Problem with the Phase-0 hand-rolled rules

Phase 0 deduces an intermediate's permutational symmetry by propagating
*declared* leaf symmetries through a **bijective index -> slot** trace
(`deduce_slot_symmetry`). That model is a limited approximation with two
failure modes:

1. **Dependent indices break bijectivity (unsound in products).** An index can
   occupy multiple slots at once, e.g. `F{a1<i1,i2>; i1}` where `i1` is both a
   ket slot and a proto-index of `a1`. The flat trace treats slots as atomic
   labels and never looks inside proto-lists, so a contraction on such an index
   silently couples slots the trace thinks are independent. In a product this
   can yield a **false positive** (a claimed symmetry the tensor lacks) --
   precisely on the proto-indexed CSV/PNO intermediates the feature targets.
   (The earlier `full_label()` fix removed only a *label collision*; it does not
   fix the incidence problem.)

2. **Emergent symmetry from topologically-equivalent tensors (incomplete).**
   `A{a1;i1} A{a2;i2}` has column symmetry `{a1,i1} <-> {a2,i2}` because the two
   `A` factors are interchangeable -- a property of the product's *graph*, not
   of any declared leaf symmetry. Declared-symmetry propagation is structurally
   blind to it (a false negative; safe but incomplete, and arguably the dominant
   origin of column/particle symmetry in built intermediates).

Root cause: symmetry of a tensor/TN is a **graph-automorphism** property, not a
declared-attribute-propagation property.

## Correct deduction: reuse the TN canonicalizer (graph-theoretic)

SeQuant already models everything needed in the `TensorNetworkV3` colored graph:
`TensorCore` colors encode tensor identity (=> identical factors give isomorphic
subgraphs => emergent symmetry), `SPBundle` vertices encode proto-index
dependencies (=> dependent indices respected), and `TensorBra`/`TensorKet`/
`TensorBraKet` distinguish bra/ket/column roles. And `canonicalize_slots`
returns a `SlotCanonicalizationMetadata` carrying a `phase` (+/-1) for the slot
permutation applied (`v3.hpp:283`).

**Deduction primitive (canonicalize-and-compare probe):** to test whether a
candidate external-slot permutation `pi` is a symmetry, canonicalize the tensor
and the `pi`-permuted tensor; if they reach the **same canonical slot-form**,
`pi` is a symmetry and the **`phase` difference is the sign**. Run a small,
structurally-motivated candidate set (column swaps, adjacent bra/ket swaps,
identical-factor exchanges); each hit is a generator + sign. Because protos and
tensor-identity are already in the graph, this subsumes **both** failure modes
above for free. (`graph->find_automorphisms(...)` -- already used in
`wick.impl.hpp:943` -- gives the whole group directly if we ever want it instead
of probing.)

**`canonicalize_slots` is sufficient.** No group-theory engine is needed for
deduction beyond what the graph canonicalizer already provides.

## Storage-sector canonicalization: NOT Butler-Portugal

Two distinct "canonicalizations" must not be conflated:

- **Symbolic double-coset** (canonicalize a tensor *expression* `S.g.D` to
  decide expression equality) -- this is Butler-Portugal / xPerm. It is the very
  thing SeQuant's TN/bliss canonicalizer was designed to **replace**; it scales
  poorly and is mono-term-only. **We never do this for storage. BP is rejected.**
- **Numeric storage-sector** (map a concrete integer index tuple to its
  canonical representative under the slot-symmetry group `G`, + sign) -- a
  per-tile runtime op.

For the symmetries we deduce (products of symmetric groups on disjoint
slot/column blocks = Young-subgroup type), the storage-canonical form is just
**sort within each symmetric block; sign = sort parity**. Classic triangular
storage, O(n log n) per tuple, no group machinery. For a general (non-block)
mono-term `G` (e.g. a pure cyclic symmetry), the canonical form is **orbit-min:
enumerate the small explicit `G` and take the minimum** -- still a for-loop
(|G| ~ 2..few hundred), still not BP/BSGS.

## Target descriptor

A **signed permutation group on external slots** (a set of generators + a `+/-1`
sign character), optionally with a **per-generator conjugation bit `kappa in
{id, *}`** for hermiticity-type symmetries (`T = conj(T o sigma)`). This is a
1-dimensional representation twisted by an optional antilinear (conjugation)
character -- i.e. the sign is a linear character `G -> {+/-1}`, and `kappa` a
second `Z2` character acting antilinearly.

- The Phase-0 bespoke 3-list (`column_groups`/`bra_groups`/`ket_groups`) is a
  **lossy projection** of this: it can represent Young-subgroup-type symmetries
  (the common case) but not general groups (e.g. cyclic-only) or conjugation.
- Add `kappa` from the start: hermiticity pervades our integrals, it is a cheap
  second character, and it is the OQ-3 gap libPerm does not cover.

## Scope: mono-term only

- **In scope (storage-actionable): mono-term symmetry** = the tensor spans a
  1-dim (`+/-1`) subrep of `G`. Admits a canonical fundamental domain =>
  triangular/canonical-sector storage.
- **Out of scope: multi-term / multidimensional-irrep symmetry** (Young/Bianchi;
  relations `sum_sigma c_sigma T o sigma = 0`). Relevant to spin-adapted
  higher-order CC (genuine `S_n` irrep multiplicity), but it reduces **rank**,
  not the index box -- a different problem, deliberately excluded here.

## libPerm

Not adopted (decided). We need neither its `canonicalize()` (block-sort /
orbit-min suffices for storage) nor group intersection at storage time. If the
descriptor ever needs to hold a general (non-Young) signed group, a libPerm-
style container becomes relevant -- but the group-intersection gap (for Sum-node
deduction) and the missing conjugation character remain its limitations.

## Interim: keep Phase 0 as a sound conservative fast-path

Until the graph-canonicalization deducer lands (Phase 0.5), the current rules
ship **guarded to be never-wrong**. The guiding principle: **guard where
soundness cannot be proven; do not guard where the result is provably only
incomplete.** A deduced group `G_claimed` is safe for a storage consumer iff
`G_claimed` is a *subset* of the true symmetry `G_true` (finer orbits => it
under-compresses => correct); only a *superset* (claiming a non-symmetry)
loses data.

- **Guard (proto-indexed externals) -- KEPT.** `deduce_slot_symmetry` returns
  empty if any participating operand/result external carries proto-indices. The
  flat index->slot trace's bijectivity precondition fails under index
  dependencies, and we cannot cheaply prove the trace stays sound, so we decline
  (potential false positive avoided). The deducer simply declines on CSV/PNO
  until Phase 0.5.
- **Repeated identical factors -- NOT guarded (deliberately).** Failure mode 2 is
  a provable false *negative*: the flat rules miss the emergent exchange
  symmetry, but the 4-tuple supplier key prevents any wrong cross-copy merge, so
  `G_claimed` is always a subset of `G_true`. A "same-core -> bail" guard would
  drop correct claims (e.g. a genuinely inherited column group in a PPL
  contraction of identical factors) with no soundness benefit. It is documented
  by a test asserting the current (sound, incomplete) behavior; the emergent
  symmetry is recovered by the Phase-0.5 deducer, not by a guard.

The net effect: every non-empty descriptor is a subset of the true symmetry, so
a consumer may safely trust it.

## Open questions for Phase 0.5

- **External coloring for symmetry:** canonicalization pins named (external)
  indices by identity; symmetry deduction instead needs externals colored **by
  type** (space x bra/ket x aux) so `Aut` can permute them. The external/internal
  distinction is configurable at the call site via `named_indices` -- confirm a
  by-type external coloring is expressible through the same color hook.
- **Sign extraction:** cleanest is for each leaf's `Antisymm` bundle to
  contribute its induced parity as a probe permutation is applied; confirm this
  is recoverable from the canonical `phase` or needs per-leaf bundle identity.
- **Conjugation (`kappa`):** how to model bra<->ket-with-conjugation in the graph
  coloring / phase so hermiticity is deduced, not just permutation symmetry.
- **Sum nodes:** the sum's symmetry is the **intersection** of the summands'
  groups. Compute via the graph of the combined expression, or as an explicit
  group intersection (the one place a real group-intersection routine would be
  needed -- and libPerm lacks it).
