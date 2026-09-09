# Kramers TRS canonicalizer — implementation plan

> For agentic workers: execute task-by-task, test-first, in
> `~/code/sequant` (round 2), building tests in `cmake-build-tests`
> (`ninja -C cmake-build-tests unit_tests-sequant`, run
> `cmake-build-tests/tests/unit/unit_tests-sequant "<tag>"`). Commit per
> task. Design: `doc/dev/specs/2026-09-02-kramers-trs-canonicalizer-design.md`.

**Goal:** the canonicalizer knows the Kramers time-reversal identity, so
only Kramers-canonical configs (and their conj+phase variants) survive
canonicalization, at the eval leaf boundary and in the symbolic equations.

**Architecture:** a `KramersSymmetry` tensor attribute; a registry-level
flavor-partner query; a single-tensor fold in `TensorBlockCanonicalizer`
whose byproduct is `CanonTransform{phase, conj}`; a network fold under
`CanonicalizeOptions::fold_kramers` with flavor-blind colours, per-component
orientation, and a `kramers_flipped_tensors`+phase byproduct.

---

### Task 1: `KramersSymmetry` attribute plumbing (DONE)

Files: `SeQuant/core/attr.hpp`, `SeQuant/core/expressions/tensor.hpp`,
`SeQuant/core/expressions/abstract_tensor.hpp`,
`SeQuant/core/io/serialization/v1/ast_conversions.hpp` (+ writer),
`tests/unit/test_conjugation.cpp` (new TEST_CASE `kramers_symmetry_attribute`).

- [ ] Test: construct `Tensor(L"C", bra{a↑_1}, ket{a↓_2<...>}, Symmetry::Nonsymm,
      BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm, KramersSymmetry::TimeReversal)`;
      REQUIRE `kramers_symmetry() == TimeReversal`; `with_slots` copies it;
      two tensors differing only in the attribute are not equal and hash
      differently; serialization roundtrip `C{a↑_1;a↓_2}:N-C-N-T`.
- [ ] Run: expect compile failure (no enum). Add the enum (`attr.hpp:37`
      area), the member near `tensor.hpp:890`, accessor near `:679`, ctor
      parameters (both main ctors and the Hermiticity ctors, default
      `Nonsymm`), `with_slots` (`:824`), `_hash_value` (`:917`),
      `static_equal`/`static_less_than` (`:940-967`), `AbstractTensor::_kramers_symmetry()`
      + free `kramers_symmetry(const AbstractTensor&)` (`abstract_tensor.hpp:416`),
      serialization letter `T`/`N` as optional 4th field.
- [ ] Run: pass. Commit `core: KramersSymmetry tensor attribute`.

### Task 2: flavor partner in the registry (DONE)

Files: `SeQuant/core/index_space_registry.hpp/.cpp`, `SeQuant/core/index.hpp`,
`SeQuant/domain/mbpt/convention.cpp`, `tests/unit/test_index.cpp` (new
TEST_CASE `kramers_partner`).

- [ ] Test (mbpt context): `registry->kramers_partner(space(L"a↑"))` is
      `a↓` and back; `kramers_partner(space(L"Κ"))` is `nullopt`;
      `Index(L"a↓_3", protos{i↑_1,i↓_2}).kramers_flipped(reg)` is
      `a↑_3<i↓_1,i↑_2>`.
- [ ] Implement: `IndexSpaceRegistry::add_kramers_partners(a, b)` + query;
      populate in `convention.cpp:56` where the spin clones are made;
      `Index::kramers_flipped(const IndexSpaceRegistry&)` (ordinal kept,
      protos mapped recursively, unflavored protos kept).
- [ ] Run: pass. Commit `core: Kramers partner spaces in the registry`.

### Task 3: single-tensor Kramers fold + eval leaf boundary (DONE)

Files: `SeQuant/core/tensor_canonicalizer.hpp/.cpp`,
`SeQuant/core/eval/eval_expr.cpp`, `tests/unit/test_canonicalize.cpp`
(`kramers_block_fold`), `tests/unit/test_eval_expr.cpp`
(`kramers_leaf_slot_identity`).

- [ ] Test 1 (`kramers_block_fold`): for `C{a↓_1;a↑_2<i↑_1,i↓_1>}` with
      TimeReversal, `TensorBlockCanonicalizer{}.apply` yields the spelling
      `C^*{a↑_1;a↓_2<i↓_1,i↑_1>}` (first flavored slot up), phase −1;
      `C{a↓_1;a↓_2<..>}` → `C^*{a↑_1;a↑_2<..>}`, phase +1; a tensor whose
      first flavored slot is up is untouched (phase +1); a Nonsymm-Kramers
      tensor is untouched; DF aux slot never flips.
- [ ] Test 2 (`kramers_leaf_slot_identity`): `EvalExpr{C↓↑}` and
      `EvalExpr{C↑↓}` have equal `hash_value()`; the former's
      `canon_transform()` is `{phase −1, conj}`, the latter trivial;
      `C↓↓` vs `C↑↑`: `{+1, conj}`; the product `f{a↓;i↑} C{...}` leaves
      compose as expected (mirror `leaf_transform_channels`).
- [ ] Implement `kramers_foldable` (TimeReversal && `_is_cnumber()` &&
      !pinned), `canonicalize_kramers` (flip all flavored slots via
      `kramers_flipped`, toggle marker, phase), call from
      `TensorBlockCanonicalizer::apply` after `canonicalize_braket`; in
      `EvalExpr::EvalExpr(Tensor const&)` (`eval_expr.cpp:206-227`) multiply
      the phase and set `conj` when the fold fired.
- [ ] Run: pass. Commit `canonicalize: single-tensor Kramers fold; eval leaves share up-row slots`.

### Task 4: network fold under `fold_kramers` (DONE — component pre-pass)

Files: `SeQuant/core/options.hpp/.cpp`, `SeQuant/core/tensor_network/v3.hpp/.cpp`,
`tests/unit/test_canonicalize.cpp` (`kramers_network_fold`).

The plan's original Test 1 (`C{a↓_1;a↑_2<P>}` with `a↑_2` external) was
invalid under the whole-tensor identity: a tensor holding an external
flavored index cannot flip at all. Replaced by:

- [x] free component: `g{i↓_1,i↓_2;a↓_1,a↓_2} t{a↓_1,a↓_2;i↓_1,i↓_2}`
      canonicalizes (fold on) onto the all-up spelling with both tensors
      marked, scalar unchanged; mixed spelling picks the up-first
      orientation.
- [x] pinned component: `f{a↓_2;a↑_1} t{i↑_1;a↓_2}` (externals `a↑_1`,
      `i↑_1`) does not flip.
- [x] default options do not fold; the fold is idempotent.
- [x] `canonicalize_slots({.fold_kramers = true})`: down spelling hashes
      equal to the MARKED up spelling, reports both ordinals in
      `kramers_flipped_tensors`; without the fold the hashes differ.
- [x] Implemented `TensorNetworkV3::kramers_orient` (see design §4);
      `CanonicalizeOptions::FoldKramers` (+`copy_and_set`),
      `CanonicalizeSlotsOptions::fold_kramers`; no vertex-painter change.
- [ ] Follow-up: serialization letter for `KramersSymmetry` (Task 1 left it
      out; tests build tensors programmatically).

### Task 5: mbpt / emission integration (DONE)

Files: `SeQuant/domain/mbpt/spinor.hpp/.cpp`, `rules/csv.cpp`, `rules/df.cpp`,
`SeQuant/core/expressions/tensor.hpp` (`set_kramers_symmetry`),
`tests/unit/test_spinor.cpp` (`kramers_symmetry_propagation`).

- [x] Test: after `closed_shell_kramers_trace` + `csv_transform` +
      `density_fit` every leaf reports `TimeReversal`; with `fold_kramers`
      the canonicalized term never has MORE down-first leaves than the input
      and canonicalization is idempotent; `kramers_internal_rebase` returns
      its input when the context folds. (The plan's "no down-first leaf
      survives" was too strong: a component such as `g{↑↓;↑↓} t{↓↑;↓↑}` has
      one down-first leaf in either orientation.)
- [x] Implement: `mark_kramers_symmetric` (deep copy + stamp) applied to
      the trace outputs (energy and CC blocks); `csv_transform` and
      `density_fit` inherit the source tensor's attribute; rebase
      early-returns under a folding context. No `OpMaker` change: marking
      at the trace output keeps non-Kramers hashes untouched.

### Task 4b: Kramers-aware braket orientation (DONE, 0ade69912)

Found by Task 6: see design §4 "Braket orientation is Kramers-aware".
Tests: `kramers_block_fold` (orientation preference), eval leaf for a
Hermitian C, dch energy term 0 reproduction in `kramers_symmetry_propagation`.

### Task 4c: symmetry-invariant orientation verdicts (DONE)

See design §4 "Orientation verdicts are symmetry-invariant". Regression:
the three dch energy twin pairs in `kramers_symmetry_propagation`.

### Task 6: MPQC measurement (dch) — DONE

- [x] `MPQC_CCK_KRAMERS_FOLD=1` (opt-in; options passed explicitly, see
      design §5) on dch default: PNS-MP2 −1.04169026891 (band), 10
      iterations, energy equation 36 terms, residual blocks 48 each.
- [x] Energy conjugate-pair fold: 36 → 20 (17 pairs + 2 self-conjugate),
      vs 24 on the rebase path. Pairing diagnostic:
      `MPQC_CCK_TRS_FOLD=1 MPQC_CCK_TRS_FOLD_TRACE=1` prints per-term
      canonical hashes of the term and its conjugate.
- [x] Default flipped (fold on for the Kramers-CSV path, MPQC
      `MPQC_CCK_NO_KRAMERS_FOLD` opts out); `kramers_internal_rebase`
      retired (SeQuant 159739551). Default-path dch: −1.0416902697, 10 it.
- [x] Residual census (exact, `[cck-eqs/census]`): blocks unchanged by the
      fold (26/46/82/50/70 non-canonical leaves of 256) — each residual term
      is one component anchored to the externals; the leftover down leaves
      are internal summations coupled to up externals → T19 is the lever.
      Twin fold on the self-conjugate block: still 0/48.
- [x] T20 (wrapped-summand eval): RESOLVED — Re/Im wrapper factors are
      transparent to the optimizer and the wrapper's inner batch axes are
      re-keyed under the summand (650b9ba92); dch folded energy 1.70 GB
      peak (was 14 GB), MPQC CSV energy fold default-on.
- [x] T19 layer 1 (bc449ff28): lazy {phase, conj, perm} views on
      `ResultTensorTA` (TA `.conj()`/scaling are lazy expressions; flat
      tensors only — `ResultTensorOfTensorTA::apply_transform` still
      materializes one ToT copy per application).
- [x] T19 layer 2 (ad05f6ce8): under `fold_kramers_eval_leaves` a
      Kramers-noncanonical leaf stores the up-row spelling in `expr()`
      (what a provider fetches) and the as-written labels in
      `canon_indices()`; `kramers_folded()`, `denoted_spelling()`.
      LANDMINE fixed on the way (7baa8bfa7, also on PR-2):
      `CanonicalizeOptions::operator==` compared only `method`, so a scoped
      context differing only in a fold flag was a silent no-op
      (`set_scoped_default_context` skips equal contexts).
      MPQC wiring: `EvalContext::fold_kramers_eval_leaves` scoped over
      optimize+binarize on the complex Kramers-CSV path
      (`MPQC_CCK_NO_KRAMERS_LEAF_FOLD` opts out); the dense f/g leaf
      builder and the CSV t-leaf reconstruction map the as-written
      annotation to the folded spelling by flipping ↑/↓. dch: iterations
      1-8 identical to the unfolded run to 1e-11, tail differs at the
      1e-9 convergence-noise level (−1.04169026663 at residual 1.6e-10),
      3.3 s/it unchanged, 1.70 GB peak; the (⇓,⇓) C block is no longer
      requested.
- [x] T19 layer 3 (MPQC): the PNS provider emits only the K_up column of
      the Kramers C blocks (2 per rank); `CSV::coefficients(rank, K_dn,
      row)` throws and `eval_csv` derives a ⇓ request as the 𝒯 image of
      the stored block ((⇓,↑) = −conj((⇑,↓)), (⇓,↓) = +conj((⇑,↑))). A
      mixed C never folds at the leaf (it is `BraKetSymmetry::Conjugate`,
      so `kramers_flavor_key` sorts its bundles and the flip ties), which
      is fine: those leaves arrive through the braket-swap orientation and
      are served from the stored column. dch: fold off reproduces the old
      trajectory to 1e-12 (−1.04169026953), fold on converges further
      (−1.04169026662, residual 2.6e-10).
- [x] T21 (numerical antisymmetrization for PNS-MP1/2, MPQC only): the
      PNS-CCD partial-Â route (`kramers_partial_A`: Â expanded over
      mixed-flavour external groups, same-flavour groups antisymmetrized
      numerically within the block) applies unchanged to the MP1 residual
      and is now its default. dch PNS-MP1: residual 48x5 -> 28/44/28/32/48
      terms, iterations identical to 1e-12, -1.04169026965 in band.

## Regression found and fixed on 2026-09-02 (evening)

- **Symptom** (caught by the user from the iteration table): HSeOH PNS-CCD
  iteration 1 took 12 s and later iterations 82 s (reference: 127 s then
  3.5 s), with wrong energies; HSeOH PNS-MP1 diverged (iteration-1 energy
  -0.2221 vs -0.299572965188). dch MP1 stayed in band, which is why every
  "certification" above missed it.
- **Two independent defects.** (1) MPQC af48aa6cf6 expanded and flattened
  the CSV flavor sums for residuals too (R2 2201 nested -> 30324 flat
  terms); fixed by confining the expansion and the per-term network
  canonicalize to the fully contracted energy. (2) SeQuant PR-2 (merged
  1bfa6dd42 + 74e596c5c): the tensor-of-tensors leaf branch took its
  antisymmetric reorder phase from `canonicalize_slots`, which only labels
  and never reorders a ToT tensor, so the stored spelling (what the
  provider serves) stayed as written while the retrieval transform negated
  it; 4 of 22 MP1 energy terms flipped sign. Localized by a two-tree bisect
  (MPQC pinned at 5891c30019, SeQuant PR-1 good / PR-2 bad), then a
  per-term energy diff. Fixed in c66a2aed2 (ToT leaves block-canonicalize
  the spelling in place and put antisymmetric bundles into the labeling's
  canonical order with the parity as the phase); mirrored on PR-2.
- **Certification rule from now on**: HSeOH PNS-MP1 (`pnsmp1_loc1e-5.json`,
  iteration 1 = -0.299572965188) and HSeOH PNS-CCD (`pnsccd_3it.json`,
  iterations 1-3 = -0.299572965188 / -0.320010602548 / -0.328438791575,
  2201 terms, 3.5 s/it) in addition to dch, for every eval/trace change.

## Review fixes, 2026-09-03

Two reviews (a SeQuant reviewer on the four T19/regression commits, an
MPQC multi-angle review) plus the follow-up investigation of the ToT leaf
phase led to these changes:

- **Slot order applied by the labeling itself.** `CanonicalizeSlotsOptions::
  apply_slot_order` makes `canonicalize_slots` permute every (anti)symmetric
  bra/ket bundle into its canonical vertex order in place; the parity comes
  from the same sort (`sort_then_replace_by_ordinals`) that defines
  `metadata.phase`, so the two cannot disagree (the hand-rolled reorder in
  the ToT leaf branch sorted by the named-index order instead and asserted
  agreement). The ToT leaf branch now just asks for it.
- **Child phases hoist into their parents.** A leaf's reorder phase is a
  transform (hash-blind), so the two spellings share one slot -- but the
  parent hashed the network of canonical child spellings while computing
  its value from the DENOTED hand-ups, so `g * t{a2,a1}` and `g * t{a1,a2}`
  shared a slot holding sign-different values. Products now fold the
  children's phases into their own transform (the spec's "multiplicatively
  hoistable" rule, previously only applied to the scalar wrap); sums hoist a
  uniform phase and salt a mixed one (`CanonTransform::phase_salt`); the
  Re/Im wrapper hoists the inner phase.
- **Sum slot identity covered only the summands before the last.** The
  prefix hashes came from `inits()` -- a lazily sliced view driven by a
  stateful lambda -- evaluated by random access in `make_sum`, so the
  n-summand node hashed n-1 summands: `A + B` and `A + C` had one hash.
  Pre-existing (5b64bdd27); `imed_hashes` is eager now, `inits` removed.
  Test `sum_slot_identity_covers_every_summand`.
- **Denoted spelling.** `EvalExpr::denoted_expr()` (public; replaces the
  file-local `denoted_spelling`) undoes the normalization channels in
  reverse order and takes the marker out for every channel that came with
  a conj (a swap, the Kramers flip): a Conjugate leaf written in its
  non-canonical orientation denotes as written, unmarked (it was marked).
- **Cleanups.** One `normalize_leaf` helper serves the flat and ToT
  branches; `kramers_flip(index_vector&, isr)` next to `kramers_flipped`;
  `ResultTensorTA::logical_array` materializes through
  `ensure_materialized` (memoized); `CanonicalizeOptions`/`SimplifyOptions`
  equality defaulted; docs on `Result::value_` (not thread-safe), `is_view`.
- **A sum's layout is part of its slot identity.** A sum hands up its
  FIRST summand's layout, so A + B and B + A must not share a slot; the
  prefix hashes are order-sensitive (first tried as a label salt, c196d43c7,
  which lost the sharing of relabeled sub-sums: 5208 -> 7955 distinct
  intermediates and a 20 GB cache high-water on HSeOH PNS-CCD).
- **The denoted marker rule of 88aaacda2 was wrong and is reverted**: the
  reviewer's "take the marker out for every channel that came with a conj"
  makes a swapped Hermitian leaf denote unmarked, which is value-correct in
  isolation but removes the conj color from the PARENT network's graph;
  HSeOH PNS-CCD then served a cached intermediate in the wrong layout in
  iteration 2 (TA range-congruence abort) while every unit suite stayed
  green. Localized by a toggle bisect (env switches for the three eval
  changes, one build, five runs) -- the denoted spelling is an identity
  convention of the parent network, not a value statement; documented on
  `EvalExpr::denoted_expr`.
- **apply_slot_order applies the NAMED-index canonical order** (4872f372e).
  The raw canonical vertex ordinals order same-color cells by the color
  hash, so a mixed-flavor antisymmetric bundle came out as `a↓,a↑` (288
  such t spellings in the HSeOH PNS-CCD trace); the PNS provider stores
  representative Kramers configurations only, and the mismatched
  configuration request aborted iteration 2 (TA range congruence in a ToT
  contraction). In apply mode the order applied -- and the phase reported
  -- is the named-index canonical order (space-major, vertex ordinal within
  a group), the order the hand-rolled reorder had used.
- Deferred: lazy views for `ResultTensorOfTensorTA` (T19 layer 1 does not
  reach the PNS hot path yet); the TA `eval_with_tiledarray/real/
  summation` test fails on this branch since before these commits
  (verified at 78841916f; PR-2's copy of the test file passes) -- the
  round-2 test file's NonHermitian deserialization plus the
  orientation-canonicalizing `tensor_to_key` fixture, to be reconciled.
