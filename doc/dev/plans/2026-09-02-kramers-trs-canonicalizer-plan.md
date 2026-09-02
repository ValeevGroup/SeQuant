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

### Task 4b: Kramers-aware braket orientation (DONE, a62e52365)

Found by Task 6: see design §4 "Braket orientation is Kramers-aware".
Tests: `kramers_block_fold` (orientation preference), eval leaf for a
Hermitian C, dch energy term 0 reproduction in `kramers_symmetry_propagation`.

### Task 4c: symmetry-invariant orientation verdicts (DONE)

See design §4 "Orientation verdicts are symmetry-invariant". Regression:
the three dch energy twin pairs in `kramers_symmetry_propagation`.

### Task 6: MPQC measurement (dch)

Files: `mpqc4` `cc/sequant.cpp` (context option), runs in
`~/code/runs-mpqc/pr2-eval-smoke/`.

- [ ] Enable `fold_kramers` in MPQC's SeQuant context for the Kramers
      path; rebuild; run dch default: energy in the −1.0416902693 band,
      iteration count unchanged.
- [ ] Census (`twin_census.py`): expect zero down-first leaves in internal
      components; report per-block counts. Energy fold: expect 36 → 18.
- [ ] Twin fold (`MPQC_CCK_TWIN_FOLD=1`, block 4): report pairs; if still
      0, run the numerical term study extension (`t' == sign·conj(perm(t))`)
      to settle value-level pairing. Record results in the design doc's
      "Scope and expectations" and in the PR-2 plan roadmap.
- [ ] Commit the MPQC context change; update memory.
