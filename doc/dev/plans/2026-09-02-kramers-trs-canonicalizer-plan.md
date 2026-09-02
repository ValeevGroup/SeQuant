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

### Task 1: `KramersSymmetry` attribute plumbing

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

### Task 2: flavor partner in the registry

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

### Task 3: single-tensor Kramers fold + eval leaf boundary

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

### Task 4: network fold under `fold_kramers`

Files: `SeQuant/core/options.hpp/.cpp`, `SeQuant/core/tensor_network/v3.hpp/.cpp`,
`SeQuant/core/tensor_network/vertex_painter.cpp`,
`SeQuant/core/expressions/expr_algorithms.cpp` (respell), tests in
`tests/unit/test_tensor_network.cpp` (`SECTION("kramers fold")` in
`tensor_network_shared`) and `tests/unit/test_canonicalize.cpp`
(`kramers_network_fold`).

- [ ] Test 1: two products differing by flipping the flavors of an
      internal component (e.g. `g{a↑_1;i↑_1;Κ_1} C{a↑_1;a↑_2<P>}` vs
      `g{a↓_1;i↓_1;Κ_1} C{a↓_1;a↑_2<P>}` with the external `a↑_2<P>`):
      under `fold_kramers` both `canonicalize_slots` runs give equal
      `hash_value()`; exactly one reports both tensors in
      `kramers_flipped_tensors` with phase per the rule.
- [ ] Test 2: canonicalize() of the down-spelled product returns the
      up-spelled product with conj markers and the phase in the scalar;
      canonicalize is idempotent; a component touching a named down index
      keeps its orientation.
- [ ] Test 3: `tn_slots_determinism` (#57) still passes with
      `fold_kramers` on: `C·C*` vs `C*·C` remain distinct.
- [ ] Implement: `fold_kramers` in `CanonicalizeOptions` (+`copy_and_set`),
      `CanonicalizeSlotsOptions` (`v3.hpp:342`), `CreateGraphOptions`
      (`v3.hpp:406`), forwarded into `canonicalize_graph` (`v3.cpp:575`);
      flavor-blind index colour (`vertex_painter.cpp:104-121`, mask the
      partner-space distinction via the registry) and core-colour
      normalization (`:29-32` pattern); component orientation after the
      canonical labeling (union-find over flavored dummies, externals fix
      their component, criterion = fewer down-first leaves then
      fingerprint); byproduct `kramers_flipped_tensors` + `phase` next to
      `conjugated_tensors` (`v3.hpp:296`, populated in the `v3.cpp:966-1012`
      walk); respell in `canonicalize` where `conjugated_tensors` is applied.
- [ ] Run: pass; run the full `[canonicalize]`, `[tensor_network]`,
      `[conjugation]`, `[eval_expr]`, `[spinor]` tags. Commit
      `canonicalize: network Kramers fold (flavor-blind colours, component orientation, conj+phase byproduct)`.

### Task 5: mbpt / emission integration

Files: `SeQuant/domain/mbpt/op.cpp`, `op_registry.hpp/.cpp`,
`rules/csv.cpp`, `rules/df.cpp`, `spinor.cpp`, `tests/unit/test_spinor.cpp`.

- [ ] Test: after `closed_shell_kramers_CC_trace` + `csv_transform` +
      `density_fit` in a Kramers context, every C/g/f/t leaf reports
      `TimeReversal`; with `fold_kramers` on, canonicalize of the energy
      expression leaves no tensor whose first flavored slot is down
      (census helper in the test); `kramers_internal_rebase` on the result
      is a no-op (returns an equal expression).
- [ ] Implement: registry default + `OpMaker` assignment; csv/df propagate
      the attribute; `kramers_internal_rebase` early-returns when the
      context has `fold_kramers` (retire after the MPQC measurement).
- [ ] Run: pass. Commit `mbpt: Kramers-symmetric ops; csv/df propagate KramersSymmetry`.

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
