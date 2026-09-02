# Lazy-Conj Eval Core (Conjugation PR 2, Plan A) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the structural `EvalOp::Adjoint` conjugation model with a
retrieval-time canonical transform `{phase, conj, braket_swap}`, giving
conj-aware cache/CSE identity (uniform-conj hoisting) and serving every
marked leaf the symbolic layer can produce.

**Architecture:** Spec at `doc/dev/specs/2026-09-01-lazy-conj-eval-design.md`.
Slot identity = canonical spelling (transform excluded); structural identity
salts child hashes with conj/swap only, with uniform-conj hoisting to the
parent transform. Plan B (separate, queued) covers Re/Im eval nodes, the
FoldConjugatePairs default flip, export conj emission, the exhaustive
symbolic-to-eval test matrix, and the MPQC smoke.

**Tech Stack:** C++20, Catch2 v3 (`build-tests-debug`, SEQUANT_ASSERT=THROW),
TiledArray/BTAS/tapp eval backends.

**Branch:** `kshitij/feature/conjugation-eval` (stacked on PR #602).
Build+test loop used by every task:
```bash
ninja -C build-tests-debug unit_tests-sequant
./build-tests-debug/tests/unit/unit_tests-sequant "<filter>"
```

---

### Task 1: `CanonTransform` value type

**Files:**
- Modify: `SeQuant/core/eval/eval_expr.hpp` (after the `EvalOp` enum, whose
  `Adjoint` enumerator survives until Task 9)
- Test: `tests/unit/test_eval_expr.cpp`

- [ ] **Step 1: Write the failing test** (append to `tests/unit/test_eval_expr.cpp`):

```cpp
TEST_CASE("canon_transform_algebra", "[EvalExpr][conj-transform]") {
  using sequant::CanonTransform;
  CanonTransform id{};
  REQUIRE(id.trivial());
  REQUIRE(id.phase == 1);
  REQUIRE_FALSE(id.conj);
  REQUIRE_FALSE(id.braket_swap);

  CanonTransform c{.phase = 1, .conj = true, .braket_swap = false};
  CanonTransform s{.phase = -1, .conj = false, .braket_swap = true};
  REQUIRE_FALSE(c.trivial());

  // composition: phases multiply, conj/swap are Z2 (xor)
  auto cs = compose(c, s);
  REQUIRE(cs.phase == -1);
  REQUIRE(cs.conj);
  REQUIRE(cs.braket_swap);
  REQUIRE(compose(c, c).trivial());  // involution
  // structural salt: conj/swap enter, phase does NOT (hoistable)
  REQUIRE(CanonTransform{.phase = -1}.structural_salt() ==
          CanonTransform{}.structural_salt());
  REQUIRE(c.structural_salt() != CanonTransform{}.structural_salt());
  REQUIRE(c.structural_salt() != s.structural_salt());
  REQUIRE(c.structural_salt() != cs.structural_salt());
}
```

- [ ] **Step 2: Run, verify FAIL** (`CanonTransform` not declared):
`./build-tests-debug/tests/unit/unit_tests-sequant canon_transform_algebra` after `ninja` — expect a compile error; that is the red state for a new type.

- [ ] **Step 3: Implement** in `SeQuant/core/eval/eval_expr.hpp`, directly after the `EvalOp` enum's closing brace:

```cpp
///
/// \brief Canonicalization byproduct mapping a node's CACHED canonical result
///        to the value the node denotes. Applied on retrieval
///        (see apply_canon_transform in eval.hpp); excluded from the node's
///        own (slot) hash, exactly as the former standalone canon_phase was.
///        conj/braket_swap DO enter the parent's structural hash via
///        structural_salt() — see the design spec's uniform-conj-hoists /
///        mixed-conj-salts rule.
struct CanonTransform {
  std::int8_t phase = 1;    ///< ±1 linear byproduct (antisymmetric reorder)
  bool conj = false;        ///< elementwise complex conjugation
  bool braket_swap = false; ///< bra<->ket transposition of the canonical slots

  [[nodiscard]] constexpr bool trivial() const noexcept {
    return phase == 1 && !conj && !braket_swap;
  }
  /// salt for the PARENT's hash combination: conj/swap only — phase is
  /// multiplicatively hoistable and never enters structural identity
  [[nodiscard]] constexpr std::size_t structural_salt() const noexcept {
    return (conj ? 1u : 0u) | (braket_swap ? 2u : 0u);
  }
  friend constexpr bool operator==(CanonTransform, CanonTransform) = default;
};

[[nodiscard]] constexpr CanonTransform compose(CanonTransform a,
                                               CanonTransform b) noexcept {
  return {static_cast<std::int8_t>(a.phase * b.phase), a.conj != b.conj,
          a.braket_swap != b.braket_swap};
}
```

- [ ] **Step 4: Run, verify PASS**, then run the full `eval_expr` case: both green.

- [ ] **Step 5: Commit** — `eval: CanonTransform value type (phase/conj/swap retrieval byproduct)`

---

### Task 2: `EvalExpr` carries the transform

**Files:**
- Modify: `SeQuant/core/eval/eval_expr.hpp` (member at :283, accessor at :252, ctor at :107)
- Modify: `SeQuant/core/eval/eval_expr.cpp` (ctor bodies, every 7-arg `EvalExpr{...}` construction)
- Test: `tests/unit/test_eval_expr.cpp`

- [ ] **Step 1: Failing test:**

```cpp
TEST_CASE("eval_expr_carries_canon_transform", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  Tensor t(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  EvalExpr ee{t};
  REQUIRE(ee.canon_transform().trivial());
  REQUIRE(ee.canon_phase() == ee.canon_transform().phase);  // compat accessor
}
```

- [ ] **Step 2: Run, verify FAIL** (no `canon_transform()`).

- [ ] **Step 3: Implement:**
  - hpp: replace `std::int8_t canon_phase_{1};` with `CanonTransform canon_transform_{};`;
    add `[[nodiscard]] CanonTransform canon_transform() const noexcept;`;
    keep `canon_phase()` returning `canon_transform_.phase`.
  - Change the 7-arg ctor signature `std::int8_t phase` → `CanonTransform transform`
    (position unchanged). In the cpp: `canon_phase_ = phase ? -1 : 1;` sites become
    `canon_transform_.phase = phase ? -1 : 1;`; every construction site passing
    `1` passes `{}`, passing `ee.canon_phase()`/`canon.phase` passes
    `CanonTransform{.phase = ...}` (make_sum: `1`→`{}`; make_prod scalar branch:
    `tl->canon_phase()` → `tl->canon_transform()`; make_prod TN branch:
    `canon.phase` → `CanonTransform{.phase = canon.phase}`; binarize(Tensor)'s
    two `make_adjoint_node` calls updated in place, deleted in Task 6).
  - Compile; chase every caller the compiler flags — no other semantic change.

- [ ] **Step 4: Run `unit_tests-sequant '~[long-tests]'` — ALL green** (pure carrier change).

- [ ] **Step 5: Commit** — `eval: EvalExpr carries CanonTransform (canon_phase_ generalized)`

---

### Task 3: slot hash = canonical spelling; structural salt helper

**Files:**
- Modify: `SeQuant/core/eval/eval_expr.cpp` (`hash_terminal_tensor` ~:346; ctor ~:180)
- Test: `tests/unit/test_eval_expr.cpp`

- [ ] **Step 1: Failing test:**

```cpp
TEST_CASE("leaf_slot_identity_is_canonical_spelling",
          "[EvalExpr][conj-transform]") {
  using namespace sequant;
  Tensor t(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  Tensor ts = t;
  ts.conjugate();
  // t and t^* SHARE one slot; the conj rides in the transform
  REQUIRE(EvalExpr{t}.hash_value() == EvalExpr{ts}.hash_value());
  REQUIRE(EvalExpr{ts}.canon_transform().conj);
  REQUIRE_FALSE(EvalExpr{t}.canon_transform().conj);
}
```

- [ ] **Step 2: Run, verify FAIL** (PR-1 salts the Nonsymm marker into the hash, and no transform is derived yet — both REQUIREs red).

- [ ] **Step 3: Implement:** in `hash_terminal_tensor` delete the Nonsymm
  marker `hash::combine`; in the flat-leaf ctor set
  `if (t.conjugated() && t.braket_symmetry() == BraKetSymmetry::Nonsymm) { canon_transform_.conj = true; }`
  and hash the unmarked spelling (copy, `conjugate()`, hash) — the marker no
  longer reaches the hash by any route. (Conjugate/Symm markers get their
  transforms in Task 4; keep their current behavior compiling.)

- [ ] **Step 4: Run — new test PASSES; `starred non-Conjugate leaves` section
  now FAILS at its hash-distinct assertion. Update that PR-1 assertion in the
  same change** (it pinned the gate model): replace

```cpp
    REQUIRE(EvalExpr{t}.hash_value() != EvalExpr{t_star}.hash_value());
```
with
```cpp
    // PR-2 slot identity: one slot, conj in the transform
    REQUIRE(EvalExpr{t}.hash_value() == EvalExpr{t_star}.hash_value());
    REQUIRE(EvalExpr{t_star}.canon_transform().conj);
```

- [ ] **Step 5: Full suite green; commit** — `eval: leaf slot identity = canonical spelling; Nonsymm marker moves to the transform`

---

### Task 4: flat-leaf boundary — fold ON, all four channels become transforms

**Files:**
- Modify: `SeQuant/core/eval/eval_expr.cpp` (flat branch of `EvalExpr(Tensor)`; `binarize(Tensor)`)
- Test: `tests/unit/test_eval_expr.cpp`

- [ ] **Step 1: Failing tests:**

```cpp
TEST_CASE("leaf_transform_channels", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  // Conjugate: folded (starred+swapped) spelling -> canonical slot + {conj,swap}
  Tensor g(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor g_folded = g;
  g_folded.conjugate();
  g_folded.adjoint();  // pure swap for Conjugate
  EvalExpr eg{g}, egf{g_folded};
  REQUIRE(eg.hash_value() == egf.hash_value());   // one slot
  REQUIRE(eg.canon_transform() == egf.canon_transform() ||
          (egf.canon_transform().conj && egf.canon_transform().braket_swap));

  // '⁺' Nonsymm adjoint: label stripped, {conj, swap}
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  Tensor t_adj = t;
  t_adj.adjoint();
  EvalExpr et{t}, eta{t_adj};
  REQUIRE(et.hash_value() == eta.hash_value());
  REQUIRE(eta.canon_transform().conj);
  REQUIRE(eta.canon_transform().braket_swap);
  REQUIRE(eta.as_tensor().label() == L"t");  // bare spelling stored

  // '⁺' + marker = pure transpose {swap}
  Tensor t_adj_star = t_adj;
  t_adj_star.conjugate();
  EvalExpr etas{t_adj_star};
  REQUIRE(etas.hash_value() == et.hash_value());
  REQUIRE_FALSE(etas.canon_transform().conj);
  REQUIRE(etas.canon_transform().braket_swap);

  // Symm marker: dropped
  Tensor s(L"s", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Symm, ColumnSymmetry::Symm);
  Tensor s_star = s;
  s_star.conjugate();
  REQUIRE(EvalExpr{s_star}.canon_transform() == EvalExpr{s}.canon_transform());
}
```

- [ ] **Step 2: Run, verify FAIL** (throws/Adjoint channels still in place).

- [ ] **Step 3: Implement** the flat branch of `EvalExpr(Tensor)` as the single
  normalization point (order matters — label first, then marker, then fold):

```cpp
    auto& t0 = expr_->as<Tensor>();
    // 1. '⁺' adjoint label (Nonsymm only): strip to the bare spelling;
    //    adjoint = conj ∘ swap
    if (!t0.label().empty() && t0.label().back() == adjoint_label) {
      t0.adjoint();  // removes the label, swaps slots back
      canon_transform_ = compose(canon_transform_,
                                 {.conj = true, .braket_swap = true});
    }
    // 2. elementwise-conjugation marker
    if (t0.conjugated()) {
      switch (t0.braket_symmetry()) {
        case BraKetSymmetry::Nonsymm:   // genuine conj: transform carries it
          canon_transform_ = compose(canon_transform_, {.conj = true});
          break;
        case BraKetSymmetry::Conjugate: // folded spelling: unfold to canonical
          canon_transform_ = compose(canon_transform_,
                                     {.conj = true, .braket_swap = true});
          break;
        case BraKetSymmetry::Symm:      // identity in value
          break;
      }
      t0.conjugate();  // expr_ stores the unmarked spelling in every case
      if (t0.braket_symmetry() == BraKetSymmetry::Conjugate) t0.adjoint();
    }
    // 3. block-canonicalize WITH the fold (default flag — the eval-boundary
    //    exception is gone); a fold performed here toggles the marker, which
    //    we immediately convert to transform bits the same way:
    auto phase = TensorBlockCanonicalizer{}.apply(t0);
    canon_transform_.phase = phase ? -1 : 1;
    if (t0.conjugated()) {  // fold byproduct from canonicalize_braket
      canon_transform_ = compose(canon_transform_,
                                 {.conj = true, .braket_swap = true});
      t0.conjugate();
      t0.adjoint();
    }
    hash_value_ = hash_terminal_tensor(t0);
    canon_indices_ = t0.const_indices() | ranges::to<index_vector>;
```
  (The Task-3 Nonsymm special case is subsumed — delete it.) Note the
  invariant this establishes: **expr_ always stores the canonical unmarked
  spelling; canon_indices_ are its slots; the transform maps it to the
  denoted value.**

- [ ] **Step 4: Run the new test + `eval_expr` + `[conjugation]` — the PR-1
  Adjoint-shape sections in `eval_expr` still pass (binarize untouched so far).**

- [ ] **Step 5: Commit** — `eval: flat leaf boundary folds to canonical spelling + CanonTransform (all four channels)`

---

### Task 5: ToT-leaf boundary — fold ON, `conjugated_tensors` consumed

**Files:**
- Modify: `SeQuant/core/eval/eval_expr.cpp` (ToT branch of `EvalExpr(Tensor)`, ~:142-157)
- Test: `tests/unit/test_eval_expr.cpp` (extend `leaf_transform_channels` with a proto-indexed tensor)

- [ ] **Step 1: Failing test** (append to `leaf_transform_channels`):

```cpp
  // ToT: same channels via canonicalize_slots' conjugated_tensors report
  Tensor ct(L"C", bra{L"a_1<i_1>"}, ket{L"i_1"}, Symmetry::Nonsymm,
            BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor ct_folded = ct;
  ct_folded.conjugate();
  ct_folded.adjoint();
  EvalExpr ec{ct}, ecf{ct_folded};
  REQUIRE(ec.hash_value() == ecf.hash_value());
  REQUIRE(ecf.canon_transform().conj);
  REQUIRE(ecf.canon_transform().braket_swap);
```

- [ ] **Step 2: Run, verify FAIL.**

- [ ] **Step 3: Implement:** in the ToT branch, drop
  `.fold_conjugate_braket = false` (fold returns to the default) and derive
  the transform from the metadata: a single-tensor network reports ordinal 0
  in `md.conjugated_tensors` iff the canonical labeling spells the leaf
  swapped+starred:

```cpp
    if (!md.conjugated_tensors.empty())
      canon_transform_ = compose(canon_transform_,
                                 {.conj = true, .braket_swap = true});
    canon_transform_.phase = md.phase;
```
  Apply the same '⁺'/marker pre-normalization as the flat branch BEFORE
  building the TensorNetwork (factor steps 1-2 of Task 4 into a private
  helper `normalize_leaf_spelling(Tensor&) -> CanonTransform` used by both
  branches).

- [ ] **Step 4: Run — green.  Step 5: Commit** — `eval: ToT leaf boundary consumes conjugated_tensors into the transform`

---

### Task 6: `binarize(Tensor)` collapses; gates and Adjoint channels deleted

**Files:**
- Modify: `SeQuant/core/eval/eval_expr.cpp` (`binarize(Tensor)`, `make_adjoint_node`)
- Test: `tests/unit/test_eval_expr.cpp`

- [ ] **Step 1:** Replace the whole body of `binarize(Tensor)` with
  `return EvalExprNode{EvalExpr{t}};` and delete `make_adjoint_node`.
  Delete the two `throw std::logic_error` gates (their job is now serving) —
  `value_oriented`'s throw stays (symbolic slot-rebuild consumers).

- [ ] **Step 2:** Rewrite the now-red PR-1 IR-shape tests in
  `tests/unit/test_eval_expr.cpp`: the `"Adjoint op"` /
  `"Adjoint op in a binarized term"` / `"starred non-Conjugate leaves"`
  sections assert LEAF nodes with the expected `canon_transform()`
  (bare label, canonical slots, `{conj,swap}` for '⁺', `{conj}` for starred
  Nonsymm, throw assertions deleted). Every expected value is fixed by
  Task 4's table.

- [ ] **Step 3:** Full suite: the remaining reds must be exactly the eval-tree
  golden tests that pinned Adjoint tree shapes (`test_eval_expr`,
  `test_export`); rebless those trees (Adjoint nodes vanish, leaves gain
  transforms) with the justification recorded in the commit message.

- [ ] **Step 4: Commit** — `eval: binarize(Tensor) = plain leaf; Adjoint IR channels and PR-1 gates removed`

---

### Task 7: hoisting + salting in Sum/Product hash combination

**Files:**
- Modify: `SeQuant/core/eval/eval_expr.cpp` (`make_sum` ~:508, `make_prod` ~:560-665)
- Test: `tests/unit/test_eval_expr.cpp`

- [ ] **Step 1: Failing test:**

```cpp
TEST_CASE("conj_hoisting_structural_identity", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  auto A = parse_expr(L"A{i_1;a_1}:N-N-S", Symmetry::Nonsymm);
  auto B = parse_expr(L"B{a_1;i_1}:N-N-S", Symmetry::Nonsymm);
  auto AB  = binarize(A * B);
  auto ABc = binarize(conjugate(A * B));  // B^*·A^* after canonicalization
  // uniform conj HOISTS: same slot hash, root transform conj
  REQUIRE(AB->hash_value() == ABc->hash_value());
  REQUIRE(ABc->canon_transform().conj);
  // mixed conj SALTS: C·C^* keeps its own slot
  auto C  = parse_expr(L"C{i_1;a_1}:N-N-S", Symmetry::Nonsymm);
  auto CC  = binarize(C * C);
  auto CCc = binarize(C * conjugate(C));
  REQUIRE(CC->hash_value() != CCc->hash_value());
}
```
  (Adapt the construction idiom to the file's existing helpers — the
  assertions are the contract.)

- [ ] **Step 2: Run, verify FAIL.**

- [ ] **Step 3: Implement** one helper used by `make_sum` and `make_prod`
  wherever child hashes are combined today:

```cpp
// Structural identity of a child: slot hash + conj/swap salt — EXCEPT when
// every sibling carries the identical conj bit, which hoists to the parent
// (uniform conj = whole-node transform; see the design spec).
// Returns the parent's hoisted transform contribution.
template <typename Nodes>
CanonTransform combine_children(size_t& h, Nodes const& children) {
  bool const all_conj =
      !ranges::empty(children) &&
      ranges::all_of(children, [](auto const& n) {
        return n->canon_transform().conj;
      });
  for (auto const& n : children) {
    auto t = n->canon_transform();
    if (all_conj) t.conj = false;           // hoisted
    hash::combine(h, n->hash_value());
    hash::combine(h, t.structural_salt());  // 0 for trivial: byte-stable
  }
  return CanonTransform{.conj = all_conj};
}
```
  In `make_prod`/`make_sum` replace the bare
  `hash::combine(h, child->hash_value())` sequences with `combine_children`,
  compose the returned hoist into the node's transform, and (product only)
  conjugate the scalar when hoisting fires. **Byte-stability check:** for a
  marker-free network every salt is 0 and every hash must equal today's —
  guard with the full suite, not by construction.

- [ ] **Step 4: Run new test + full suite. Step 5: Commit** — `eval: uniform-conj hoisting + mixed-conj salting in structural hashes`

---

### Task 8: `Result::apply_transform` + backends

**Files:**
- Modify: `SeQuant/core/eval/result.hpp` (base virtual beside `adjoint` :293 and `mult_by_phase` :352; `ResultScalar` :461)
- Modify: `SeQuant/core/eval/backends/tiledarray/result.hpp` (:535 and ToT :748), `backends/btas/result.hpp` (:395), `backends/tapp/result.hpp` (:251)
- Test: `tests/unit/test_eval_ta.cpp`

- [ ] **Step 1: Failing test** (complex TA array; follow the fixture idiom of
  `TEST_CASE("eval_with_tiledarray")`):

```cpp
TEST_CASE("result_apply_transform_ta", "[eval][conj-transform]") {
  using sequant::CanonTransform;
  using sequant::eval_result;
  using sequant::ResultPtr;
  using ZArray = TA::DistArray<TA::Tensor<std::complex<double>>>;
  using ResultZ = sequant::ResultTensorTA<ZArray>;
  auto& world = TA::get_default_world();

  TA::TiledRange tr{{0, 2, 4}, {0, 3, 6}};
  ZArray R(world, tr);
  R.fill_random();
  world.gop.fence();

  ResultPtr res = eval_result<ResultZ>(R);
  std::array<std::any, 2> ann{std::string{"i,a"}, std::string{"a,i"}};

  // full transform: -conj(R^T)
  auto got = res->apply_transform(
      CanonTransform{.phase = -1, .conj = true, .braket_swap = true}, ann);
  ZArray ref;
  ref("a,i") = -1.0 * R("i,a").conj();
  world.gop.fence();
  ZArray diff;
  diff("a,i") = got->get<ZArray>()("a,i") - ref("a,i");
  REQUIRE(TA::norm2(diff) < 1e-12);

  // trivial transform: identical ResultPtr, no copy
  std::array<std::any, 2> ann_id{std::string{"i,a"}, std::string{"i,a"}};
  auto same = res->apply_transform(CanonTransform{}, ann_id);
  REQUIRE(same.get() == res.get());
}
```
  (Adjust the result-class and norm helper names to the file's existing
  usage — `ResultTensorTA` and the TA idioms above are the ones
  `test_eval_ta.cpp` already exercises; the assertions are the contract.)

- [ ] **Step 2: FAIL** (no `apply_transform`).

- [ ] **Step 3: Implement:**
  - base `Result`:
```cpp
  [[nodiscard]] virtual ResultPtr apply_transform(
      CanonTransform t, std::array<std::any, 2> const& ann) const {
    // default: compose existing pieces (adjoint retired in Task 9 makes the
    // tensor backends override this with a fused expression)
    throw detail::unimplemented_method("apply_transform");
  }
```
  - `ResultScalar<T>`: ignore annotations; `phase * (t.conj ? conj(value) : value)`
    (with `conj` a no-op for real `T` via `if constexpr`).
  - TA (both flat and ToT): one fused expression, conj elided for real
    numeric_type exactly as today's `adjoint` does, `braket_swap=false` uses
    `ann[1] == ann[0]`:
```cpp
    result(post) = t.phase == 1
        ? maybe_conj(get<ArrayT>()(pre))
        : t.phase * maybe_conj(get<ArrayT>()(pre));
```
  - BTAS/tapp: same composition using their existing `adjoint`/`mult_by_phase`
    building blocks.

- [ ] **Step 4: PASS.  Step 5: Commit** — `eval: Result::apply_transform (fused phase/conj/swap on retrieval)`

---

### Task 9: `evaluate` retrieves through the transform; Adjoint retired

**Files:**
- Modify: `SeQuant/core/eval/eval.hpp` (`apply_phase` :528 → `apply_canon_transform`; delete `NeedLeftAdj` :653-690; aliasing checks :674, :743-746, :817, :878, :1400)
- Modify: `SeQuant/core/eval/eval_expr.hpp` (delete `EvalOp::Adjoint`, `is_adjoint`)
- Modify: `SeQuant/core/export/export.hpp` (:224 Adjoint case — emit the swap-only expression exactly as today's documented real-field limitation, now keyed off a leaf transform rather than a node)
- Test: existing suites are the net.

- [ ] **Step 1:** `apply_phase` → `apply_canon_transform`: trivial transform
  returns `res` unchanged; else build the annotation pair
  `{node canonical annot, node denoted annot}` (denoted = canonical with
  bra/ket groups swapped when `braket_swap`; reuse the annot machinery the
  Adjoint dispatcher used) and call `res->apply_transform(...)`. Trace block
  keeps `MultByPhase`-style logging with the transform pretty-printed.

- [ ] **Step 2:** Delete the `NeedLeftAdj` stage + `EvalOp::Adjoint` + every
  `is_adjoint` consumer; update the five `canon_phase() != 1` aliasing checks
  to `!canon_transform().trivial()`.

- [ ] **Step 3:** Full suite (`'~[long-tests]'`) — including the PR-1
  eval-cache identity regressions (`eval_expr_conjugation_marker_identity`,
  `conjugate eval fold`, `ta_tot_adjoint_end_to_end` — the latter re-blessed
  to the transform model in Task 6). Long tests too: run once with no filter.

- [ ] **Step 4: Commit** — `eval: retrieval through apply_canon_transform; EvalOp::Adjoint retired`

---

### Task 10: conjugated `Variable` / `Power` leaves

**Files:**
- Modify: `SeQuant/core/eval/eval_expr.cpp` (`EvalExpr(Variable)` :196, `EvalExpr(Power)` :202)
- Test: `tests/unit/test_eval_expr.cpp`

- [ ] **Step 1: Failing test:**

```cpp
TEST_CASE("conjugated_scalar_leaves", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  Variable x{L"x"};
  Variable xs = x;
  xs.conjugate();
  REQUIRE(EvalExpr{x}.hash_value() == EvalExpr{xs}.hash_value());
  REQUIRE(EvalExpr{xs}.canon_transform().conj);
  REQUIRE(EvalExpr{x}.canon_transform().trivial());
}
```
  (Mirror for a conjugated `Power` whose base is a Variable.)

- [ ] **Step 2: FAIL.  Step 3:** in both ctors: if marked, hash the unmarked
  clone and set `canon_transform_.conj = true` (`conj(b^n) = conj(b)^n`).
  `ResultScalar::apply_transform` from Task 8 already serves it.

- [ ] **Step 4: PASS.  Step 5: Commit** — `eval: conjugated Variable/Power leaves served via the transform`

---

### Task 11: end-to-end conj evaluation + cache-reuse battery

**Files:**
- Test: `tests/unit/test_eval_ta.cpp`

- [ ] **Step 1:** Numeric tests on random COMPLEX data (fixture idiom of
  `eval_with_tiledarray`), each asserting against a hand-built TA reference.
  The reuse assertions all use one instrumented leaf yielder — build it once:

```cpp
// counts how many times evaluation asks for each leaf label; a cache hit
// leaves the count unchanged
struct CountingYielder {
  std::map<std::wstring, int> counts;
  std::map<std::wstring, ResultPtr> data;  // pre-built complex TA arrays
  ResultPtr operator()(EvalExprNode const& leaf) {
    auto lbl = std::wstring(leaf->as_tensor().label());
    ++counts[lbl];
    return data.at(lbl);
  }
};
```

  Tests (each its own SECTION, each with a TA-expression reference):
  1. per-symmetry serving: Hermitian `g` (Conjugate) evaluated from either
     orientation spelling == reference; starred Nonsymm `t^*` == `conj(t)`;
     '⁺' == conj-transpose; '⁺'+marker == transpose — e.g.

```cpp
    auto node = binarize(parse_expr(L"t^*{i_1;a_1}:N-N-S"));
    auto got = evaluate(node, yielder, cache);
    ZArray ref;
    ref("i,a") = data_t("i,a").conj();
    // compare got vs ref as in result_apply_transform_ta
```

  2. **Kramers hoisting reuse:** evaluate `A·B`; record
     `counts_after_first = yielder.counts`; evaluate `conjugate(A·B)` with
     the SAME cache; `REQUIRE(yielder.counts == counts_after_first)` (pure
     cache hit) and result == `conj(A·B)` reference;
  3. **mixed-term CSE:** with `A·B` cached, evaluate `A^*·B^*·C`;
     `REQUIRE(yielder.counts[L"A"] == 1)` and same for B (intermediate hit;
     only C newly yielded); result matches the no-reuse reference computed
     on a fresh cache;
  4. uniform-conj SUM of products (the \mathcal{T}-partner shape):
     `conjugate(A·B + C·D)` after `A·B + C·D` — zero new yields + equality.

- [ ] **Step 2:** All green; run everything including `[long-tests]` once.

- [ ] **Step 3: Commit** — `eval: conj evaluation + hoisting cache-reuse battery`

---

### Task 12: wrap-up

- [ ] `clang-format --style=file -i` on every touched file; full suite; commit any fallout.
- [ ] Update `doc/dev/specs/2026-09-01-lazy-conj-eval-design.md` if any
  design detail shifted during implementation (record the delta, don't
  rewrite history).
- [ ] Push `kshitij/feature/conjugation-eval`; open the DRAFT PR against
  `kshitij/feature/conjugation-symbolic` titled
  "Lazy-conj eval: retrieval-time canonical transform (draft for evaluation)"
  with the spec inlined in the description and the golden-churn justification.
  **Only on explicit user instruction to push/open.**

---

## Extended in-draft scope (decided 2026-09-01): Tasks 13-16

Owner decision: the exact `fold_conjugate_pairs` ({s, s*} -> 2 Re(s),
{s, -s*} -> 2i Im(s)) is the wanted fold; `fold_conjugate_pairs_of_real_sum`
is fragile (unverifiable reality assertion; difference pairs silently
unfolded) and is to be retired. Since that needs eval-side Re/Im ingestion,
the following Plan-B items move INTO this draft PR.

### Task 13: Re/Im eval nodes (IR + backends)

**Files:** `SeQuant/core/eval/eval_expr.hpp/.cpp`, `SeQuant/core/eval/eval.hpp`,
`SeQuant/core/eval/result.hpp` + `backends/{tiledarray,btas,tapp}/result.hpp`,
tests `test_eval_expr.cpp`, `test_eval_ta.cpp`.

- Failing IR test first: `binarize` of the fold emission
  `Constant(2) * RealPart(s)` (s = closed-contraction scalar network) yields a
  Product node whose RealPart factor is a UNARY node -- `EvalOp::RealPart`,
  `ResultType::Scalar`, `Constant{1}` sentinel right child (the retired
  Adjoint's FullBinaryNode pattern), hash = inner subtree hash combined with
  the EvalOp -- over the SHARED inner subtree (inner slots identical to the
  ones other consumers of s's pieces use). Same for ImagPart.
- `EvalOp::RealPart`/`ImagPart` enumerators; dispatcher case in
  `binarize(ExprPtr,...)` (replacing the "unsupported expression" throw for
  these types); recursive `binarize(inner)` under the wrapper.
- evaluate: unary compute branch `result = f.left->real_part(...)` /
  `imag_part(...)` (right sentinel ignored), cache/store machinery unchanged.
- `Result::real_part()/imag_part()` virtuals, default-throwing like
  `apply_transform`; overrides: `ResultScalar` (std::real/imag),
  TA (elementwise), BTAS, tapp.
- Numeric tests: `Re(A) + i*Im(A) == A` on random complex data; inner-slot
  sharing proven with the counting yielder; both suites green
  (`build-tests-debug` THROW + `build-ta-tests`).

### Task 14: Flip `FoldConjugatePairs` default to Yes

`options.hpp` member default `No -> Yes` (the parked "until the evaluation
layer understands them" condition is now met -- update that comment); fix
suite-wide fallout; verify `2*Re(A)` from the auto-fold evaluates equal to
`A + A^*` computed with the fold off.

### Task 15: Deprecate `fold_conjugate_pairs_of_real_sum`

`[[deprecated]]` with rationale in `expr_algorithms.hpp`; migrate or
warning-silence its direct unit tests; removal deferred until no caller
remains (MPQC switches in Task 16).

### Task 16: MPQC CC path to the exact fold + smoke re-certification

On the MPQC smoke setup (`smoke/pr2-eval` SeQuant +
`kshitij/refactor/kramer-pairs-separation`): switch
`src/mpqc/chemistry/qc/lcao/cc/sequant.cpp` from
`fold_conjugate_pairs_of_real_sum` to `fold_conjugate_pairs` (same
`swap_spin` conjugate_op), rebuild, re-run the h2o/dch decks against the
recorded baselines; capture the eval-dump proof of Re-node evaluation with
shared inner slots; resolve the open dch ~3e-9 old-vs-new-binary delta;
revert the temporary diagnostics (eval.hpp cache trace, eval_expr.cpp leaf
trace, MPQC fold-skip gate) before the draft goes up.

## Plan B (queued, separate document)

Remaining after the in-draft pull: export conj emission via `wrap_conj`
(tensor conj in generated code), and the exhaustive
`test_eval_conjugation.cpp` symbolic-to-eval matrix mirroring every PR-1
`test_conjugation.cpp` case. Written after this PR lands.

---

## Execution deviations (recorded 2026-09-01, tasks 13-15)

- The default flip (T14) needed hardening well beyond the one-line change:
  - simplify() gates the fold on `expr->is_cnumber()` -- operator-carrying
    intermediates head back into Wick, which does not ingest Re/Im wrappers
    (unfolded, a worker-thread TN ctor threw and terminated the process).
  - Wrapped summands from DIFFERENT simplify passes hold conjugate-related
    inners (for a closed c-number network the adjoint IS the conjugate);
    the fold now merges Re/Im-wrapped summands at entry AND exit over
    canonical representatives (Re(x*) = Re(x), Im(x*) = -Im(x)), so
    e.g. +c Re(X) - c Re(X^+) cancels exactly.
  - RealPart/ImagPart canonicalize their inner in place (+ REAL-scalar
    hoist via the byproduct contract); Product::canonicalize_impl resets
    its memoized hash when a subfactor mutates in place (first mutating
    non-tensor factor canonicalization ever -- the stale-hash self-check
    fired).
  - has_tensor sees through the wrappers (top level and product factors).
  - Folding before the canonicalize stage was tried and REVERTED: it
    preempts the Symm-braket collapse (a real-field a - a^T pair became
    2i Im instead of 0).
  - Diagnosis was repeatedly misled by (a) a file-local has_tensor lambda
    in test_mbpt_cc shadowing the core function, and (b) wide fwprintf
    probes silently no-opping on a byte-oriented stderr. Test fallout:
    fold-shape assertions compare canonically; UCC energy term counts
    halve per folded pair (46->23, 20->14, 74->41); re_im_evaluation pins
    canonicalizer + context (order-independent).
- T15 landed as staged: [[deprecated]] + pragma-wrapped self-test.

## Execution deviations (recorded 2026-09-01, tasks 4/6/7)

- Tasks 4, 6, and the product half of 7 landed as ONE green unit: the leaf
  ctor's marker-clearing changes binarize's observable behavior, so the
  plan's "IR-shape tests still pass" prediction at T4 was wrong.
- The plan's T4 step-2 "unfold" (marker => {conj, braket_swap} + slot ops)
  was WRONG: it collapses a starred-canonical spelling's transform onto the
  plain spelling and re-aliases C with C^*. Implemented rule: the marker
  composes a PURE {conj} bit (slots untouched); orientation deltas come from
  the step-3 fold alone.
- collect_tensor_factors re-materializes each leaf's DENOTED spelling
  (swap + marker) for TN building and slot counting -- the canonical
  respelling otherwise moves indices across bra/ket slots (strict-braket
  assert) and erases the marker coloring that keeps mixed-conj products
  identity-distinct.
- Uniform-conj hoisting is implemented at product roots (marker stripping on
  the collected TN copies + {conj} on the node); sum-level hoisting and the
  cache-reuse tests remain in T7/T11.

## Execution deviations (recorded 2026-09-02, T16 MPQC smoke)

- The MPQC CC-path fold (process_equations) needed THREE upstream repairs
  before pairing worked on the Kramers-CSV energy:
  (1) expand + flatten ALL CSV flavor sums before the korbit rebase
  (residuals included; 48 flat terms/block replace 12 nested);
  (2) rebase the fully contracted energy too, so a member's TRS partner is
  its plain elementwise conjugate -- conjugate_op = sequant::conjugate on
  the CSV path, mbpt::swap_spin on the non-CSV path (unrebased members
  pair by all-flipped relabeling);
  (3) round2-side kr_flavor guard: a spin-bit-free IndexSpace (the DF aux)
  was misread as flavored under IGNOREd asserts (to_spin(qns) on zero spin
  bits) and the rebase minted spurious flavored aux dummies, defeating
  every pairing op. Measured after repairs: dch 24/36 paired -> 12 Re
  (eq0 36 -> 24 terms), h2o 8/10.
- optimize_impl was OPAQUE to Re/Im wrappers (returned untouched -> naive
  inner contraction order). Fixed with a see-through case + regression
  test (test_optimize.cpp "optimize sees through Re/Im wrappers").
- OPEN: a Re-wrapped ToT/CSV scalar summand's inner root evaluates through
  a MATERIALIZING DeNest einsum instead of the plain summand's scalar
  reduction (binarize_re_im's inner binarize lacks the ResultExpr root
  treatment): ~14 GB peak / OOM on dch vs 1.7 GB unfolded. Flat-TA wrapper
  eval is certified (h2o). The MPQC CSV-energy fold therefore ships opt-in
  (MPQC_CCK_TRS_FOLD=1) until the wrapper's inner root is routed through
  the same result-expression treatment as an unwrapped scalar term.
- Certified after all changes: h2o CCk -0.13510773505222942 (fold on),
  dch PNS-MP2 -1.04169026886 (default, fold auto-off; within the known
  X2C thread scatter of the old reference).

## Follow-up roadmap (2026-09-02)

- **T18 -- Kramers-canonical configs only (symbolic layer): RESOLVED on the
  Kramers round-2 branch** by a time-reversal-aware canonicalizer
  (`KramersSymmetry` tensor attribute, registry partner spaces, network
  fold under `CanonicalizeOptions::fold_kramers`; design and plan under
  `doc/dev/{specs,plans}/2026-09-02-kramers-trs-canonicalizer-*.md`).
  Measured on dch: energy conjugate-pair fold 36 -> 20 terms (floor);
  residual blocks are externally anchored and unchanged.
- **T19 -- MPQC serving-level aliasing.** Serve only the up-row C arrays;
  down-row and conj-orientation leaves become CanonTransform views of the
  same slot. Use TA's lazy `.conj()` inside the consuming expression rather
  than materializing the transformed copy on retrieval; this also enables
  the eval-leaf Kramers fold (`fold_kramers_eval_leaves`).
- **T20 -- wrapped-summand CSV eval fix.** Route binarize_re_im's inner root
  through the ResultExpr scalar-head treatment so a Re-wrapped ToT summand
  reduces instead of materializing (14 GB defect); then flip the CSV
  energy fold to default-on (MPQC_CCK_TRS_FOLD gate removed).
