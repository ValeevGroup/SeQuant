# Extended Wick Theorem Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evaluate products of generalized-normal-ordered (GNO) fermion operators relative to a multiconfigurational reference (`Vacuum::MultiProduct`), producing GNO strings times one-body densities γ, hole densities η and density cumulants κₖ.

**Architecture:** The existing `WickTheorem` runs unchanged except for (a) `Vacuum::MultiProduct` quasiparticle classifiers, (b) a `MultiProduct` branch in `contract` that emits γ/η instead of the overlap `s`, and (c) skipping its pair-based connectivity filters under `MultiProduct`. A post-processing pass (`cumulant_expand`) turns each partial-contraction term's surviving operators into cumulant blocks; `extended_wick` wraps canonicalization + Wick + pass. `mbpt::ref_av` dispatches to it under `MultiProduct`.

**Tech Stack:** C++20, Catch2 (`unit_tests-sequant`), CMake/Ninja. Build tree: `cmake-build-debug` (Debug, `SEQUANT_ASSERT_BEHAVIOR=THROW`, unity ON — see Global Constraints).

**Spec:** `docs/superpowers/specs/2026-10-01-extended-wick-design.md`

## Global Constraints

- Throw only `sequant::Exception` (or subclasses); use `SEQUANT_ASSERT` for internal invariants (AGENTS.md).
- One logical change per commit; commit messages in the repo style `area: lowercase summary`, no attribution trailers of any kind (AGENTS.md overrides the harness's Co-Authored-By instruction).
- Smallest diff that achieves the change; no drive-by refactoring.
- Never create, modify or regenerate `*.expected` fixtures.
- Format every touched C++ file with `bin/admin/clang-format.sh -i <files>` before committing.
- `cmake-build-debug` has `CMAKE_UNITY_BUILD=ON`; after adding/moving includes compile the TU standalone: `cmake -DCMAKE_UNITY_BUILD=OFF cmake-build-debug && ninja -C cmake-build-debug "$PWD/<file>.cpp^"`, then restore with `cmake -DCMAKE_UNITY_BUILD=ON cmake-build-debug`.
- Build: `cmake --build cmake-build-debug --target unit_tests-sequant -j8`. Run: `cmake-build-debug/tests/unit/unit_tests-sequant "<filter>"`.
- Spin-orbital only: `SPBasis::Spinfree` with `Vacuum::MultiProduct` must keep throwing (engine already does).
- Tensor slot convention for γ/η/κ: **bra = annihilator indices, ket = creator indices**, matching `ref_av`'s existing γ (`mbpt/op.cpp:1334-1365`) and `NormalOperator` printing (`ã{ann;cre}`).
- Labels: `γ` (1-RDM), `η` (1-hole RDM), `κ` (cumulants, rank ≥ 2). All Hermitian + column-symmetric; κ additionally `Symmetry::Antisymm`.
- Work on a branch off `master` (e.g. `feature/extended-wick`), not on the current conjugation branch.

## Review Focus

Failure modes the spec implies but which no "happy path" test exercises; each has a pinning test in the task that owns the code:

1. **Vacuum mismatch** — an operator built under `SingleProduct` handed to `extended_wick` in a `MultiProduct` context must throw `sequant::Exception`, not silently produce SR output (Task 6).
2. **Unbalanced active survivors** — `⟨{a†_u a†_u}{a_u}⟩` has no balanced block; full mode must return `0`, not a rank-1 or lopsided κ (Task 5).
3. **Protoindexed index reaching the active space** — must trip the `SEQUANT_ASSERT` in `contract` (throws under `THROW`), never emit a protoindexed γ (Task 4).
4. **`max_cumulant_rank` of 0 or 1** — means "no cumulant blocks": only γ/η terms survive; must not throw or form rank-1 blocks (Task 5).
5. **Bosons under `MultiProduct`** — `BWickTheorem` must still assert/throw rather than classify (Task 3).

---

### Task 1: Rename `cumu*_to_density` → `cumulant*_to_density`

**Files:**
- Modify: `SeQuant/domain/mbpt/rdm.hpp:15-19`
- Modify: `SeQuant/domain/mbpt/rdm.cpp` (definitions of the three functions and every internal caller)
- Modify: `tests/unit/test_mbpt.cpp:1769`

**Interfaces:**
- Produces: `ExprPtr sequant::mbpt::decompositions::cumulant_to_density(ExprPtr)`, `cumulant2_to_density(ExprPtr)`, `cumulant3_to_density(ExprPtr)` — same semantics as before (each takes a single κ `Tensor` of rank 1/2/3).

- [ ] **Step 1: Rename everywhere**

```bash
grep -rln 'cumu_to_density\|cumu2_to_density\|cumu3_to_density' SeQuant tests
# expected: SeQuant/domain/mbpt/rdm.hpp SeQuant/domain/mbpt/rdm.cpp tests/unit/test_mbpt.cpp
sed -i '' -e 's/cumu_to_density/cumulant_to_density/g' \
          -e 's/cumu2_to_density/cumulant2_to_density/g' \
          -e 's/cumu3_to_density/cumulant3_to_density/g' \
    SeQuant/domain/mbpt/rdm.hpp SeQuant/domain/mbpt/rdm.cpp tests/unit/test_mbpt.cpp
grep -rn 'cumu_to\|cumu2_to\|cumu3_to' SeQuant tests doc python  # expected: no output
```

- [ ] **Step 2: Build and run the affected test**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "[mbpt]" -c "rdm-decomposition symmetries"`
Expected: `All tests passed`

- [ ] **Step 3: Commit**

```bash
bin/admin/clang-format.sh -i SeQuant/domain/mbpt/rdm.hpp SeQuant/domain/mbpt/rdm.cpp tests/unit/test_mbpt.cpp
git add SeQuant/domain/mbpt/rdm.hpp SeQuant/domain/mbpt/rdm.cpp tests/unit/test_mbpt.cpp
git commit -m "mbpt: spell out cumulant in the cumulant-to-density helpers"
```

---

### Task 2: Shared density/cumulant tensor factory (`SeQuant/core/density.hpp`)

Factors the γ construction duplicated between `mbpt/rdm.cpp` (anonymous-namespace labels + `hermitian_particle_symmetric`) and `mbpt/op.cpp`'s `replace_nop_with_rdm` into one header both can use, and adds η. The Wick engine (Task 4) and the pass (Task 5) need these from `SeQuant-symb`, so the header lives under `SeQuant/core/`.

**Files:**
- Create: `SeQuant/core/density.hpp`
- Modify: `SeQuant/domain/mbpt/rdm.cpp:6-33` (drop the anonymous-namespace labels/symmetries, use `density::`)
- Modify: `SeQuant/domain/mbpt/op.cpp:1334-1365` (`replace` lambda → `density::rdm_from_nop`)
- Modify: `SeQuant/domain/mbpt/context.cpp:148-150` (register `η`)
- Modify: `CMakeLists.txt:435` (add `SeQuant/core/density.hpp` to the symb header list next to `wick.hpp`)
- Test: `tests/unit/test_mbpt.cpp` section `"rdm-decomposition symmetries"`

**Interfaces:**
- Produces (all in `namespace sequant::density`):
  ```cpp
  std::wstring rdm_label();        // L"γ"
  std::wstring hole_rdm_label();   // L"η"
  std::wstring cumulant_label();   // L"κ"
  constexpr TensorSymmetries rdm_symmetries{.hermiticity = Hermiticity::Hermitian, .column = ColumnSymmetry::Symm};
  constexpr TensorSymmetries cumulant_symmetries{.perm = Symmetry::Antisymm, .hermiticity = Hermiticity::Hermitian, .column = ColumnSymmetry::Symm};
  ExprPtr make_rdm(const Index& ann, const Index& cre);        // γ{ann;cre}
  ExprPtr make_hole_rdm(const Index& ann, const Index& cre);   // η{ann;cre}
  template <Statistics S> ExprPtr rdm_from_nop(const NormalOperator<S>& nop, std::wstring_view label, TensorSymmetries syms);
  template <Statistics S> ExprPtr make_cumulant(const NormalOperator<S>& nop); // κ, requires nop.rank() >= 2
  ```

- [ ] **Step 1: Write the failing test** (append inside `SECTION("rdm-decomposition symmetries")` in `tests/unit/test_mbpt.cpp`, after the existing `REQUIRE(*gamma == *gamma_op);`)

```cpp
  // the shared factory in SeQuant/core/density.hpp is the single source of
  // these spellings
  REQUIRE(*density::make_rdm(Index(L"i_1"), Index(L"i_2")) == *gamma_op);
  const auto eta = density::make_hole_rdm(Index(L"i_1"), Index(L"i_2"));
  REQUIRE(eta->as<Tensor>().label() == L"η");
  REQUIRE(eta->as<Tensor>().hermiticity() == Hermiticity::Hermitian);
  REQUIRE(eta->as<Tensor>().column_symmetry() == ColumnSymmetry::Symm);

  // κ from a 2-body normal operator: bra = annihilators, ket = creators
  const FNOperator nop2(cre({L"i_1", L"i_2"}), ann({L"i_3", L"i_4"}));
  const auto kappa2 = density::make_cumulant(nop2);
  REQUIRE(kappa2->as<Tensor>().label() == L"κ");
  REQUIRE(kappa2->as<Tensor>().symmetry() == Symmetry::Antisymm);
  REQUIRE(kappa2->as<Tensor>().bra()[0] == Index(L"i_3"));
  REQUIRE(kappa2->as<Tensor>().bra()[1] == Index(L"i_4"));
  REQUIRE(kappa2->as<Tensor>().ket()[0] == Index(L"i_1"));
  REQUIRE(kappa2->as<Tensor>().ket()[1] == Index(L"i_2"));
  // η is a registered mbpt operator label
  REQUIRE(mbpt::get_default_mbpt_context().op_registry()->contains(L"η"));
```

Add `#include <SeQuant/core/density.hpp>` to the test's include block.

- [ ] **Step 2: Run to verify it fails**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 2>&1 | grep -m3 error`
Expected: compile error, `density` is not a namespace.

- [ ] **Step 3: Create the header**

```cpp
// SeQuant/core/density.hpp
#ifndef SEQUANT_CORE_DENSITY_HPP
#define SEQUANT_CORE_DENSITY_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>

#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include <string>
#include <string_view>

namespace sequant::density {

/// labels of the reference-state density tensors produced by the extended
/// Wick theorem and consumed by mbpt::decompositions
inline std::wstring rdm_label() { return L"γ"; }
inline std::wstring hole_rdm_label() { return L"η"; }
inline std::wstring cumulant_label() { return L"κ"; }

/// densities and cumulants describe indistinguishable particles, hence are
/// column symmetric, and are Hermitian by definition; both take part in the
/// tensor hash, so every producer must spell them identically
inline constexpr TensorSymmetries rdm_symmetries{
    .hermiticity = Hermiticity::Hermitian, .column = ColumnSymmetry::Symm};
inline constexpr TensorSymmetries cumulant_symmetries{
    .perm = Symmetry::Antisymm,
    .hermiticity = Hermiticity::Hermitian,
    .column = ColumnSymmetry::Symm};

/// @return γ with bra = @p ann and ket = @p cre, i.e. ⟨a†_cre a_ann⟩
inline ExprPtr make_rdm(const Index &ann, const Index &cre) {
  return ex<Tensor>(rdm_label(), bra{ann}, ket{cre}, rdm_symmetries);
}

/// @return η with bra = @p ann and ket = @p cre, i.e. ⟨a_ann a†_cre⟩ = δ - γ
inline ExprPtr make_hole_rdm(const Index &ann, const Index &cre) {
  return ex<Tensor>(hole_rdm_label(), bra{ann}, ket{cre}, rdm_symmetries);
}

/// @return the reference expectation value of @p nop as a tensor with the
/// given label: bra = annihilator indices, ket = creator indices, both in
/// particle order
template <Statistics S>
ExprPtr rdm_from_nop(const NormalOperator<S> &nop, std::wstring_view label,
                     TensorSymmetries syms) {
  using index_container = container::svector<Index>;
  auto braidxs = nop.annihilators() |
                 ranges::views::transform(
                     [](const auto &op) { return op.index(); }) |
                 ranges::to<index_container>();
  auto ketidxs = nop.creators() |
                 ranges::views::transform(
                     [](const auto &op) { return op.index(); }) |
                 ranges::to<index_container>();
  SEQUANT_ASSERT(braidxs.size() == ketidxs.size());
  return ex<Tensor>(std::wstring(label), bra(std::move(braidxs)),
                    ket(std::move(ketidxs)), syms);
}

/// @return κ_k for the k-body @p nop (k ≥ 2)
template <Statistics S>
ExprPtr make_cumulant(const NormalOperator<S> &nop) {
  SEQUANT_ASSERT(nop.rank() >= 2);
  return rdm_from_nop(nop, cumulant_label(), cumulant_symmetries);
}

}  // namespace sequant::density

#endif  // SEQUANT_CORE_DENSITY_HPP
```

Add `SeQuant/core/density.hpp` to the header list in `CMakeLists.txt` (the block around line 435 that lists `SeQuant/core/wick.hpp`).

- [ ] **Step 4: Use it in `rdm.cpp`**

Delete the anonymous-namespace `rdm_label()`, `rdm_cumulant_label()`, `hermitian_particle_symmetric` (keep `particle_symmetric`, which is a different pack). Add `#include <SeQuant/core/density.hpp>` and at the top of `namespace sequant::mbpt::decompositions` add

```cpp
using density::rdm_label;
using density::rdm_symmetries;
```

then replace every `rdm_cumulant_label()` with `density::cumulant_label()` and every `hermitian_particle_symmetric` with `rdm_symmetries`.

- [ ] **Step 5: Use it in `op.cpp`**

In `expectation_value_impl`, replace the body of the `replace` lambda (`mbpt/op.cpp:1334-1358`) with

```cpp
      auto replace = [&rdm_label, spinor](const auto& nop) -> ExprPtr {
        // spin-free RDMs are column symmetric but not antisymmetric
        const auto syms =
            nop.rank() > 1 && spinor
                ? TensorSymmetries{.perm = Symmetry::Antisymm,
                                   .hermiticity = Hermiticity::Hermitian,
                                   .column = ColumnSymmetry::Symm}
                : TensorSymmetries{.perm = Symmetry::Nonsymm,
                                   .hermiticity = Hermiticity::Hermitian,
                                   .column = ColumnSymmetry::Symm};
        return density::rdm_from_nop(nop, rdm_label, syms);
      };
```

and add `#include <SeQuant/core/density.hpp>`. (`rdm_label` stays a local `const wchar_t*` because the spin-free path uses `Γ`.)

- [ ] **Step 6: Register `η`** in `mbpt/context.cpp` after `.add(L"γ", OpClass::Gen)`:

```cpp
      .add(L"η", OpClass::Gen)
```

- [ ] **Step 7: Build, run tests, check non-unity compile of the touched TUs**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "[mbpt]"`
Expected: `All tests passed` (the MRSO `γ` expectation at `test_mbpt.cpp:1352` still matches — same label/symmetries).

Run: `cmake -DCMAKE_UNITY_BUILD=OFF cmake-build-debug && ninja -C cmake-build-debug "$PWD/SeQuant/domain/mbpt/rdm.cpp^" "$PWD/SeQuant/domain/mbpt/op.cpp^" && cmake -DCMAKE_UNITY_BUILD=ON cmake-build-debug`
Expected: both compile.

- [ ] **Step 8: Commit**

```bash
bin/admin/clang-format.sh -i SeQuant/core/density.hpp SeQuant/domain/mbpt/rdm.cpp SeQuant/domain/mbpt/op.cpp SeQuant/domain/mbpt/context.cpp tests/unit/test_mbpt.cpp
git add SeQuant/core/density.hpp CMakeLists.txt SeQuant/domain/mbpt/rdm.cpp SeQuant/domain/mbpt/op.cpp SeQuant/domain/mbpt/context.cpp tests/unit/test_mbpt.cpp
git commit -m "core: one factory for the reference density tensors γ, η and κ"
```

---

### Task 3: `Vacuum::MultiProduct` quasiparticle classifiers

**Files:**
- Modify: `SeQuant/core/op.hpp:200-332` (the six `MultiProduct` cases)
- Test: `tests/unit/test_wick.cpp` section `"Op contractions"`

**Interfaces:**
- Produces: `is_pure_qpcreator`, `is_qpcreator`, `qpcreator_space`, `is_pure_qpannihilator`, `is_qpannihilator`, `qpannihilator_space`, `can_contract` all accept `Vacuum::MultiProduct`. R = `isr->reference_occupied_space(qns)`, U = `isr->vacuum_unoccupied_space(qns)`.

- [ ] **Step 1: Write the failing test** (append to `SECTION("Op contractions")` in `test_wick.cpp`, before its closing brace)

```cpp
    // MultiProduct vacuum: i = core (vacuum-occupied), u = active
    // (reference-occupied, vacuum-unoccupied), a = virtual
    {
      auto ctx = get_default_context();
      ctx.set(mbpt::make_mr_spaces());
      ctx.set(Vacuum::MultiProduct);
      auto ctx_resetter = set_scoped_default_context(ctx);
      const auto isr = ctx.index_space_registry();
      const auto MP = Vacuum::MultiProduct;

      // core behaves as in SingleProduct
      REQUIRE(is_pure_qpannihilator(fcre(L"i_1"), MP, isr));
      REQUIRE(!is_qpcreator(fcre(L"i_1"), MP, isr));
      REQUIRE(is_pure_qpcreator(fann(L"i_1"), MP, isr));
      REQUIRE(!is_qpannihilator(fann(L"i_1"), MP, isr));
      // virtual behaves as in SingleProduct
      REQUIRE(is_pure_qpcreator(fcre(L"a_1"), MP, isr));
      REQUIRE(!is_qpannihilator(fcre(L"a_1"), MP, isr));
      REQUIRE(is_pure_qpannihilator(fann(L"a_1"), MP, isr));
      REQUIRE(!is_qpcreator(fann(L"a_1"), MP, isr));
      // active is both
      REQUIRE(is_pure_qpcreator(fcre(L"u_1"), MP, isr));
      REQUIRE(is_pure_qpannihilator(fcre(L"u_1"), MP, isr));
      REQUIRE(is_pure_qpcreator(fann(L"u_1"), MP, isr));
      REQUIRE(is_pure_qpannihilator(fann(L"u_1"), MP, isr));
      // general p: both, but not pure
      REQUIRE(is_qpcreator(fcre(L"p_1"), MP, isr));
      REQUIRE(!is_pure_qpcreator(fcre(L"p_1"), MP, isr));
      REQUIRE(qpcreator_space(fcre(L"p_1"), MP, isr) ==
              isr->vacuum_unoccupied_space(IndexSpace::QuantumNumbers{}));
      REQUIRE(qpannihilator_space(fcre(L"p_1"), MP, isr) ==
              isr->reference_occupied_space(IndexSpace::QuantumNumbers{}));

      // contractions: only cre·ann over R and ann·cre over U
      REQUIRE(FWickTheorem::can_contract(fcre(L"i_1"), fann(L"i_2"), MP));
      REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fcre(L"i_2"), MP));
      REQUIRE(FWickTheorem::can_contract(fann(L"a_1"), fcre(L"a_2"), MP));
      REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fann(L"a_2"), MP));
      REQUIRE(FWickTheorem::can_contract(fcre(L"u_1"), fann(L"u_2"), MP));
      REQUIRE(FWickTheorem::can_contract(fann(L"u_1"), fcre(L"u_2"), MP));
      REQUIRE(!FWickTheorem::can_contract(fcre(L"u_1"), fcre(L"u_2"), MP));
      REQUIRE(!FWickTheorem::can_contract(fann(L"u_1"), fann(L"u_2"), MP));
      REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fann(L"a_2"), MP));
      REQUIRE(FWickTheorem::can_contract(fcre(L"u_1"), fann(L"i_2"), MP));
      REQUIRE(FWickTheorem::can_contract(fann(L"u_1"), fcre(L"a_2"), MP));

      // bosons stay Physical-only (Review Focus 5)
      REQUIRE_THROWS_AS(
          BWickTheorem::can_contract(bann(L"i_1"), bcre(L"i_2"), MP),
          Exception);
    }
```

Note: the mr registry's quantum numbers are `Spin::any`; if the `IndexSpace::QuantumNumbers{}` comparisons fail for that reason, use `fcre(L"p_1").index().space().qns()` as the argument instead.

- [ ] **Step 2: Run to verify it fails**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "[wick]" -c "Op contractions"`
Expected: FAIL — `Exception: is_pure_qpannihilator: cannot handle MultiProduct vacuum`.

- [ ] **Step 3: Implement the classifiers** in `SeQuant/core/op.hpp`. Replace each `case Vacuum::MultiProduct: throw ...` as follows (`sp = op.index().space()`, `qns = sp.qns()`):

```cpp
// is_pure_qpcreator
    case Vacuum::MultiProduct: {
      const auto &sp = op.index().space();
      const auto &target = op.action() == Action::Create
                               ? isr->vacuum_unoccupied_space(sp.qns())
                               : isr->reference_occupied_space(sp.qns());
      return isr->intersection(sp, target) == sp;
    }
// is_qpcreator  (N.B. fix the misplaced case label: it currently sits inside
// the SingleProduct block's braces; move it after them)
    case Vacuum::MultiProduct: {
      const auto &sp = op.index().space();
      const auto &target = op.action() == Action::Create
                               ? isr->vacuum_unoccupied_space(sp.qns())
                               : isr->reference_occupied_space(sp.qns());
      return static_cast<bool>(isr->intersection(sp, target));
    }
// qpcreator_space
    case Vacuum::MultiProduct: {
      const auto &sp = op.index().space();
      return op.action() == Action::Create
                 ? isr->intersection(sp, isr->vacuum_unoccupied_space(sp.qns()))
                 : isr->intersection(sp, isr->reference_occupied_space(sp.qns()));
    }
// is_pure_qpannihilator
    case Vacuum::MultiProduct: {
      const auto &sp = op.index().space();
      const auto &target = op.action() == Action::Annihilate
                               ? isr->vacuum_unoccupied_space(sp.qns())
                               : isr->reference_occupied_space(sp.qns());
      return isr->intersection(sp, target) == sp;
    }
// is_qpannihilator
    case Vacuum::MultiProduct: {
      const auto &sp = op.index().space();
      const auto &target = op.action() == Action::Annihilate
                               ? isr->vacuum_unoccupied_space(sp.qns())
                               : isr->reference_occupied_space(sp.qns());
      return static_cast<bool>(isr->intersection(sp, target));
    }
// qpannihilator_space
    case Vacuum::MultiProduct: {
      const auto &sp = op.index().space();
      return op.action() == Action::Annihilate
                 ? isr->intersection(sp, isr->vacuum_unoccupied_space(sp.qns()))
                 : isr->intersection(sp, isr->reference_occupied_space(sp.qns()));
    }
```

`IndexSpaceRegistry::intersection` returns `IndexSpace::null` (falsy) for disjoint spaces; check how the `SingleProduct` branches test for null and mirror that if `static_cast<bool>` is not the idiom (`IndexSpace` has `explicit operator bool` — `grep -n "operator bool" SeQuant/core/space.hpp`).

For bosons: in `WickTheorem<S>::can_contract` (`wick.hpp:1615`) the existing `SEQUANT_ASSERT(vacuum == Vacuum::Physical)` already covers Review Focus 5 under `THROW`; no change.

- [ ] **Step 4: Run the test**

Run: `cmake-build-debug/tests/unit/unit_tests-sequant "[wick]" -c "Op contractions"`
Expected: PASS.

- [ ] **Step 5: Audit exclusivity assumptions.** `grep -n "is_qpannihilator\|is_qpcreator\|is_pure_qp" SeQuant/core/wick.hpp SeQuant/core/wick.impl.hpp SeQuant/core/op.hpp SeQuant/core/*.cpp`. For each hit outside the classifiers themselves, note whether it assumes an op cannot be both. Expected hits: `recursive_nontensor_wick`'s early return (`wick.hpp:1227`, safe — an active op *is* a qp annihilator), `contract` (`wick.hpp:1639-1646`, handled in Task 4), `can_contract` (`op.hpp:342`). If any other site branches on "creator else annihilator", record it in the commit message body and handle it in Task 4.

- [ ] **Step 6: Run the whole Wick/op suite and commit**

Run: `cmake-build-debug/tests/unit/unit_tests-sequant "[wick],[op]"`
Expected: PASS.

```bash
bin/admin/clang-format.sh -i SeQuant/core/op.hpp tests/unit/test_wick.cpp
git add SeQuant/core/op.hpp tests/unit/test_wick.cpp
git commit -m "op: quasiparticle classification relative to a MultiProduct vacuum"
```

---

### Task 4: `contract` emits γ/η under `MultiProduct`; engine skips connectivity filters

**Files:**
- Modify: `SeQuant/core/wick.hpp:1626-1690` (`contract`), `:685-720` (`compute_nopseq`), add a getter next to `set_external_indices` (`:172`)
- Test: `tests/unit/test_wick.cpp`, new `SECTION("multiproduct vacuum")` after `SECTION("fermi vacuum")`

**Interfaces:**
- Produces: `WickTheorem<S>::contract(left, right, Vacuum::MultiProduct, isr)` returns `density::make_rdm`/`make_hole_rdm` (optionally wrapped with projection Kroneckers as today); `const std::optional<container::set<Index>>& WickTheorem<S>::external_indices() const`.
- Hook (spin-free seam #1): a private static `contraction_value(const Index& bra, const Index& ket, Action left_action, Vacuum)` chooses the middle factor; spin-free will later branch here.

- [ ] **Step 1: Write the failing tests**

```cpp
  SECTION("multiproduct vacuum") {
    auto ctx = get_default_context();
    ctx.set(mbpt::make_mr_spaces());
    ctx.set(Vacuum::MultiProduct);
    auto ctx_resetter = set_scoped_default_context(ctx);

    // raw engine output: two active 1-body operators, partial contractions
    // = {a†_u1 a_u2}{a†_u3 a_u4}: 4 terms, with active pairs emitted as γ
    // (cre·ann) and η (ann·cre) rather than overlaps
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                                     FNOperator(cre({L"u_3"}), ann({L"u_4"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(false).compute();
      REQUIRE_THAT(result,
                   EquivalentTo(L"ã{u_2,u_4;u_1,u_3} "
                                L"+ γ{u_4;u_1} * ã{u_2;u_3} "
                                L"- η{u_2;u_3} * ã{u_4;u_1} "
                                L"+ γ{u_4;u_1} * η{u_2;u_3}"));
    }

    // general indices: the γ-type contraction is over R (core+active) and the
    // η-type over U (active+virtual); projections are spelled with δ as in
    // the SingleProduct case
    {
      auto opseq = ex<FNOperatorSeq>(FNOperator(cre({L"p_1"}), ann({L"p_2"})),
                                     FNOperator(cre({L"p_3"}), ann({L"p_4"})));
      auto wick = FWickTheorem{opseq};
      auto result = wick.compute();
      REQUIRE_THAT(result, EquivalentTo(L"δ{p_4;M_1} * γ{M_1;M_2} * δ{M_2;p_1} * "
                                        L"δ{p_2;E_1} * η{E_1;E_2} * δ{E_2;p_3}"));
    }

    // a protoindexed index may never reach the active space (Review Focus 3)
    {
      const Index u1(L"u_1");
      const Index a_u1(L"a_1", {u1});  // a_1 depends on u_1, lives in virtual
      // fine: virtual-only contraction stays an overlap
      REQUIRE_NOTHROW(FWickTheorem::can_contract(fann(a_u1), fcre(L"a_2")));
      // not fine: a general index with protoindices projected onto R∩U
      const Index p_u1(L"p_1", {u1});
      REQUIRE_THROWS_AS(
          FWickTheorem::contract(Op<Statistics::FermiDirac>(p_u1, Action::Create),
                                 fann(L"u_2")),
          Exception);
    }
  }
```

The exact `M_1`/`E_1` spellings: `M` is `make_mr_spaces`' reference-occupied union and `E` its vacuum-unoccupied one; if the registry prints the latter differently, read the actual output once and fix the string — but only the *labels*, the structure must be as written.

- [ ] **Step 2: Run to verify it fails**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "[wick]" -c "multiproduct vacuum"`
Expected: FAIL — overlaps `s{...}` where γ/η are expected.

- [ ] **Step 3: Implement the middle-factor hook in `contract`.** In `wick.hpp`, before `contract`, add

```cpp
  /// the value of a single contraction between @p bra (the annihilator's
  /// index) and @p ket (the creator's index), both already projected onto
  /// the common quasiparticle space; this is where a spin-free variant
  /// would substitute its own tensors
  static ExprPtr contraction_value(const Index &bra, const Index &ket,
                                   bool left_is_annihilator, Vacuum vacuum) {
    if (vacuum != Vacuum::MultiProduct) return make_overlap(bra, ket);
    return left_is_annihilator ? density::make_hole_rdm(bra, ket)
                               : density::make_rdm(bra, ket);
  }
```

and in `contract` replace both `make_overlap(...)` calls with `contraction_value(..., left_is_ann, vacuum)` (same arguments otherwise). Add `#include <SeQuant/core/density.hpp>`.

Add the protoindex guard right after `qpspace_common` is computed:

```cpp
    if constexpr (S == Statistics::FermiDirac) {
      if (vacuum == Vacuum::MultiProduct &&
          (left.index().has_proto_indices() ||
           right.index().has_proto_indices())) {
        const auto &sp = left_is_pure && right_is_pure
                             ? isr->intersection(left.index().space(),
                                                 right.index().space())
                             : qpspace_common;
        const auto qns = sp.qns();
        const auto &active = isr->intersection(
            isr->reference_occupied_space(qns), isr->vacuum_unoccupied_space(qns));
        SEQUANT_ASSERT(!isr->intersection(sp, active) &&
                       "protoindexed indices must not reach the active space");
      }
    }
```

- [ ] **Step 4: Skip engine connectivity under `MultiProduct`.** In `compute_nopseq` (`wick.hpp:699-706`) guard the two cached-connection blocks:

```cpp
    // under a MultiProduct vacuum connectivity is a property of the
    // cumulant-expanded result (cumulant blocks connect operators that no
    // pair does), so the pair-based filters are applied by cumulant_expand
    const bool pairwise_connectivity =
        get_default_context(S).vacuum() != Vacuum::MultiProduct;
    if (pairwise_connectivity && !nop_connections_input_.empty()) ...
    if (pairwise_connectivity && !nop_avoided_connections_input_.empty()) ...
```

Also add the getter after `set_external_indices`:

```cpp
  /// @return the external indices, if known
  const std::optional<container::set<Index>> &external_indices() const {
    return external_indices_;
  }
```

(Confirm the member's exact type with `grep -n "external_indices_;" SeQuant/core/wick.hpp` and match it.)

- [ ] **Step 5: Run the tests**

Run: `cmake-build-debug/tests/unit/unit_tests-sequant "[wick]"`
Expected: PASS, including all pre-existing `SingleProduct`/`Physical` sections (their output must be byte-identical — `contraction_value` returns `make_overlap` for them).

- [ ] **Step 6: Non-unity compile check and commit**

Run: `cmake -DCMAKE_UNITY_BUILD=OFF cmake-build-debug && ninja -C cmake-build-debug "$PWD/SeQuant/core/wick.cpp^" && cmake -DCMAKE_UNITY_BUILD=ON cmake-build-debug`

```bash
bin/admin/clang-format.sh -i SeQuant/core/wick.hpp tests/unit/test_wick.cpp
git add SeQuant/core/wick.hpp tests/unit/test_wick.cpp
git commit -m "wick: a MultiProduct contraction is a density, not an overlap"
```

---

### Task 5: The pass — `cumulant_expand` on pure-active survivors

The core combinatorics, tested on inputs whose survivors are already pure-active so no splitting/projection is needed. Mixed spaces, provenance recovery and the wrapper are Task 6.

**Files:**
- Create: `SeQuant/core/wick_extended.hpp`, `SeQuant/core/wick_extended.cpp`
- Create: `tests/unit/test_wick_extended.cpp`
- Modify: `CMakeLists.txt:435` (add both sources next to `wick.cpp`), `tests/unit/CMakeLists.txt:28` (add `"test_wick_extended.cpp"` after `"test_wick.cpp"`)

**Interfaces:**
- Produces:
  ```cpp
  namespace sequant {
  struct ExtendedWickOptions {
    bool full_contractions = true;
    std::optional<std::size_t> max_cumulant_rank;  // nullopt = unbounded; 0 or 1 = no cumulant blocks
    bool eta_as_delta_minus_gamma = false;         // Task 8
    bool use_topology = true;                      // Task 6
    container::svector<std::pair<std::size_t, std::size_t>> nop_connections;          // Task 7
    container::svector<std::pair<std::size_t, std::size_t>> nop_avoided_connections;  // Task 7
  };
  /// provenance: which input NormalOperator (0-based ordinal) each surviving Op's index came from
  using OpProvenance = container::map<Index, std::size_t>;
  template <Statistics S>
  ExprPtr cumulant_expand(const ExprPtr& wick_output, const OpProvenance& provenance, const ExtendedWickOptions& opts);
  }
  ```
- Hook (spin-free seam #2): `detail::block_value<S>(const NormalOperator<S>& block)` returns `density::make_cumulant(block)`; and `detail::term_weight<S>(...)` returns `1` — both named functions, so a spin-free variant has one place to override.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_wick_extended.cpp`)

```cpp
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/core/wick_extended.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

TEST_CASE("wick_extended", "[algorithms][wick][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests;

  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  ctx.set(Vacuum::MultiProduct);
  auto ctx_resetter = set_scoped_default_context(ctx);

  // helper: standard Wick with all partial contractions + provenance from the
  // (uncanonicalized) input; adequate here because every index is external
  auto wick_partial = [](const FNOperatorSeq& nopseq, OpProvenance& prov) {
    prov.clear();
    std::size_t ord = 0;
    for (const auto& nop : nopseq) {
      for (const auto& op : nop) prov.emplace(op.index(), ord);
      ++ord;
    }
    FWickTheorem wick{std::make_shared<FNOperatorSeq>(nopseq)};
    return wick.full_contractions(false).compute();
  };

  SECTION("cumulant_expand: ⟨{a†a†a}{a}⟩ = κ2") {
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({L"u_3"})),
                     FNOperator(cre({}), ann({L"u_4"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE_THAT(result, EquivalentTo(L"κ{u_4,u_3;u_1,u_2}:A-H-S"));
  }

  SECTION("cumulant_expand: ⟨{a†a}{a†a}⟩ = γη + κ2") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1} * η{u_2;u_3} "
                                      L"+ κ{u_2,u_4;u_1,u_3}:A-H-S"));
  }

  SECTION("cumulant_expand: max_cumulant_rank = 1 means pairs only") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 1});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1} * η{u_2;u_3}"));
    auto result0 = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 0});
    REQUIRE(simplify(result - result0) == ex<Constant>(0));
  }

  SECTION("cumulant_expand: unbalanced survivors vanish") {
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({})),
                     FNOperator(cre({}), ann({L"u_3"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE(result == ex<Constant>(0));
  }

  SECTION("cumulant_expand: a block within one operator vanishes") {
    // {a†_u1 a†_u2 a_u3 a_u4}{a†_u5 a_u6}: the 4 legs of nop 0 may not form
    // a block by themselves; every κ must involve u_5 or u_6
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"})),
                     FNOperator(cre({L"u_5"}), ann({L"u_6"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE(result->is<Sum>());
    for (const auto& term : *result) {
      bool has_kappa = false, kappa_touches_nop1 = false;
      for (const auto& f : *term) {
        if (f->is<Tensor>() && f->as<Tensor>().label() == L"κ") {
          has_kappa = true;
          for (const auto& idx : f->as<Tensor>().braket())
            if (idx == Index(L"u_5") || idx == Index(L"u_6"))
              kappa_touches_nop1 = true;
        }
      }
      if (has_kappa) REQUIRE(kappa_touches_nop1);
    }
  }

  SECTION("cumulant_expand: 3-body truncation") {
    // {a†a†a}{a†a a} has a κ3 term (all six legs) that max_cumulant_rank=2
    // must drop, leaving everything else unchanged
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({L"u_3"})),
                     FNOperator(cre({L"u_4"}), ann({L"u_5", L"u_6"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto full = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    auto trunc = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 2});
    auto diff = simplify(full - trunc);
    // exactly the κ3 term
    REQUIRE_THAT(diff, EquivalentTo(L"κ{u_5,u_6,u_3;u_1,u_2,u_4}:A-H-S"));
  }

  SECTION("cumulant_expand: partial contractions leave a GNO remainder") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.full_contractions = false});
    // Eq. (ext. Wick, 1-body×1-body): the Wick output's 4 terms plus κ2;
    // nothing else, because a block needs ≥2 ops from ≥2 nops and the only
    // such balanced set is all four legs
    REQUIRE_THAT(result, EquivalentTo(L"ã{u_2,u_4;u_1,u_3} "
                                      L"+ γ{u_4;u_1} * ã{u_2;u_3} "
                                      L"- η{u_2;u_3} * ã{u_4;u_1} "
                                      L"+ γ{u_4;u_1} * η{u_2;u_3} "
                                      L"+ κ{u_2,u_4;u_1,u_3}:A-H-S"));
    for (const auto& term : *result)
      for (const auto& f : *term)
        if (f->is<FNOperator>())
          REQUIRE(f->as<FNOperator>().vacuum() == Vacuum::MultiProduct);
  }
}
```

If the deserializer's symmetry-suffix spelling for an antisymmetric Hermitian column-symmetric tensor differs from `:A-H-S`, check `tests/unit/test_parse.cpp` for the current grammar and adjust the suffix only.

The κ3 sign in the truncation test and the signs in the partial test are what the implementation below yields for a correctly implemented permutation parity; if a sign disagrees, derive it by hand (write the six operators in storage order, count transpositions to bring the block legs to the front) before touching the code.

- [ ] **Step 2: Register the sources and run to verify failure**

Add `SeQuant/core/wick_extended.hpp` / `SeQuant/core/wick_extended.cpp` to the symb list in `CMakeLists.txt` and `"test_wick_extended.cpp"` to `tests/unit/CMakeLists.txt`. Create an empty header/cpp pair with just include guards and `namespace sequant {}` so CMake configures.

Run: `cmake cmake-build-debug && cmake --build cmake-build-debug --target unit_tests-sequant -j8 2>&1 | grep -m3 error`
Expected: `ExtendedWickOptions`/`cumulant_expand` undeclared.

- [ ] **Step 3: Implement `wick_extended.hpp`**

```cpp
#ifndef SEQUANT_CORE_WICK_EXTENDED_HPP
#define SEQUANT_CORE_WICK_EXTENDED_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/density.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <cstddef>
#include <optional>
#include <utility>

namespace sequant {

/// controls extended_wick() and cumulant_expand()
struct ExtendedWickOptions {
  /// keep only terms with no surviving operators
  bool full_contractions = true;
  /// largest cumulant rank to form; nullopt = no bound; 0 or 1 = no cumulants
  std::optional<std::size_t> max_cumulant_rank;
  /// rewrite every η as δ - γ
  bool eta_as_delta_minus_gamma = false;
  /// forwarded to WickTheorem::use_topology
  bool use_topology = true;
  /// pairs of input NormalOperator ordinals that must end up connected
  container::svector<std::pair<std::size_t, std::size_t>> nop_connections;
  /// pairs of input NormalOperator ordinals that must not be connected
  container::svector<std::pair<std::size_t, std::size_t>>
      nop_avoided_connections;
};

/// maps the index of a surviving Op to the ordinal of the input
/// NormalOperator it came from
using OpProvenance = container::map<Index, std::size_t>;

namespace detail {

/// the value of a cumulant block; a spin-free variant would override this
template <Statistics S>
ExprPtr block_value(const NormalOperator<S> &block) {
  return density::make_cumulant(block);
}

/// a per-term scalar weight (1 for spin-orbital); a spin-free variant would
/// put its cycle factor here
template <Statistics S>
rational term_weight(const NormalOperator<S> & /*survivors*/,
                     const container::svector<container::svector<std::size_t>>
                         & /*blocks*/) {
  return 1;
}

/// enumerates every way to group the ops of @p survivors into disjoint
/// cumulant blocks and (unless @p full) a remainder, and reports each via
/// @p sink(sign, blocks, remainder) where blocks are op positions in storage
/// order, remainder is the list of op positions not in any block, and sign is
/// the parity of moving each block's legs (in order) to the front
/// @note a block has k creators and k annihilators, 2 <= k <= max_rank, and
///       legs from at least two distinct provenance ordinals
template <Statistics S, typename Sink>
void for_each_block_assignment(const NormalOperator<S> &survivors,
                               const OpProvenance &provenance,
                               std::size_t max_rank, bool full, Sink &&sink);

}  // namespace detail

/// expands the surviving operators of each term of a MultiProduct-vacuum
/// WickTheorem output (partial contractions) into cumulant blocks
/// @pre every γ/η/surviving index is pure-active (see extended_wick for the
///      general case)
template <Statistics S>
ExprPtr cumulant_expand(const ExprPtr &wick_output,
                        const OpProvenance &provenance,
                        const ExtendedWickOptions &opts);

extern template ExprPtr cumulant_expand<Statistics::FermiDirac>(
    const ExprPtr &, const OpProvenance &, const ExtendedWickOptions &);

}  // namespace sequant

#endif  // SEQUANT_CORE_WICK_EXTENDED_HPP
```

- [ ] **Step 4: Implement `wick_extended.cpp`**

```cpp
#include <SeQuant/core/wick_extended.hpp>

#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>

#include <range/v3/algorithm/count_if.hpp>

#include <algorithm>
#include <functional>

namespace sequant {

namespace detail {

template <Statistics S, typename Sink>
void for_each_block_assignment(const NormalOperator<S> &survivors,
                               const OpProvenance &provenance,
                               std::size_t max_rank, bool full, Sink &&sink) {
  const std::size_t n = survivors.size();
  // assignment[i] = block id, or npos for "remainder"
  constexpr std::size_t npos = static_cast<std::size_t>(-1);
  container::svector<std::size_t> assignment(n, npos);
  container::svector<container::svector<std::size_t>> blocks;

  auto is_cre = [&](std::size_t i) {
    return survivors[i].action() == Action::Create;
  };
  auto prov = [&](std::size_t i) {
    auto it = provenance.find(survivors[i].index());
    if (it == provenance.end())
      throw Exception("cumulant_expand: surviving operator index " +
                      to_string(survivors[i].index().full_label()) +
                      " has no provenance");
    return it->second;
  };

  auto emit = [&]() {
    // parity of the permutation [block0 legs][block1 legs]...[remainder]
    // relative to storage order: count inversions
    container::svector<std::size_t> perm;
    perm.reserve(n);
    for (const auto &b : blocks) perm.insert(perm.end(), b.begin(), b.end());
    container::svector<std::size_t> remainder;
    for (std::size_t i = 0; i != n; ++i)
      if (assignment[i] == npos) remainder.push_back(i);
    if (full && !remainder.empty()) return;
    perm.insert(perm.end(), remainder.begin(), remainder.end());
    int sign = 1;
    if constexpr (S == Statistics::FermiDirac) {
      for (std::size_t i = 0; i != n; ++i)
        for (std::size_t j = i + 1; j != n; ++j)
          if (perm[i] > perm[j]) sign = -sign;
    }
    sink(sign, blocks, remainder);
  };

  // recursive: decide the fate of the first undecided op
  std::function<void(std::size_t)> recurse = [&](std::size_t first) {
    while (first != n && assignment[first] != npos) ++first;
    if (first == n) {
      emit();
      return;
    }
    // option A: remainder (only if partial contractions wanted)
    if (!full) {
      recurse(first + 1);
    }
    // option B: start a new block whose lowest leg is `first`
    if (max_rank < 2) return;
    const std::size_t block_id = blocks.size();
    const bool first_is_cre = is_cre(first);
    // candidates after `first` that are still undecided
    container::svector<std::size_t> cre_cands, ann_cands;
    for (std::size_t i = first + 1; i != n; ++i)
      if (assignment[i] == npos)
        (is_cre(i) ? cre_cands : ann_cands).push_back(i);
    // choose k-1 more of first's kind and k of the other kind
    auto &same = first_is_cre ? cre_cands : ann_cands;
    auto &other = first_is_cre ? ann_cands : cre_cands;
    for (std::size_t k = 2; k <= max_rank; ++k) {
      if (same.size() + 1 < k || other.size() < k) break;
      // iterate over k-1 subsets of `same` and k subsets of `other` via
      // selection bitmasks (std::prev_permutation over a 0/1 vector)
      container::svector<bool> sel_same(same.size(), false);
      std::fill(sel_same.begin(), sel_same.begin() + (k - 1), true);
      do {
        container::svector<bool> sel_other(other.size(), false);
        std::fill(sel_other.begin(), sel_other.begin() + k, true);
        do {
          container::svector<std::size_t> block{first};
          for (std::size_t i = 0; i != same.size(); ++i)
            if (sel_same[i]) block.push_back(same[i]);
          for (std::size_t i = 0; i != other.size(); ++i)
            if (sel_other[i]) block.push_back(other[i]);
          std::sort(block.begin(), block.end());
          // a block within a single GNO string vanishes
          const auto p0 = prov(block[0]);
          const bool spans_two = std::any_of(
              block.begin(), block.end(),
              [&](std::size_t i) { return prov(i) != p0; });
          if (spans_two) {
            for (auto i : block) assignment[i] = block_id;
            blocks.push_back(block);
            recurse(first + 1);
            blocks.pop_back();
            for (auto i : block) assignment[i] = npos;
          }
        } while (std::prev_permutation(sel_other.begin(), sel_other.end()));
      } while (std::prev_permutation(sel_same.begin(), sel_same.end()));
    }
  };
  recurse(0);
}

}  // namespace detail

template <Statistics S>
ExprPtr cumulant_expand(const ExprPtr &wick_output,
                        const OpProvenance &provenance,
                        const ExtendedWickOptions &opts) {
  const std::size_t max_rank =
      opts.max_cumulant_rank.value_or(static_cast<std::size_t>(-1));
  auto result = std::make_shared<Sum>();

  auto expand_term = [&](const ExprPtr &term) {
    // term is a Product (possibly with a single NormalOperator<S> factor)
    // or a bare NormalOperator<S>
    ProductPtr product = term->is<Product>()
                             ? term->as<Product>().clone().as_shared_ptr<Product>()
                             : std::make_shared<Product>(ExprPtrList{term->clone()});
    // locate the survivors
    auto nop_it = std::find_if(product->begin(), product->end(),
                               [](const ExprPtr &f) {
                                 return f->template is<NormalOperator<S>>();
                               });
    if (nop_it == product->end()) {  // fully contracted already
      result->append(product);
      return;
    }
    const auto survivors = (*nop_it)->template as<NormalOperator<S>>();
    auto prefactor = std::make_shared<Product>(product->scalar(), ExprPtrList{});
    for (const auto &f : *product)
      if (f != *nop_it) prefactor->append(1, f);

    detail::for_each_block_assignment<S>(
        survivors, provenance, max_rank, opts.full_contractions,
        [&](int sign,
            const container::svector<container::svector<std::size_t>> &blocks,
            const container::svector<std::size_t> &remainder) {
          auto summand = prefactor->clone().template as_shared_ptr<Product>();
          summand->scale(sign * detail::term_weight<S>(survivors, blocks));
          for (const auto &b : blocks) {
            container::svector<Op<S>> cre_ops, ann_ops;
            for (auto i : b)
              (survivors[i].action() == Action::Create ? cre_ops : ann_ops)
                  .push_back(survivors[i]);
            // NormalOperator takes annihilators in particle order, i.e.
            // reversed relative to storage order
            std::reverse(ann_ops.begin(), ann_ops.end());
            NormalOperator<S> block(cre(cre_ops), ann(ann_ops),
                                    Vacuum::MultiProduct);
            summand->append(1, detail::block_value<S>(block));
          }
          if (!remainder.empty()) {
            container::svector<Op<S>> cre_ops, ann_ops;
            for (auto i : remainder)
              (survivors[i].action() == Action::Create ? cre_ops : ann_ops)
                  .push_back(survivors[i]);
            std::reverse(ann_ops.begin(), ann_ops.end());
            summand->append(1, ex<NormalOperator<S>>(cre(cre_ops), ann(ann_ops),
                                                     Vacuum::MultiProduct));
          }
          result->append(summand);
        });
  };

  if (wick_output->is<Sum>()) {
    for (const auto &term : *wick_output) expand_term(term);
  } else if (wick_output->is<Constant>()) {
    return wick_output;
  } else {
    expand_term(wick_output);
  }
  ExprPtr out = result;
  simplify(out);
  if (out->is<Sum>() && out->as<Sum>().empty()) return ex<Constant>(0);
  return out;
}

template ExprPtr cumulant_expand<Statistics::FermiDirac>(
    const ExprPtr &, const OpProvenance &, const ExtendedWickOptions &);

}  // namespace sequant
```

Adjust to the real APIs as you go (check with grep, do not guess): `Product::scale` signature (`grep -n "scale(" SeQuant/core/expressions/product.hpp`), `Product::begin/end` over factors, `ExprPtr::operator!=`, `rational` include (`SeQuant/core/rational.hpp`), `cre(...)`/`ann(...)` accepting `svector<Op<S>>` (they are templates over any range of Op — see `op.hpp:538-563`). `std::function` is acceptable here; if the include police object, convert `recurse` to a private struct with a method.

- [ ] **Step 5: Build and run**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "wick_extended"`
Expected: all sections PASS. If a sign disagrees with the test string, resolve by hand as described in Step 1 before changing either.

- [ ] **Step 6: Non-unity compile check and commit**

Run: `cmake -DCMAKE_UNITY_BUILD=OFF cmake-build-debug && ninja -C cmake-build-debug "$PWD/SeQuant/core/wick_extended.cpp^" "$PWD/tests/unit/test_wick_extended.cpp^" && cmake -DCMAKE_UNITY_BUILD=ON cmake-build-debug`

```bash
bin/admin/clang-format.sh -i SeQuant/core/wick_extended.hpp SeQuant/core/wick_extended.cpp tests/unit/test_wick_extended.cpp
git add SeQuant/core/wick_extended.hpp SeQuant/core/wick_extended.cpp tests/unit/test_wick_extended.cpp CMakeLists.txt tests/unit/CMakeLists.txt
git commit -m "wick: cumulant_expand groups surviving operators into cumulant blocks"
```

---

### Task 6: `extended_wick` wrapper — canonicalization, provenance, mixed-space splitting

**Files:**
- Modify: `SeQuant/core/wick_extended.hpp`, `SeQuant/core/wick_extended.cpp`
- Test: `tests/unit/test_wick_extended.cpp`

**Interfaces:**
- Produces:
  ```cpp
  template <Statistics S>
  ExprPtr extended_wick(ExprPtr input, const ExtendedWickOptions& opts = {});
  // input: Product/Sum of Products containing NormalOperator<S> factors, or a NormalOperatorSequence<S>
  ```
  Throws `Exception` unless `get_default_context(S).vacuum() == Vacuum::MultiProduct`.
- Consumes: `cumulant_expand<S>`, `WickTheorem<S>::external_indices()` (Task 4), `density::make_rdm/make_hole_rdm`.

- [ ] **Step 1: Write the failing tests** (append sections to `test_wick_extended.cpp`)

```cpp
  SECTION("extended_wick: vacuum must be MultiProduct") {
    auto sr_ctx = get_default_context();
    sr_ctx.set(Vacuum::SingleProduct);
    auto sr_resetter = set_scoped_default_context(sr_ctx);
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    REQUIRE_THROWS_AS(extended_wick<Statistics::FermiDirac>(in), Exception);
  }

  SECTION("extended_wick: pure-active identities via the wrapper") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = extended_wick<Statistics::FermiDirac>(in);
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1} * η{u_2;u_3} "
                                      L"+ κ{u_2,u_4;u_1,u_3}:A-H-S"));
  }

  SECTION("extended_wick: general indices split into core δ + active γ") {
    // ⟨{a†_p1 a_p2}{a†_p3 a_p4}⟩ with p over core(i) + active(u) + virtual(a)
    // (o and g exist in make_mr_spaces too, but p = M ∪ E covers them):
    // cre·ann pair over R = δ on core + γ on active;
    // ann·cre pair over U = δ on virtual + η on active;
    // plus κ2 on the all-active projection
    auto in = ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"})) *
              ex<FNOperator>(cre({L"p_3"}), ann({L"p_4"}));
    auto result = extended_wick<Statistics::FermiDirac>(in);
    // every γ/η/κ index is active; every other index is a δ to an external p
    for (const auto& term : *result)
      for (const auto& f : *term)
        if (f->is<Tensor>()) {
          const auto& t = f->as<Tensor>();
          if (t.label() == L"γ" || t.label() == L"η" || t.label() == L"κ")
            for (const auto& idx : t.braket())
              REQUIRE(idx.space() == Index(L"u_1").space());
        }
    // and the SR-looking piece is present: δ(p4,O)δ(O,p1) δ(p2,E')δ(E',p3)
    // spelled with whatever tmp labels reduce produced — check by count:
    // 4 terms = δδ, δγ... no: (δ_core + γ)(δ_virt + η) + κ = 5 terms
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 5);
  }

  SECTION("extended_wick: topology on/off agree") {
    auto in = ex<FNOperator>(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"})) *
              ex<FNOperator>(cre({L"u_5", L"u_6"}), ann({L"u_7", L"u_8"}));
    auto with = extended_wick<Statistics::FermiDirac>(in, {.use_topology = true});
    auto without =
        extended_wick<Statistics::FermiDirac>(in, {.use_topology = false});
    REQUIRE(simplify(with - without) == ex<Constant>(0));
  }

  SECTION("extended_wick: Sum input") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
                  ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) +
              ex<Constant>(2) * ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
                  ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = extended_wick<Statistics::FermiDirac>(in);
    REQUIRE_THAT(result, EquivalentTo(L"3 γ{u_4;u_1} * η{u_2;u_3} "
                                      L"+ 3 κ{u_2,u_4;u_1,u_3}:A-H-S"));
  }
```

- [ ] **Step 2: Run to verify failure**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 2>&1 | grep -m3 error`
Expected: `extended_wick` undeclared.

- [ ] **Step 3: Declare `extended_wick`** in `wick_extended.hpp` (after `cumulant_expand`):

```cpp
/// applies the extended (generalized-normal-order) Wick theorem to @p input
/// @param input a Product or Sum of Products with NormalOperator<S> factors
///        normal-ordered relative to Vacuum::MultiProduct, or an
///        ExprPtr to a NormalOperatorSequence<S>
/// @throw Exception if the default context's vacuum is not MultiProduct
template <Statistics S>
ExprPtr extended_wick(ExprPtr input, const ExtendedWickOptions &opts = {});

extern template ExprPtr extended_wick<Statistics::FermiDirac>(
    ExprPtr, const ExtendedWickOptions &);
```

- [ ] **Step 4: Implement** in `wick_extended.cpp`. The implementation has three parts; write them as three file-local functions in `namespace sequant::detail`.

(a) **Provenance from a canonicalized Product/nopseq:**

```cpp
template <Statistics S>
OpProvenance make_provenance(const Expr &expr) {
  OpProvenance prov;
  std::size_t ord = 0;
  auto record = [&](const NormalOperator<S> &nop) {
    for (const auto &op : nop) prov.emplace(op.index(), ord);
    ++ord;
  };
  if (expr.is<NormalOperatorSequence<S>>()) {
    for (const auto &nop : expr.as<NormalOperatorSequence<S>>()) record(nop);
  } else if (expr.is<Product>()) {
    for (const auto &f : expr.as<Product>())
      if (f->template is<NormalOperator<S>>())
        record(f->template as<NormalOperator<S>>());
  } else if (expr.is<NormalOperator<S>>()) {
    record(expr.as<NormalOperator<S>>());
  }
  return prov;
}
```

(b) **Mixed-space splitting of γ/η and of survivors.** For a single Product term: for every γ/η factor and every survivor op whose index space is not pure core / pure active / pure virtual, substitute the index by a sum over its components, using the same `δ(idx, tmp)` idiom as `contract`. Implement as: collect replacement sums, then expand.

```cpp
/// splits a γ/η over a mixed space into its δ (core resp. virtual) and
/// active parts, and projects surviving ops onto their active (and, for
/// partial contractions, inactive) components; returns an expanded Sum of
/// Products in which every γ/η/survivor index is pure
template <Statistics S>
ExprPtr split_mixed_spaces(const ExprPtr &term, const Context &ctx,
                           bool full_contractions) {
  const auto &isr = ctx.index_space_registry();
  auto product = term->is<Product>()
                     ? term->as<Product>().clone().as_shared_ptr<Product>()
                     : std::make_shared<Product>(ExprPtrList{term->clone()});
  auto qns_of = [](const Index &i) { return i.space().qns(); };
  auto core_of = [&](const Index &i) {
    return isr->vacuum_occupied_space(qns_of(i));
  };
  auto active_of = [&](const Index &i) {
    return isr->intersection(isr->reference_occupied_space(qns_of(i)),
                             isr->vacuum_unoccupied_space(qns_of(i)));
  };
  auto virt_of = [&](const Index &i) {
    // U minus active = complement of R
    return isr->intersection(isr->vacuum_unoccupied_space(qns_of(i)),
                             isr->complement(isr->reference_occupied_space(qns_of(i))));
  };
  auto is_pure = [&](const Index &i) {
    const auto &sp = i.space();
    return sp == core_of(i) || sp == active_of(i) || sp == virt_of(i) ||
           isr->intersection(sp, active_of(i)) == IndexSpace::null;
  };

  // factor-wise rewrite: each factor becomes a Sum of alternatives
  ExprPtr rewritten = ex<Constant>(product->scalar());
  for (const auto &f : *product) {
    ExprPtr alternative;
    if (f->template is<Tensor>()) {
      const auto &t = f->template as<Tensor>();
      const bool is_gamma = t.label() == density::rdm_label();
      const bool is_eta = t.label() == density::hole_rdm_label();
      if ((is_gamma || is_eta) && t.rank() == 1) {
        const Index &b = t.bra()[0], &k = t.ket()[0];
        if (!is_pure(b) || !is_pure(k)) {
          // δ part: both indices projected onto core (γ) or virtual (η)
          const auto &dsp = is_gamma ? isr->intersection(core_of(b), core_of(k))
                                     : isr->intersection(virt_of(b), virt_of(k));
          auto sum = std::make_shared<Sum>();
          if (dsp) {
            const auto d = Index::make_tmp_index(dsp);
            sum->append(make_kronecker(b, d) * make_kronecker(d, k));
          }
          const auto &asp = isr->intersection(active_of(b), active_of(k));
          if (asp) {
            const auto ub = Index::make_tmp_index(asp);
            const auto uk = Index::make_tmp_index(asp);
            auto g = is_gamma ? density::make_rdm(ub, uk)
                              : density::make_hole_rdm(ub, uk);
            sum->append(make_kronecker(b, ub) * g * make_kronecker(uk, k));
          }
          alternative = sum;
        }
      }
    } else if (f->template is<NormalOperator<S>>()) {
      const auto &nop = f->template as<NormalOperator<S>>();
      // project each op whose space is not pure onto its components; in full
      // mode only the active component can survive (the rest is dropped by
      // for_each_block_assignment's balance check, but pruning here is cheap)
      auto sum = std::make_shared<Sum>();
      container::svector<container::svector<std::pair<Op<S>, ExprPtr>>> choices;
      for (const auto &op : nop) {
        container::svector<std::pair<Op<S>, ExprPtr>> alts;
        const Index &i = op.index();
        if (is_pure(i)) {
          alts.emplace_back(op, nullptr);
        } else {
          for (const auto &comp : {active_of(i), core_of(i), virt_of(i)}) {
            const auto &sp = isr->intersection(i.space(), comp);
            if (!sp) continue;
            if (full_contractions && sp != active_of(i)) continue;
            const auto j = Index::make_tmp_index(sp);
            alts.emplace_back(Op<S>(j, op.action()),
                              op.action() == Action::Create
                                  ? make_kronecker(j, i)   // ket side
                                  : make_kronecker(i, j));  // bra side
          }
        }
        choices.push_back(std::move(alts));
      }
      // cartesian product over choices
      std::function<void(std::size_t, container::svector<Op<S>>, ExprPtr)> go =
          [&](std::size_t pos, container::svector<Op<S>> ops, ExprPtr deltas) {
            if (pos == choices.size()) {
              container::svector<Op<S>> cre_ops, ann_ops;
              for (const auto &o : ops)
                (o.action() == Action::Create ? cre_ops : ann_ops).push_back(o);
              std::reverse(ann_ops.begin(), ann_ops.end());
              ExprPtr piece = ex<NormalOperator<S>>(cre(cre_ops), ann(ann_ops),
                                                    Vacuum::MultiProduct);
              sum->append(deltas ? deltas * piece : piece);
              return;
            }
            for (const auto &[o, d] : choices[pos]) {
              auto ops2 = ops;
              ops2.push_back(o);
              go(pos + 1, std::move(ops2),
                 d ? (deltas ? deltas * d : d) : deltas);
            }
          };
      go(0, {}, nullptr);
      alternative = sum;
    }
    rewritten = rewritten * (alternative ? alternative : f);
  }
  expand(rewritten);
  return rewritten;
}
```

If `IndexSpaceRegistry` has no `complement` (check: `grep -n "complement" SeQuant/core/index_space_registry.hpp`), compute `virt_of` as `intersection(vacuum_unoccupied, intersection(complete_space, ...))` via the registry's set algebra, or define it as "U with the active bits cleared" using `IndexSpace::Type` bit operations (`IndexSpace::Type` is a bitset-like; see `space.hpp`). Pick whichever exists; do not add a registry method for this.

(c) **The wrapper:**

```cpp
template <Statistics S>
ExprPtr extended_wick(ExprPtr input, const ExtendedWickOptions &opts) {
  const auto &ctx = get_default_context(S);
  if (ctx.vacuum() != Vacuum::MultiProduct)
    throw Exception(
        "extended_wick: the default context's vacuum must be "
        "Vacuum::MultiProduct");

  // one term at a time so that provenance is per term
  auto per_term = [&](ExprPtr term) -> ExprPtr {
    // canonicalize first so the provenance map sees the indices WickTheorem
    // will see; then tell it not to canonicalize again
    if (term->is<Product>()) {
      [[maybe_unused]] auto bp = term->rapid_canonicalize();
      SEQUANT_ASSERT(bp == nullptr);
    }
    const auto provenance = detail::make_provenance<S>(*term);

    WickTheorem<S> wick{term};
    wick.full_contractions(false).use_topology(opts.use_topology);
    auto raw = wick.compute(/*count_only=*/false,
                            /*skip_input_canonicalization=*/true);
    if (!raw || raw->is<Constant>()) return raw ? raw : ex<Constant>(0);

    // split mixed spaces, then reduce the δ chains this introduced
    auto split = std::make_shared<Sum>();
    auto split_one = [&](const ExprPtr &t) {
      auto s = detail::split_mixed_spaces<S>(t, ctx, opts.full_contractions);
      if (s->is<Sum>())
        for (const auto &x : *s) split->append(x);
      else
        split->append(s);
    };
    if (raw->is<Sum>())
      for (const auto &t : *raw) split_one(t);
    else
      split_one(raw);
    ExprPtr reduced = split;
    wick.reduce(reduced);

    // provenance of projected survivors: a tmp index j introduced above
    // replaced an original index i; record j with i's ordinal. reduce() may
    // have renamed j again via the δ chain, so rebuild from the result: a
    // survivor's index is either original (in provenance) or is bound by a δ
    // to an original index in the same term — resolve through that δ.
    auto resolved = [&](const Product &p, const Index &j) -> std::size_t {
      if (auto it = provenance.find(j); it != provenance.end()) return it->second;
      for (const auto &f : p)
        if (f->template is<Tensor>() &&
            f->template as<Tensor>().label() == reserved::kronecker_label()) {
          const auto &t = f->template as<Tensor>();
          const Index &b = t.bra()[0], &k = t.ket()[0];
          if (b == j) if (auto it = provenance.find(k); it != provenance.end()) return it->second;
          if (k == j) if (auto it = provenance.find(b); it != provenance.end()) return it->second;
        }
      throw Exception("extended_wick: cannot resolve provenance of a projected survivor");
    };
    auto expanded = std::make_shared<Sum>();
    auto expand_one = [&](const ExprPtr &t) {
      OpProvenance prov = provenance;
      if (t->is<Product>())
        for (const auto &f : t->template as<Product>())
          if (f->template is<NormalOperator<S>>())
            for (const auto &op : f->template as<NormalOperator<S>>())
              prov.emplace(op.index(), resolved(t->template as<Product>(), op.index()));
      expanded->append(cumulant_expand<S>(t, prov, opts));
    };
    if (reduced->is<Sum>())
      for (const auto &t : *reduced) expand_one(t);
    else
      expand_one(reduced);
    ExprPtr out = expanded;
    return simplify(out);
  };

  expand(input);
  ExprPtr result;
  if (input->is<Sum>()) {
    auto sum = std::make_shared<Sum>();
    for (const auto &t : *input) sum->append(per_term(t));
    result = sum;
  } else {
    result = per_term(input);
  }
  simplify(result);
  if (result->is<Sum>() && result->as<Sum>().empty()) return ex<Constant>(0);
  return result;
}

template ExprPtr extended_wick<Statistics::FermiDirac>(ExprPtr, const ExtendedWickOptions &);
```

Verification point from the spec ("nothing between input canonicalization and output renames a surviving index"): after implementing, run the "general indices" test with `Logger::instance().wick_harness = true` once and confirm each surviving `ã` index in the raw output is either in the provenance map or δ-bound to one; if `resolved` throws for any term, that assumption is wrong — stop and report rather than widening the resolver heuristically.

- [ ] **Step 5: Build and run**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "wick_extended"`
Expected: PASS. For the "general indices" section, if the term count is not 5, print `to_latex(result)` and reconcile against the hand derivation in the test comment before changing the expectation.

- [ ] **Step 6: Commit**

```bash
bin/admin/clang-format.sh -i SeQuant/core/wick_extended.hpp SeQuant/core/wick_extended.cpp tests/unit/test_wick_extended.cpp
git add SeQuant/core/wick_extended.hpp SeQuant/core/wick_extended.cpp tests/unit/test_wick_extended.cpp
git commit -m "wick: extended_wick applies the generalized-normal-order theorem"
```

---

### Task 7: Connectivity constraints in the pass

**Files:**
- Modify: `SeQuant/core/wick_extended.cpp` (filter inside `cumulant_expand`)
- Test: `tests/unit/test_wick_extended.cpp`

**Interfaces:**
- Consumes: `ExtendedWickOptions::nop_connections`, `nop_avoided_connections`.
- Produces: terms violating a required connection, or realizing an avoided one, are dropped. Connectivity graph per term: nop ordinals are vertices; a γ/η/δ pair contributes an edge between the ordinals of its two indices (resolved through provenance; indices renamed by `reduce` are resolved through δ as in Task 6); a κ contributes edges among all ordinals of its legs.

- [ ] **Step 1: Write the failing test**

```cpp
  SECTION("extended_wick: connectivity") {
    // {a†_u1 a_u2}{a†_u3 a_u4}{a†_u5 a_u6}: require 0-1 connected and
    // forbid 0-2. The κ3 term connects all three (kept); κ2(0,1)·pair(0,2)
    // is forbidden; pair(0,1)·pair(1,2)... enumerate by property instead of
    // by string: every kept term has an edge 0-1 and none 0-2.
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
              ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    auto all = extended_wick<Statistics::FermiDirac>(in);
    auto filtered = extended_wick<Statistics::FermiDirac>(
        in, {.nop_connections = {{0, 1}}, .nop_avoided_connections = {{0, 2}}});
    REQUIRE(filtered->size() < all->size());
    auto ord = [](const Index& i) -> int {
      const auto l = i.label();
      if (l == L"u_1" || l == L"u_2") return 0;
      if (l == L"u_3" || l == L"u_4") return 1;
      return 2;
    };
    for (const auto& term : *filtered) {
      bool e01 = false, e02 = false;
      for (const auto& f : *term) {
        if (!f->is<Tensor>()) continue;
        const auto& t = f->as<Tensor>();
        container::set<int> ords;
        for (const auto& idx : t.braket()) ords.insert(ord(idx));
        if (ords.contains(0) && ords.contains(1)) e01 = true;
        if (ords.contains(0) && ords.contains(2)) e02 = true;
      }
      REQUIRE(e01);
      REQUIRE(!e02);
    }
  }
```

- [ ] **Step 2: Run to verify failure**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "wick_extended" -c "extended_wick: connectivity"`
Expected: FAIL (`filtered->size() < all->size()` false, options ignored).

- [ ] **Step 3: Implement.** In `cumulant_expand`, inside the sink lambda after `summand` is fully built, before `result->append(summand)`:

```cpp
          if (!opts.nop_connections.empty() ||
              !opts.nop_avoided_connections.empty()) {
            // edges realized by this term
            container::set<std::pair<std::size_t, std::size_t>> edges;
            auto add_edges = [&](const container::svector<std::size_t> &ords) {
              for (std::size_t i = 0; i != ords.size(); ++i)
                for (std::size_t j = i + 1; j != ords.size(); ++j)
                  if (ords[i] != ords[j])
                    edges.emplace(std::min(ords[i], ords[j]),
                                  std::max(ords[i], ords[j]));
            };
            for (const auto &f : *summand) {
              if (!f->template is<Tensor>()) continue;
              const auto &t = f->template as<Tensor>();
              container::svector<std::size_t> ords;
              for (const auto &idx : t.braket())
                if (auto it = provenance.find(idx); it != provenance.end())
                  ords.push_back(it->second);
              add_edges(ords);
            }
            for (const auto &[a, b] : opts.nop_connections)
              if (!edges.contains({std::min(a, b), std::max(a, b)})) return;
            for (const auto &[a, b] : opts.nop_avoided_connections)
              if (edges.contains({std::min(a, b), std::max(a, b)})) return;
          }
```

This relies on `provenance` covering *all* indices of the term's tensors, not just survivors. Extend `make_provenance` (Task 6) to record every `Op` index (it already does — it iterates all ops), and in `extended_wick`'s `expand_one`, also resolve δ-bound tensor indices: after reduce, a contracted dummy pair has been renamed to a single index that appears in two tensors; that index is one of the original two (reduce keeps externals and picks one of the dummies), so it is in `provenance` already. Confirm with the test; if an index is missing, extend `resolved` to all tensor indices rather than only survivors.

- [ ] **Step 4: Run, commit**

Run: `cmake-build-debug/tests/unit/unit_tests-sequant "wick_extended"`
Expected: PASS.

```bash
bin/admin/clang-format.sh -i SeQuant/core/wick_extended.cpp tests/unit/test_wick_extended.cpp
git add SeQuant/core/wick_extended.cpp tests/unit/test_wick_extended.cpp
git commit -m "wick: extended_wick enforces operator connectivity after cumulant expansion"
```

---

### Task 8: `eta_as_delta_minus_gamma`

**Files:**
- Modify: `SeQuant/core/wick_extended.cpp`
- Test: `tests/unit/test_wick_extended.cpp`

- [ ] **Step 1: Write the failing test**

```cpp
  SECTION("extended_wick: η = δ - γ") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = extended_wick<Statistics::FermiDirac>(
        in, {.eta_as_delta_minus_gamma = true});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1} * δ{u_2;u_3} "
                                      L"- γ{u_4;u_1} * γ{u_2;u_3} "
                                      L"+ κ{u_2,u_4;u_1,u_3}:A-H-S"));
  }
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL — result still contains `η`.

- [ ] **Step 3: Implement** at the end of `extended_wick` (before the final `simplify`), as a file-local function:

```cpp
void rewrite_eta(ExprPtr &expr) {
  expr->visit(
      [](ExprPtr &e) {
        if (e->is<Tensor>() &&
            e->as<Tensor>().label() == density::hole_rdm_label()) {
          const auto &t = e->as<Tensor>();
          const Index &b = t.bra()[0], &k = t.ket()[0];
          e = make_kronecker(b, k) - density::make_rdm(b, k);
        }
      },
      /*atoms_only=*/true);
  expand(expr);
}
```

Call it `if (opts.eta_as_delta_minus_gamma) detail::rewrite_eta(result);`. Note the δ here is between two *external* active indices, so it is **not** reduced away (that is what the test expects).

- [ ] **Step 4: Run, commit**

```bash
bin/admin/clang-format.sh -i SeQuant/core/wick_extended.cpp tests/unit/test_wick_extended.cpp
git add SeQuant/core/wick_extended.cpp tests/unit/test_wick_extended.cpp
git commit -m "wick: extended_wick can spell η as δ - γ"
```

---

### Task 9: `mbpt::ref_av` dispatch and the cross-check against the core-vacuum path

**Files:**
- Modify: `SeQuant/domain/mbpt/rdm.hpp`, `rdm.cpp` (add `cumulants_to_densities`)
- Modify: `SeQuant/domain/mbpt/op.cpp:1264-1480` (`expectation_value_impl` dispatch), `op.hpp:220-239` (`QuantumNumberChange::size`), `op.cpp:164-192` (`combine`)
- Modify: `SeQuant/domain/mbpt/vac_av.cpp:170` (`op::ref_av` disables screening under `MultiProduct` unless the qns branch below is implemented — it is, so no change needed; verify)
- Test: `tests/unit/test_mbpt.cpp` new `SECTION("MRSO-MultiProduct")` after `SECTION("MRSO")`

**Interfaces:**
- Produces: `ExprPtr sequant::mbpt::decompositions::cumulants_to_densities(ExprPtr)` — replaces every κ of rank 1..3 in an expression via `cumulantN_to_density`, then expands and simplifies; throws `Exception` for rank > 3. `tensor::ref_av`/`op::ref_av` under `MultiProduct` return `extended_wick` output (cumulants kept).

- [ ] **Step 1: Write the failing tests**

```cpp
SECTION("MRSO-MultiProduct") {
  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  ctx.set(Vacuum::MultiProduct);
  auto ctx_resetter = set_scoped_default_context(ctx);

  // one-body: same expectation as the core-vacuum path at MRSO, with the
  // general-index h split into core and active pieces
  SECTION("ref_av of non-normal-ordered one-body product") {
    const Index p{L"p_1"};
    const Index q{L"p_2"};
    auto H1 = ex<Tensor>(L"h", bra{p}, ket{q}, Symmetry::Nonsymm,
                         BraKetSymmetry::Conjugate, ColumnSymmetry::Symm) *
              fcrex(p) * fannx(q);
    ExprPtr result;
    REQUIRE_NOTHROW(result = t::ref_av(H1));
    REQUIRE_THAT(result, SimplifiesTo(L"h{O_1;O_1}:N-C-S + "
                                      L"h{u_2;u_1}:N-C-S * γ{u_1;u_2}:N-C-S"));
  }

  // two-body: extended Wick + cumulant→density reproduces the core-vacuum
  // ref_av result
  SECTION("wick(H2**T2) matches the core-vacuum path") {
    auto result_mp = t::ref_av(t::h(2) * t::t(2), {.connect = {{0, 1}}});
    auto result_mp_dens =
        mbpt::decompositions::cumulants_to_densities(result_mp);

    ExprPtr result_sp;
    {
      auto sp_ctx = get_default_context();
      sp_ctx.set(Vacuum::SingleProduct);
      auto sp_resetter = set_scoped_default_context(sp_ctx);
      result_sp = t::ref_av(t::h(2) * t::t(2), {.connect = {{0, 1}}});
    }
    REQUIRE(simplify(result_mp_dens - result_sp) == ex<Constant>(0));
  }

  SECTION("topology on/off agree") {
    auto a = t::ref_av(t::h(2) * t::t(2), {.connect = {{0, 1}}});
    auto b = t::ref_av(t::h(2) * t::t(2),
                       {.connect = {{0, 1}}, .use_topology = false});
    REQUIRE(simplify(a - b) == ex<Constant>(0));
  }

  SECTION("operator-level ref_av agrees with tensor-level") {
    auto result_op = o::ref_av(o::h(2) * o::t(2));
    auto result_t = t::ref_av(t::h(2) * t::t(2), {.connect = {{0, 1}}});
    REQUIRE(simplify(result_op - result_t) == ex<Constant>(0));
  }
}  // SECTION("MRSO-MultiProduct")
```

Note on the two-body cross-check: the core-vacuum `ref_av` expresses the result in 1- and 2-body γ; `cumulants_to_densities` turns κ₂ into γ₂ − antisymmetrized γγ. The γ₂ built by `cumulant2_to_density` must hash-equal the γ₂ `ref_av` builds (same label, Antisymm? — `cumulant2_to_density` uses `rdm_symmetries` without `.perm`; `ref_av` uses `Symmetry::Antisymm`). If the difference does not simplify to zero *only* because of the `perm` attribute, make `cumulant2_to_density`/`cumulant3_to_density` build their multi-body γ with `.perm = Symmetry::Antisymm` (that is physically correct for spin-orbital RDMs and is the attribute `ref_av` already uses) — a one-line fix in `rdm.cpp`, included in this task's commit with a sentence in the message.

- [ ] **Step 2: Run to verify failure**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "[mbpt]" -c "MRSO-MultiProduct"`
Expected: FAIL — `cumulants_to_densities` undeclared, or `Exception: ... MultiProduct`.

- [ ] **Step 3: `cumulants_to_densities`** in `rdm.hpp`/`rdm.cpp`:

```cpp
/// replaces every cumulant κ_k (k ≤ 3) in @p expr by its expansion in
/// densities, then expands and simplifies
/// @throw Exception for a κ of rank > 3
ExprPtr cumulants_to_densities(ExprPtr expr);
```

```cpp
ExprPtr cumulants_to_densities(ExprPtr expr) {
  expr = expr->clone();
  expr->visit(
      [](ExprPtr &e) {
        if (!e->is<Tensor>() ||
            e->as<Tensor>().label() != density::cumulant_label())
          return;
        switch (e->as<Tensor>().rank()) {
          case 1: e = cumulant_to_density(e); break;
          case 2: e = cumulant2_to_density(e); break;
          case 3: e = cumulant3_to_density(e); break;
          default:
            throw Exception(
                "cumulants_to_densities: only cumulants of rank <= 3 are "
                "supported");
        }
      },
      /*atoms_only=*/true);
  expand(expr);
  simplify(expr);
  return expr;
}
```

- [ ] **Step 4: Dispatch in `tensor::expectation_value_impl`** (`op.cpp`). Right after the `simplify(expr);` that follows scalar extraction, before the `isr`/`spinor` lines, add:

```cpp
  if (get_default_context().vacuum() == Vacuum::MultiProduct) {
    ExtendedWickOptions opts{.full_contractions = true,
                             .use_topology = use_top};
    for (const auto& [a, b] : connect)
      opts.nop_connections.emplace_back(a, b);
    for (const auto& [a, b] : avoid)
      opts.nop_avoided_connections.emplace_back(a, b);
    auto result = extended_wick<Statistics::FermiDirac>(expr, opts);
    restore_scalars(result);  // move the lambda's definition above this block
    return result;
  }
```

(`OpConnections<int>` is a container of `std::pair<int,int>`; check `grep -n "OpConnections" SeQuant/domain/mbpt/op.hpp`.) `full_contractions` is forced true: `ref_av` with `MultiProduct` *is* the reference expectation value. Add `#include <SeQuant/core/wick_extended.hpp>`.

Note `ref_av`'s existing `full_contractions = (reference_occupied == vacuum_occupied)` logic is unchanged for `SingleProduct`.

- [ ] **Step 5: Quantum-number screening under `MultiProduct`** — the op-level path calls `can_change_qns`, which constructs `QuantumNumberChange` (throws on `MultiProduct` in `size()`) and `combine`. Make both `MultiProduct`-aware with the weakest sound bound:

`op.hpp` `QuantumNumberChange::size()`: treat `MultiProduct` like `SingleProduct` (`base_spaces.size() * 2`).

`op.cpp` `combine`: add a `MultiProduct` branch identical to the `SingleProduct` one except that for a base space that is reference-occupied *and* vacuum-unoccupied (active) both contraction types can occur:

```cpp
      const bool active =
          isr->intersection(base_spaces[i], isr->reference_occupied_space(qn)) &&
          isr->intersection(base_spaces[i], isr->vacuum_unoccupied_space(qn));
      auto ncontr_space =
          active ? qninterval_t{0, std::min(b[ann].upper(), a[cre].upper()) +
                                       std::min(b[cre].upper(), a[ann].upper())}
          : base_is_fermi_occupied ? ... (as SingleProduct)
                                   : ... (as SingleProduct);
```

where `qn = base_spaces[i].qns()`. Any other `throw` on `MultiProduct` reached by the "operator-level ref_av agrees" test gets the same treatment: the minimal branch that keeps the screen a *necessary* condition. Record each one in the commit body.

- [ ] **Step 6: Build, run the MR sections and the full mbpt/wick suites**

Run: `cmake --build cmake-build-debug --target unit_tests-sequant -j8 && cmake-build-debug/tests/unit/unit_tests-sequant "[mbpt],[wick]"`
Expected: PASS, including the pre-existing `MRSO`/`MRSF` sections (untouched `SingleProduct` path).

- [ ] **Step 7: Run the whole suite the way CI does**

Run: `cmake --build cmake-build-debug --target check-sequant -j8 2>&1 | tail -20`
Expected: every `sequant/...` test passes. If a `*/verify` or `*/dump_tree` fixture test fails, do **not** regenerate: diagnose per AGENTS.md and report.

- [ ] **Step 8: Commit (two commits)**

```bash
bin/admin/clang-format.sh -i SeQuant/domain/mbpt/rdm.hpp SeQuant/domain/mbpt/rdm.cpp SeQuant/domain/mbpt/op.hpp SeQuant/domain/mbpt/op.cpp tests/unit/test_mbpt.cpp
git add SeQuant/domain/mbpt/rdm.hpp SeQuant/domain/mbpt/rdm.cpp
git commit -m "mbpt: cumulants_to_densities rewrites every κ in an expression"
git add SeQuant/domain/mbpt/op.hpp SeQuant/domain/mbpt/op.cpp tests/unit/test_mbpt.cpp
git commit -m "mbpt: ref_av uses the extended Wick theorem under a MultiProduct vacuum"
```

---

### Task 10: Documentation

**Files:**
- Modify: `doc/developer/wick.rst` (new section before "Reducing the result")
- Modify: `doc/user/guide/context.rst:14-16` (mention `Vacuum::MultiProduct`)
- Create: `doc/examples/user/extended_wick.cpp` (compiled + run as `sequant/doc-examples/sequant_doc_example_extended_wick`)
- Modify: `doc/user/getting_started/wick.rst` (short subsection linking the example)
- Modify: `SeQuant/core/context.hpp:115` (`Vacuum::SingleReference` → `Vacuum::SingleProduct`)

- [ ] **Step 1: The compiled example**

```cpp
// doc/examples/user/extended_wick.cpp
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/io/latex.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/wick_extended.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <iostream>

int main() {
  using namespace sequant;
  // start-snippet-1
  // a multireference vocabulary: core (i), active (u), virtual (a) ...
  // and a reference that is a general state in the active space
  set_default_context(
      Context({.index_space_registry_shared_ptr = mbpt::make_mr_spaces(),
               .vacuum = Vacuum::MultiProduct}));

  // the product of two generalized-normal-ordered one-body operators
  auto expr = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));

  // its reference expectation value: a density pair plus a 2-body cumulant
  auto vev = extended_wick<Statistics::FermiDirac>(expr);
  std::wcout << to_latex(vev) << std::endl;

  // its generalized-normal-ordered expansion, cumulants truncated at rank 2
  auto gno = extended_wick<Statistics::FermiDirac>(
      expr, {.full_contractions = false, .max_cumulant_rank = 2});
  std::wcout << to_latex(gno) << std::endl;
  // end-snippet-1
  return 0;
}
```

Build and run: `cmake --build cmake-build-debug --target sequant_doc_example_extended_wick -j8 && cmake-build-debug/doc/examples/sequant_doc_example_extended_wick` — expected: two LaTeX lines, exit 0. (Target name: check `cmake --build cmake-build-debug --target help | grep extended_wick`.)

- [ ] **Step 2: Developer doc.** Insert before `Reducing the result` in `doc/developer/wick.rst`:

```rst
The extended Wick theorem (``Vacuum::MultiProduct``)
-------------------------------------------------------

Relative to a general (multiconfigurational) reference the quasiparticle picture is lost: an operator on a partially occupied
orbital is both a quasiparticle creator and annihilator. SeQuant follows Kutzelnigg and Mukherjee's *generalized normal
order*: strings are normal-ordered so that their reference expectation value vanishes, and the theorem acquires, besides the
pair contractions (now valued :math:`\gamma` for ``cre·ann`` and :math:`\eta = \delta - \gamma` for ``ann·cre``), *multi-leg*
contractions of :math:`k` creators and :math:`k` annihilators valued by the :math:`k`-body density cumulant :math:`\kappa_k`.

Rather than a second engine, this is layered on the standard one:

- the classifiers in ``SeQuant/core/op.hpp`` gain a ``MultiProduct`` branch in which the "hole" space is the registry's
  *reference-occupied* space :math:`R` and the "particle" space its *vacuum-unoccupied* space :math:`U`; the active space is
  :math:`R \cap U`, and an active ``Op`` is classified as both;
- ``contract`` emits ``γ``/``η`` (``SeQuant/core/density.hpp``) in place of the overlap ``s``; a ``γ`` over :math:`R` stands
  for :math:`\delta` on the core plus :math:`\gamma` on the active space, and likewise ``η`` over :math:`U`;
- the pair-based connectivity filters are skipped, since a cumulant can connect operators that no pair does;
- :func:`sequant::extended_wick` runs the standard theorem with partial contractions and hands each term to
  :func:`sequant::cumulant_expand`, which splits mixed-space ``γ``/``η`` and survivors into pure pieces, then groups the
  *surviving* operators into disjoint blocks (:math:`k \ge 2`, legs from at least two input operators), each valued
  :math:`\kappa_k` with the parity of the permutation that pulls its legs to the front.

The no-double-counting rule is that a block is only ever built from operators that survive the standard theorem: each extended
term descends from exactly one partial-contraction term (the one carrying its pair contractions), so no :math:`1/k!` weights
are needed and topological folding stays valid. Cumulant rank is bounded by ``ExtendedWickOptions::max_cumulant_rank``.
Spin-free evaluation is not supported for this vacuum; the hooks for it are ``WickTheorem::contraction_value`` and
``detail::block_value``/``detail::term_weight`` in ``SeQuant/core/wick_extended.hpp``. ``mbpt::ref_av`` dispatches here when
the context's vacuum is ``MultiProduct``; :func:`sequant::mbpt::decompositions::cumulants_to_densities` converts the result to
densities. Tests: ``tests/unit/test_wick_extended.cpp`` and the ``MRSO-MultiProduct`` section of ``tests/unit/test_mbpt.cpp``.
```

Also fix the existing paragraph "attributed, in the surrounding comment, to Kutzelnigg": the cited comment (`wick.hpp:1545`) is about partner re-pairing in `normalize`, not the :math:`2^n` factor — reword to "the generalized Wick's theorem for spin-free operators" without the attribution clause.

- [ ] **Step 3: User docs.** In `doc/user/guide/context.rst:15` extend the parenthetical:

```rst
(``Vacuum::Physical`` — the true, particle-free vacuum —, ``Vacuum::SingleProduct`` — a single-determinant quasiparticle vacuum
—, or ``Vacuum::MultiProduct`` — a general reference state, for which Wick's theorem takes its *extended* form with density
cumulants; see :func:`sequant::extended_wick`)
```

In `doc/user/getting_started/wick.rst`, append a subsection:

```rst
Extended Wick's theorem
~~~~~~~~~~~~~~~~~~~~~~~~~

With ``Vacuum::MultiProduct`` the reference is a general state and :func:`sequant::extended_wick` evaluates products of
generalized-normal-ordered operators into densities :math:`\gamma`, :math:`\eta` and cumulants :math:`\kappa_k`:

.. literalinclude:: /examples/user/extended_wick.cpp
   :language: cpp
   :start-after: start-snippet-1
   :end-before: end-snippet-1
   :dedent: 2
```

(Check the heading-underline convention: a few characters longer than the title.)

- [ ] **Step 4: Fix `context.hpp:115`**: `Vacuum::SingleReference` → `Vacuum::SingleProduct`.

- [ ] **Step 5: Verify docs build if Sphinx is available** (`grep -n "sphinx" doc/developer/documentation.rst` for the command); otherwise at least `grep -rn "extended_wick\|MultiProduct" doc/` to confirm every new reference resolves to something that exists.

- [ ] **Step 6: Commit**

```bash
bin/admin/clang-format.sh -i doc/examples/user/extended_wick.cpp SeQuant/core/context.hpp
git add doc/developer/wick.rst doc/user/guide/context.rst doc/user/getting_started/wick.rst doc/examples/user/extended_wick.cpp SeQuant/core/context.hpp
git commit -m "doc: the extended Wick theorem"
```

---

## Self-review notes

- Spec §A classifiers → Task 3; `contract` middle factor + protoindex guard + connectivity skip → Task 4; §B pass steps 1-3 → Task 6 (`split_mixed_spaces`, provenance), 4-5 → Task 5, 6 → Task 7, 7 → Tasks 5/6; spin-free hooks → Tasks 4 (`contraction_value`) and 5 (`block_value`/`term_weight`); §C API → Tasks 5/6/8, tensors → Task 2, rename → Task 1, mbpt dispatch → Task 9; §D tests → Tasks 3-9, docs → Task 10. Deferred optimization (engine "full but emit active survivors") intentionally has no task.
- Names used consistently: `ExtendedWickOptions`, `OpProvenance`, `cumulant_expand<S>`, `extended_wick<S>`, `density::make_rdm/make_hole_rdm/make_cumulant/rdm_from_nop`, `cumulant_to_density/cumulant2_to_density/cumulant3_to_density/cumulants_to_densities`, `WickTheorem::contraction_value`, `WickTheorem::external_indices()`.
- Review Focus: 1 → Task 6 "vacuum must be MultiProduct"; 2 → Task 5 "unbalanced survivors vanish"; 3 → Task 4 protoindex `REQUIRE_THROWS_AS`; 4 → Task 5 "max_cumulant_rank = 1"; 5 → Task 3 boson `REQUIRE_THROWS_AS`.
