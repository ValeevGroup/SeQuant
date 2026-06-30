#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/slot_symmetry.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <memory>

TEST_CASE("slot_symmetry", "[slot_symmetry]") {
  using namespace sequant;

  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  SECTION("default descriptor is empty") {
    SlotSymmetry ss{};
    REQUIRE(ss.empty());
  }

  SECTION("operator== on two default descriptors") {
    SlotSymmetry ss1{};
    SlotSymmetry ss2{};
    REQUIRE(ss1 == ss2);
  }

  SECTION("carrier present on a Nonsymm leaf EvalExpr - empty") {
    // Explicitly Nonsymm in every axis (deserialize defaults ColumnSymmetry to
    // Symm, which is no longer empty after leaf translation).
    Tensor t{L"t",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    EvalExpr ee{t};
    REQUIRE(ee.slot_symmetry().empty());
  }

  SECTION("non-empty SlotSymmetry not equal to empty") {
    SlotSymmetry empty{};

    SlotSymmetry nonempty{};
    nonempty.bra_groups.push_back(
        SlotSymmetry::SlotGroup{container::svector<std::size_t>{0, 1}, 1});

    REQUIRE(!(empty == nonempty));
    REQUIRE(!nonempty.empty());
  }

  // ---- Task 0.2: leaf descriptor translation ----

  SECTION("leaf closed-shell g{a,b;i,j} (ColumnSymmetry::Symm) -> column grp") {
    // Closed-shell two-electron integral: column-symmetric, bra/ket Nonsymm.
    Tensor g{L"g",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    EvalExpr ee{g};
    auto const& ss = ee.slot_symmetry();

    REQUIRE(!ss.empty());
    REQUIRE(ss.column_groups.size() == 1);
    REQUIRE(ss.column_groups[0].cols == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.column_groups[0].sign == 1);
    REQUIRE(ss.bra_groups.empty());
    REQUIRE(ss.ket_groups.empty());
  }

  SECTION("leaf bra-antisymm -> bra_group sign -1, no column/ket group") {
    // Symmetry::Antisymm implies ColumnSymmetry::Symm (the invariant). The
    // descriptor records both: a column group AND an antisymmetric bra/ket
    // group. The plan's bra-antisymm acceptance row is exercised by the
    // bra_group sign here; a separate fully-bra-only (no column) case is
    // produced via product deduction in a later task.
    Tensor t{L"t",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Antisymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    EvalExpr ee{t};
    auto const& ss = ee.slot_symmetry();

    REQUIRE(ss.bra_groups.size() == 1);
    REQUIRE(ss.bra_groups[0].slots == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.bra_groups[0].sign == -1);
    REQUIRE(ss.ket_groups.size() == 1);
    REQUIRE(ss.ket_groups[0].slots == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.ket_groups[0].sign == -1);
  }

  SECTION("leaf fully Nonsymm -> empty descriptor") {
    Tensor t{L"f",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    EvalExpr ee{t};
    REQUIRE(ee.slot_symmetry().empty());
  }

  SECTION("from_leaf_tensor matches the leaf EvalExpr descriptor") {
    Tensor g{L"g",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    EvalExpr ee{g};
    REQUIRE(from_leaf_tensor(g) == ee.slot_symmetry());
  }

  // ---- Task 0.3: product column-group inheritance (PPL / giant) ----

  SECTION("PPL g{a,b;c,d} t{c,d;i,j} -> 2-column group {0,1} sign +1") {
    // Both factors column-symmetric; the contraction pairs (c,d) symmetrically,
    // so the result column symmetry is inherited.
    Tensor g{L"g",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"a_3", L"a_4"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    Tensor t{L"t",
             bra(IndexList{L"a_3", L"a_4"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(g) * ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    auto const& ss = (*node).slot_symmetry();

    REQUIRE(ss.column_groups.size() == 1);
    REQUIRE(ss.column_groups[0].cols == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.column_groups[0].sign == 1);
    REQUIRE(ss.bra_groups.empty());
    REQUIRE(ss.ket_groups.empty());
  }

  SECTION("scalar * tensor inherits the tensor operand descriptor") {
    Tensor g{L"g",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    // 1/2 * g{a,b;i,j}: scalar*tensor product node keeps g's column group.
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Constant>(rational{1, 2}) * ex<Tensor>(g));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    auto const& ss = (*node).slot_symmetry();
    REQUIRE(ss.column_groups.size() == 1);
    REQUIRE(ss.column_groups[0].cols == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.column_groups[0].sign == 1);
  }

  SECTION("product of two Nonsymm factors -> empty descriptor") {
    Tensor f{L"f",
             bra(IndexList{L"a_1"}),
             ket(IndexList{L"a_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    Tensor h{L"h",
             bra(IndexList{L"a_2"}),
             ket(IndexList{L"i_1"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(f) * ex<Tensor>(h));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE((*node).slot_symmetry().empty());
  }

  // ---- Task 0.4: column-group break (g.f gate negative) ----

  SECTION("g{a,b;i,k} f{k,j} -> no column group (gate negative)") {
    // g is a 2-column column-symmetric factor; f is a 1-column Fock-like
    // factor.  One ket index of g (a_3) is contracted with f's bra.  The
    // result r{a_1,a_2; i_1,i_2} has 2 columns, but column 1's ket (i_2)
    // traces to f which carries no ColumnGroup (ncols < 2 guard in
    // from_leaf_tensor), so column_grouped == false and the break rule fires.
    Tensor g{L"g",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"a_3"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    Tensor f{L"f",
             bra(IndexList{L"a_3"}),
             ket(IndexList{L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(g) * ex<Tensor>(f));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE((*node).slot_symmetry().column_groups.empty());
  }
}
