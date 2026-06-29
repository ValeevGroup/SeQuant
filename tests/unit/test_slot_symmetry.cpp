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
}
