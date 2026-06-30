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

  // ---- Task 0.5: Sum intersection ----

  SECTION("intersect: same 2-column group in both -> retained") {
    SlotSymmetry a, b;
    a.column_groups.push_back(
        SlotSymmetry::ColumnGroup{container::svector<std::size_t>{0, 1}, 1});
    b.column_groups.push_back(
        SlotSymmetry::ColumnGroup{container::svector<std::size_t>{0, 1}, 1});
    auto result = sequant::intersect(a, b);
    REQUIRE(result.column_groups.size() == 1);
    REQUIRE(result.column_groups[0].cols ==
            container::svector<std::size_t>{0, 1});
    REQUIRE(result.column_groups[0].sign == 1);
  }

  SECTION("intersect: group in a but b is empty -> dropped") {
    SlotSymmetry a, b;
    a.column_groups.push_back(
        SlotSymmetry::ColumnGroup{container::svector<std::size_t>{0, 1}, 1});
    // b has no column groups
    auto result = sequant::intersect(a, b);
    REQUIRE(result.column_groups.empty());
    REQUIRE(result.empty());
  }

  SECTION("Sum of PPL+PPL products -> column group {0,1} retained (positive)") {
    // Two PPL-style products: g1*t1 and g2*t2, both yielding r{a_1,a_2;i_1,i_2}
    // with a 2-column group.  Their sum must also carry the column group.
    Tensor g1{L"g",
              bra(IndexList{L"a_1", L"a_2"}),
              ket(IndexList{L"a_3", L"a_4"}),
              Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm,
              ColumnSymmetry::Symm};
    Tensor t1{L"t",
              bra(IndexList{L"a_3", L"a_4"}),
              ket(IndexList{L"i_1", L"i_2"}),
              Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm,
              ColumnSymmetry::Symm};
    Tensor g2{L"g2",
              bra(IndexList{L"a_1", L"a_2"}),
              ket(IndexList{L"a_5", L"a_6"}),
              Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm,
              ColumnSymmetry::Symm};
    Tensor t2{L"t2",
              bra(IndexList{L"a_5", L"a_6"}),
              ket(IndexList{L"i_1", L"i_2"}),
              Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm,
              ColumnSymmetry::Symm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(g1) * ex<Tensor>(t1) +
                         ex<Tensor>(g2) * ex<Tensor>(t2));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    auto const& ss = (*node).slot_symmetry();
    REQUIRE(ss.column_groups.size() == 1);
    REQUIRE(ss.column_groups[0].cols == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.column_groups[0].sign == 1);
  }

  SECTION("Sum of PPL+g.f products -> column group empty (negative)") {
    // First summand: PPL g{a_1,a_2;a_3,a_4} * t{a_3,a_4;i_1,i_2}
    //   -> 2-column group {0,1}
    // Second summand: g2{a_1,a_2;i_1,a_5} * f{a_5;i_2} (g.f-like)
    //   -> no column group (ket of column 1 traces to f, not column-grouped)
    // Sum result must have empty column group.
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
    Tensor g2{L"g2",
              bra(IndexList{L"a_1", L"a_2"}),
              ket(IndexList{L"i_1", L"a_5"}),
              Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm,
              ColumnSymmetry::Symm};
    Tensor f{L"f",
             bra(IndexList{L"a_5"}),
             ket(IndexList{L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(g) * ex<Tensor>(t) +
                         ex<Tensor>(g2) * ex<Tensor>(f));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE((*node).slot_symmetry().column_groups.empty());
  }

  // ---- Task 0.6: bra-only / ket-only group inheritance ----

  SECTION("antisymm bra group inherited whole into result bra (primary)") {
    // A{a_1,a_2; i_3,i_4} Antisymm contracted on ket with B{i_3,i_4; i_1}
    // Nonsymm. Result r{a_1,a_2; i_1}: bra traces whole to A's antisymm bra
    // group -> result bra_group {0,1} sign -1; ket rank 1 -> no ket group;
    // min(2,1)=1 -> no column group.
    Tensor A{L"A",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_3", L"i_4"}),
             Symmetry::Antisymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    Tensor B{L"B",
             bra(IndexList{L"i_3", L"i_4"}),
             ket(IndexList{L"i_1"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(A) * ex<Tensor>(B));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    auto const& ss = (*node).slot_symmetry();
    REQUIRE(ss.bra_groups.size() == 1);
    REQUIRE(ss.bra_groups[0].slots == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.bra_groups[0].sign == -1);
    REQUIRE(ss.ket_groups.empty());
    REQUIRE(ss.column_groups.empty());
  }

  SECTION("antisymm bra group partially contracted -> no bra_group (break)") {
    // A{a_1,a_2,a_3; i_3,i_4} Antisymm: bra_group {0,1,2} sign -1.
    // B{i_3,i_4; a_3,i_1} Nonsymm: contracts i_3,i_4 (A ket) and a_3 (A
    // bra). Result r{a_1,a_2; i_1}: only a_1,a_2 survive in result bra, not
    // a_3. Whole-group guard fires: no bra_group emitted.
    Tensor A{L"A",
             bra(IndexList{L"a_1", L"a_2", L"a_3"}),
             ket(IndexList{L"i_3", L"i_4"}),
             Symmetry::Antisymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    Tensor B{L"B",
             bra(IndexList{L"i_3", L"i_4"}),
             ket(IndexList{L"a_3", L"i_1"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(A) * ex<Tensor>(B));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    auto const& ss = (*node).slot_symmetry();
    REQUIRE(ss.bra_groups.empty());
  }
}
