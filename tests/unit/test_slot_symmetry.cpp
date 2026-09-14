#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/slot_symmetry.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <memory>
#include <vector>

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

  // ---- Task 0.7: n-column generalization + maximal-subset ----

  SECTION("3-column triples: g*t both ColumnSymm -> ColumnGroup {0,1,2}") {
    // g{a_1,a_2,a_3; a_4,a_5,a_6} ColumnSymm * t{a_4,a_5,a_6; i_1,i_2,i_3}
    // ColumnSymm: all 3 result columns inherit symmetrically -> {0,1,2} sign
    // +1.
    Tensor g{L"g",
             bra(IndexList{L"a_1", L"a_2", L"a_3"}),
             ket(IndexList{L"a_4", L"a_5", L"a_6"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    Tensor t{L"t",
             bra(IndexList{L"a_4", L"a_5", L"a_6"}),
             ket(IndexList{L"i_1", L"i_2", L"i_3"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(g) * ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    auto const& ss = (*node).slot_symmetry();
    REQUIRE(ss.column_groups.size() == 1);
    REQUIRE(ss.column_groups[0].cols ==
            container::svector<std::size_t>{0, 1, 2});
    REQUIRE(ss.column_groups[0].sign == 1);
    REQUIRE(ss.bra_groups.empty());
    REQUIRE(ss.ket_groups.empty());
  }

  SECTION(
      "sub-group: g{a1,a2,a3;a4,a5,a6} * f{a6;a7} -> ColumnGroup {0,1} only"
      " (column 2 ket from ColumnNonsymm operand)") {
    // g has ColumnGroup {0,1,2} (ColumnSymm); f is ColumnSymmetry::Nonsymm
    // (no column group). Contraction on a_6. Result bra=[a_1,a_2,a_3],
    // ket=[a_4,a_5,a_7] (all virtual, ascending label order: a_4<a_5<a_7).
    // Columns 0=(a_1,a_4) and 1=(a_2,a_5) both trace to g (column-grouped);
    // column 2=(a_3,a_7): a_7 traces to f (NOT column-grouped) -> col 2
    // skipped. Maximal-subset emits ColumnGroup {0,1} only; not {0,1,2}; not
    // empty.
    Tensor g{L"g",
             bra(IndexList{L"a_1", L"a_2", L"a_3"}),
             ket(IndexList{L"a_4", L"a_5", L"a_6"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    Tensor f{L"f",
             bra(IndexList{L"a_6"}),
             ket(IndexList{L"a_7"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(g) * ex<Tensor>(f));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    auto const& ss = (*node).slot_symmetry();
    REQUIRE(ss.column_groups.size() == 1);
    REQUIRE(ss.column_groups[0].cols == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.column_groups[0].sign == 1);
  }

  // ---- Task 0.8 Part A: adjoint() free function ----

  SECTION("adjoint(s): bra_group moves to ket_group, column_group preserved") {
    // Input: one bra_group {0,1} sign -1, one column_group {0,1} sign +1,
    // no ket_group. Output: column_group preserved, bra_groups empty,
    // ket_group {0,1} sign -1 (the original bra_group).
    SlotSymmetry s;
    s.column_groups.push_back(
        SlotSymmetry::ColumnGroup{container::svector<std::size_t>{0, 1}, 1});
    s.bra_groups.push_back(
        SlotSymmetry::SlotGroup{container::svector<std::size_t>{0, 1}, -1});

    auto r = sequant::adjoint(s);

    REQUIRE(r.column_groups.size() == 1);
    REQUIRE(r.column_groups[0].cols == container::svector<std::size_t>{0, 1});
    REQUIRE(r.column_groups[0].sign == 1);
    REQUIRE(r.bra_groups.empty());
    REQUIRE(r.ket_groups.size() == 1);
    REQUIRE(r.ket_groups[0].slots == container::svector<std::size_t>{0, 1});
    REQUIRE(r.ket_groups[0].sign == -1);
  }

  SECTION("adjoint(adjoint(s)) == s (involution)") {
    SlotSymmetry s;
    s.column_groups.push_back(
        SlotSymmetry::ColumnGroup{container::svector<std::size_t>{0, 1}, 1});
    s.bra_groups.push_back(
        SlotSymmetry::SlotGroup{container::svector<std::size_t>{0, 1}, -1});
    REQUIRE(sequant::adjoint(sequant::adjoint(s)) == s);
  }

  SECTION("Adjoint EvalExpr node carries swapped descriptor") {
    // t{a_1,a_2; i_1}: bra_rank=2 >= 2, ket_rank=1 < 2.
    // from_leaf_tensor: bra_group {0,1} sign -1 (Antisymm), no ket_group,
    // no column_group (ncols = min(2,1) = 1 < 2).
    // After adjoint: ket_group {0,1} sign -1, no bra_group, no column_group.
    Tensor t{L"t",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1"}),
             Symmetry::Antisymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    Tensor t_adj = t;
    t_adj.adjoint();
    // Adjoint marker must have been applied (BraKetSymmetry::Nonsymm).
    REQUIRE(!t_adj.label().empty());
    REQUIRE(t_adj.label().back() == adjoint_label);

    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t_adj));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    REQUIRE(tree->op_type() == EvalOp::Adjoint);

    auto const& ss = tree->slot_symmetry();
    REQUIRE(!ss.empty());
    // bra_groups: ket_groups of bare = empty (ket_rank=1 < 2).
    REQUIRE(ss.bra_groups.empty());
    // ket_groups: bra_groups of bare = [{0,1}, -1].
    REQUIRE(ss.ket_groups.size() == 1);
    REQUIRE(ss.ket_groups[0].slots == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.ket_groups[0].sign == -1);
    // column_groups: preserved from bare = empty (ncols=1 < 2).
    REQUIRE(ss.column_groups.empty());
  }

  // ---- Task 0.8 Part B: CSE descriptor probe ----
  //
  // Observation: the CSE round-trips through Expr (to_expr -> cse_placeholder
  // rebuild -> binarize). The definition tree is re-binarized from the
  // original expression so its root carries the deduced descriptor. The
  // reference nodes in parent trees are fresh Nonsymm leaves (built as
  // ColumnSymmetry::Nonsymm in common_subexpression_elimination.hpp) and
  // therefore have empty descriptors.
  //
  // This is a Phase-1 concern: consumers must look up the definition tree to
  // learn the CSE intermediate's symmetry. No production fix in Phase 0.

  SECTION(
      "CSE: definition tree root has descriptor; reference leaf is empty"
      " (Phase-1 TODO)") {
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    IndexSpaceRegistry registry;
    registry.add("a", 0b01);
    registry.add("i", 0b10);
    *get_default_context().mutable_index_space_registry() = registry;

    // g{a_1,a_2; a_3,a_4} * t{a_3,a_4; i_1,i_2}: column-symmetric PPL product.
    // The product tree carries ColumnGroup {0,1} at its root.
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

    auto binarizer = [](auto&& expr) {
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      return binarize(std::forward<decltype(expr)>(expr));
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    };

    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    EvalNode<EvalExpr> tree1 = binarize(ex<Tensor>(g) * ex<Tensor>(t));
    EvalNode<EvalExpr> tree2 = binarize(ex<Tensor>(g) * ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    // Pre-CSE: both roots have non-empty descriptor (column group {0,1}).
    REQUIRE(!tree1->slot_symmetry().empty());
    REQUIRE(tree1->slot_symmetry().column_groups.size() == 1u);

    std::vector<EvalNode<EvalExpr>> exprs;
    exprs.push_back(std::move(tree1));
    exprs.push_back(std::move(tree2));

    opt::eliminate_common_subexpressions(exprs, binarizer);

    // After CSE: exprs[0] is the definition tree (CSE1 = g*t). The definition
    // is built by re-running binarize on the original product expression, so
    // the root node retains the deduced ColumnGroup {0,1} descriptor.
    REQUIRE(exprs.size() == 3u);
    REQUIRE(!exprs[0]->slot_symmetry().empty());
    REQUIRE(exprs[0]->slot_symmetry().column_groups.size() == 1u);

    // exprs[1] and exprs[2] are the modified reference trees. Each is now a
    // Product node whose effective content is a CSE placeholder leaf. The
    // placeholder is constructed with ColumnSymmetry::Nonsymm (see
    // common_subexpression_elimination.hpp), so from_leaf_tensor returns
    // empty. The root descriptor of the reference tree is therefore empty.
    //
    // Phase-1 TODO: propagate the definition tree's descriptor to consumers
    // that reference the CSE intermediate as a leaf.
    REQUIRE(exprs[1]->slot_symmetry().empty());
    REQUIRE(exprs[2]->slot_symmetry().empty());
  }

  // ---- Review-findings regression tests (C1, I1, I2, Symm+1) ----

  SECTION(
      "C1: 3-tensor outer-product keeps two column groups separate"
      " (no false group merge)") {
    // L{a_1,a_2;a_3,a_4} and R{a_5,a_6;a_7,a_8} both ColumnSymm.
    // L*R is an outer product with TWO column groups: {0,1} (from L, cols
    // (a_1,a_3)/(a_2,a_4)) and {2,3} (from R, cols (a_5,a_7)/(a_6,a_8)).
    // P{a_3,a_4,a_7,a_8;i_1,i_2,i_3,i_4} ColumnSymm contracts L*R's ket.
    //
    // Bug (pre-fix): deduce_slot_symmetry for (L*R)*P sees all four L*R
    // columns as "column_grouped=true" (bool flattened from ANY group), so
    // all four cluster under the single (bra_supplier=LR, ket_supplier=P)
    // pair -> one merged ColumnGroup {0,1,2,3}.
    //
    // Fix: the 4-tuple key (bra_supplier, bra_group_idx, ket_supplier,
    // ket_group_idx) keeps columns from L*R group 0 separate from group 1 ->
    // two ColumnGroups {0,1} and {2,3}.
    Tensor L{L"L",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"a_3", L"a_4"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    Tensor R{L"R",
             bra(IndexList{L"a_5", L"a_6"}),
             ket(IndexList{L"a_7", L"a_8"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    Tensor P{L"P",
             bra(IndexList{L"a_3", L"a_4", L"a_7", L"a_8"}),
             ket(IndexList{L"i_1", L"i_2", L"i_3", L"i_4"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(L) * ex<Tensor>(R) * ex<Tensor>(P));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    auto const& ss = (*node).slot_symmetry();

    // With the C1 fix: two separate ColumnGroups {0,1} and {2,3}.
    // Without the fix: one merged ColumnGroup {0,1,2,3} (false positive).
    REQUIRE(ss.column_groups.size() == 2);
    for (auto const& cg : ss.column_groups) {
      REQUIRE(cg.cols.size() == 2);
      REQUIRE(cg.sign == 1);
    }
  }

  SECTION("I1 merge_indices: scalar*tensor result descriptor is empty") {
    // Under merge_indices mode, make_tensor_wo_symmetries puts all indices
    // into aux (bra_rank=0, ket_rank=0). The passthrough of the operand's
    // descriptor (whose slot positions refer to non-aux bra/ket) would be
    // meaningless. The fix guards on rank match and returns empty instead.
    Tensor g{L"g",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Symm};
    BinarizationOptions opts;
    opts.merge_indices = true;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Constant>(rational{1, 2}) * ex<Tensor>(g),
                         IndexSet{}, opts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    // Without the fix the column group from g would be copied verbatim
    // (column positions 0,1 now meaninglessly refer to aux slots).
    REQUIRE((*node).slot_symmetry().empty());
  }

  SECTION("Guard 1: proto-indexed externals -> empty descriptor (soundness)") {
    // A product whose operands carry proto-indexed externals (here a_1<i_1>,
    // a_1<i_2>) violates the flat index->slot trace's bijectivity assumption:
    // a proto-index also lives inside another slot, coupling slots the trace
    // treats as independent, which can produce a WRONG symmetry. Rather than
    // trust the model where its precondition fails, deduce_slot_symmetry
    // declines (empty) on any proto-indexed participant. (The full_label()
    // keying stays for the non-proto paths; this SECTION formerly asserted an
    // inherited bra_group, before the guard was added -- an intended change.)
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    IndexSpaceRegistry registry;
    registry.add("a", 0b01);
    registry.add("i", 0b10);
    *get_default_context().mutable_index_space_registry() = registry;

    Index i1(L"i_1"), i2(L"i_2");
    Index a1_i1(L"a_1", {i1});  // a_1<i_1>
    Index a1_i2(L"a_1", {i2});  // a_1<i_2>

    Tensor L_t{L"L",
               bra(Index::index_vector{a1_i1, a1_i2}),
               ket(IndexList{L"a_3", L"a_4"}),
               Symmetry::Antisymm,
               BraKetSymmetry::Nonsymm,
               ColumnSymmetry::Symm};
    Tensor R_t{L"R",
               bra(IndexList{L"a_3", L"a_4"}),
               ket(IndexList{L"a_5"}),
               Symmetry::Nonsymm,
               BraKetSymmetry::Nonsymm,
               ColumnSymmetry::Nonsymm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(L_t) * ex<Tensor>(R_t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    // Guard 1: proto-indexed externals -> deduction declines.
    REQUIRE((*node).slot_symmetry().empty());
  }

  SECTION("repeated identical factors: emergent symmetry MISSED (documented)") {
    // A{a_1;i_1} * A{a_2;i_2}: the two factors are the same tensor core, so the
    // outer product has an emergent column symmetry {a_1,i_1} <-> {a_2,i_2}
    // (topological equivalence of the two column bundles). The flat rules
    // cannot see it, so the descriptor is empty here. This is a SAFE false
    // negative -- the deduced group is always a subset of the true symmetry,
    // never a superset -- so it is left UNGUARDED (a guard would only drop
    // correct claims elsewhere). The emergent symmetry awaits the
    // graph-canonicalization deducer (Phase 0.5). This SECTION documents the
    // known incompleteness.
    Tensor A1{L"A",
              bra(IndexList{L"a_1"}),
              ket(IndexList{L"i_1"}),
              Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm,
              ColumnSymmetry::Symm};
    Tensor A2{L"A",
              bra(IndexList{L"a_2"}),
              ket(IndexList{L"i_2"}),
              Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm,
              ColumnSymmetry::Symm};
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(A1) * ex<Tensor>(A2));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE((*node).slot_symmetry().empty());
  }

  SECTION("Symm leaf: bra_group and ket_group with sign +1") {
    // Symmetry::Symm -> sign +1 in both bra_group and ket_group.
    // This path was previously untested (only Antisymm sign -1 was covered).
    Tensor t{L"t",
             bra(IndexList{L"a_1", L"a_2"}),
             ket(IndexList{L"i_1", L"i_2"}),
             Symmetry::Symm,
             BraKetSymmetry::Nonsymm,
             ColumnSymmetry::Nonsymm};
    EvalExpr ee{t};
    auto const& ss = ee.slot_symmetry();

    REQUIRE(ss.bra_groups.size() == 1);
    REQUIRE(ss.bra_groups[0].slots == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.bra_groups[0].sign == 1);
    REQUIRE(ss.ket_groups.size() == 1);
    REQUIRE(ss.ket_groups[0].slots == container::svector<std::size_t>{0, 1});
    REQUIRE(ss.ket_groups[0].sign == 1);
    REQUIRE(ss.column_groups.empty());
  }
}
