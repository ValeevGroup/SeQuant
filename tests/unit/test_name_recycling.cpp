#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/optimize/name_recycling.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace sequant;

using Tree = EvalNode<EvalExpr>;

const io::serialization::DeserializationOptions nonsymm_opts{
    .def_perm_symm = Symmetry::Nonsymm,
    .def_braket_symm = BraKetSymmetry::Nonsymm,
    .def_col_symm = ColumnSymmetry::Nonsymm};

std::vector<Tree> parse_trees(const std::vector<std::wstring>& inputs) {
  std::vector<Tree> trees;
  for (const std::wstring& current : inputs) {
    trees.push_back(binarize(deserialize<ResultExpr>(current, nonsymm_opts)));
  }
  return trees;
}

std::vector<ResultExpr> parse_results(const std::vector<std::wstring>& inputs) {
  std::vector<ResultExpr> results;
  for (const std::wstring& current : inputs) {
    results.push_back(deserialize<ResultExpr>(current, nonsymm_opts));
  }
  return results;
}

std::vector<ResultExpr> to_results(const std::vector<Tree>& trees) {
  std::vector<ResultExpr> results;
  for (const Tree& current : trees) {
    if (current->is_tensor()) {
      results.emplace_back(current->as_tensor(), to_expr(current));
    } else {
      results.emplace_back(current->as_variable(), to_expr(current));
    }
  }
  return results;
}

auto binarizer = [](auto&& expr) {
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  return binarize(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
};

/// Considers every result whose label starts with one of the given prefixes as
/// recyclable
auto recyclable_prefixes(std::vector<std::wstring> prefixes) {
  return [prefixes = std::move(prefixes)](std::size_t, const Tree& tree) {
    const std::wstring label(tree->is_tensor() ? tree->as_tensor().label()
                                               : tree->as_variable().label());
    return std::ranges::any_of(prefixes, [&](const std::wstring& prefix) {
      return label.starts_with(prefix);
    });
  };
}

}  // namespace

TEST_CASE("intermediate name recycling", "[optimize]") {
  using namespace sequant;
  using opt::NameRecyclingOptions;
  using opt::ShapeCompatibility;

  IndexBasisRegistry registry;
  registry.add("a", 0b01);
  registry.add("i", 0b10);
  auto ctx_resetter = set_scoped_default_context(
      get_default_context_snapshot().set(std::move(registry)));

  auto check = [](const std::vector<std::wstring>& inputs,
                  const std::vector<std::wstring>& outputs,
                  const NameRecyclingOptions<Tree>& opts,
                  std::size_t expected_values, std::size_t expected_names) {
    CAPTURE(inputs);
    std::vector<Tree> trees = parse_trees(inputs);

    opt::NameRecyclingResult stats =
        opt::recycle_intermediate_names(trees, binarizer, opts);

    REQUIRE(to_results(trees) == parse_results(outputs));
    REQUIRE(stats.num_values == expected_values);
    REQUIRE(stats.num_names == expected_names);
  };

  SECTION("scalars") {
    const NameRecyclingOptions<Tree> opts{.is_recyclable =
                                              recyclable_prefixes({L"X"})};

    // Disjoint lifetimes
    check({L"X1 = A B", L"R1 = X1 C", L"X2 = D E", L"R2 = X2 F"},
          {L"X1 = A B", L"R1 = X1 C", L"X1 = D E", L"R2 = X1 F"}, opts, 2, 1);

    // Overlapping lifetimes
    check({L"X1 = A B", L"X2 = D E", L"R = X1 X2"},
          {L"X1 = A B", L"X2 = D E", L"R = X1 X2"}, opts, 2, 2);

    // Last use of X1 is in the definition of X2
    check({L"X1 = A B", L"X2 = X1 C", L"R = X2 D"},
          {L"X1 = A B", L"X2 = X1 C", L"R = X2 D"}, opts, 2, 2);

    // Names are compacted
    check({L"X3 = A B", L"R1 = X3 C", L"X7 = D E", L"R2 = X7 F"},
          {L"X1 = A B", L"R1 = X1 C", L"X1 = D E", L"R2 = X1 F"}, opts, 2, 1);
  }

  SECTION("accumulation") {
    // By default, later definitions add to the existing value
    check({L"X1 = A B", L"X1 = C D", L"R1 = X1 E", L"X2 = F G", L"R2 = X2 H"},
          {L"X1 = A B", L"X1 = C D", L"R1 = X1 E", L"X1 = F G", L"R2 = X1 H"},
          {.is_recyclable = recyclable_prefixes({L"X"})}, 2, 1);

    // A value that is also written by a non-recyclable tree keeps its name
    check({L"X1 = A B", L"X1 = C D", L"R1 = X1 E", L"X3 = F G", L"R2 = X3 H"},
          {L"X1 = A B", L"X1 = C D", L"R1 = X1 E", L"X2 = F G", L"R2 = X2 H"},
          {.is_recyclable =
               [](std::size_t pos, const Tree& tree) {
                 return pos != 1 && recyclable_prefixes({L"X"})(pos, tree);
               }},
          1, 1);
  }

  SECTION("redefinitions") {
    const NameRecyclingOptions<Tree> opts{
        .is_recyclable = recyclable_prefixes({L"CSE"}),
        .starts_new_value = [](std::size_t, const Tree&) { return true; }};

    // The same label is used for unrelated values
    check({L"CSE1 = A B", L"CSE2 = C D", L"R1 = CSE1 CSE2", L"CSE2 = E F",
           L"R2 = CSE2 G"},
          {L"CSE1 = A B", L"CSE2 = C D", L"R1 = CSE1 CSE2", L"CSE1 = E F",
           L"R2 = CSE1 G"},
          opts, 3, 2);

    // Redefinition while the previous value is still in use in the same tree
    check({L"CSE1 = A B", L"CSE1 = CSE1 C", L"R = CSE1 D"},
          {L"CSE1 = A B", L"CSE2 = CSE1 C", L"R = CSE2 D"}, opts, 2, 2);
  }

  SECTION("tensors") {
    const std::vector<std::wstring> inputs = {
        L"I1{;;a1,i1} = A{a1;i2} B{i2;i1}",
        L"R1{a1;i1} = I1{;;a1,i2} F{i2;i1}",
        L"I2{;;i1,a1} = C{i1;a2} D{a2;a1}",
        L"R2{a1;i1} = I2{;;i1,a2} G{a2;a1}",
    };

    // Same block
    check({L"I1{;;a1,i1} = A{a1;i2} B{i2;i1}",
           L"R1{a1;i1} = I1{;;a1,i2} F{i2;i1}",
           L"I2{;;a1,i1} = C{a1;a2} D{a2;i1}",
           L"R2{a1;i1} = I2{;;a2,i1} G{a2;a1}"},
          {L"I1{;;a1,i1} = A{a1;i2} B{i2;i1}",
           L"R1{a1;i1} = I1{;;a1,i2} F{i2;i1}",
           L"I1{;;a1,i1} = C{a1;a2} D{a2;i1}",
           L"R2{a1;i1} = I1{;;a2,i1} G{a2;a1}"},
          {.is_recyclable = recyclable_prefixes({L"I"})}, 2, 1);

    // Same index spaces in different order. Tensors of different blocks are
    // different objects and may thus share the same label.
    check(inputs,
          {L"I1{;;a1,i1} = A{a1;i2} B{i2;i1}",
           L"R1{a1;i1} = I1{;;a1,i2} F{i2;i1}",
           L"I1{;;i1,a1} = C{i1;a2} D{a2;a1}",
           L"R2{a1;i1} = I1{;;i1,a2} G{a2;a1}"},
          {.is_recyclable = recyclable_prefixes({L"I"})}, 2, 2);
    check(inputs,
          {L"I1{;;a1,i1} = A{a1;i2} B{i2;i1}",
           L"R1{a1;i1} = I1{;;a1,i2} F{i2;i1}",
           L"I1{;;a1,i1} = C{i1;a2} D{a2;a1}",
           L"R2{a1;i1} = I1{;;a2,i1} G{a2;a1}"},
          {.is_recyclable = recyclable_prefixes({L"I"}),
           .shape_compatibility = ShapeCompatibility::IndexMultiset},
          2, 1);

    // Index spaces may only be permuted within bra, ket and aux
    check(
        {
            L"I1{a1;i1} = A{a1;i2} B{i2;i1}",
            L"R1{a1;i1} = I1{a1;i2} F{i2;i1}",
            L"I2{i1;a1} = C{i1;a2} D{a2;a1}",
            L"R2{a1;i1} = I2{i1;a2} G{a2;a1}",
        },
        {
            L"I1{a1;i1} = A{a1;i2} B{i2;i1}",
            L"R1{a1;i1} = I1{a1;i2} F{i2;i1}",
            L"I1{i1;a1} = C{i1;a2} D{a2;a1}",
            L"R2{a1;i1} = I1{i1;a2} G{a2;a1}",
        },
        {.is_recyclable = recyclable_prefixes({L"I"}),
         .shape_compatibility = ShapeCompatibility::IndexMultiset},
        2, 2);

    // Different index spaces
    check(
        {
            L"I1{;;a1,i1} = A{a1;i2} B{i2;i1}",
            L"R1{a1;i1} = I1{;;a1,i2} F{i2;i1}",
            L"I2{;;a1,a2} = C{a1;i1} D{i1;a2}",
            L"R2{a1;i1} = I2{;;a1,a2} G{a2;i1}",
        },
        {
            L"I1{;;a1,i1} = A{a1;i2} B{i2;i1}",
            L"R1{a1;i1} = I1{;;a1,i2} F{i2;i1}",
            L"I1{;;a1,a2} = C{a1;i1} D{i1;a2}",
            L"R2{a1;i1} = I1{;;a1,a2} G{a2;i1}",
        },
        {.is_recyclable = recyclable_prefixes({L"I"}),
         .shape_compatibility = ShapeCompatibility::IndexMultiset},
        2, 2);

    // Overlapping lifetimes with different blocks
    check(
        {
            L"I3{;;a1,i1} = A{a1;i2} B{i2;i1}",
            L"I4{;;a1,a2} = C{a1;i1} D{i1;a2}",
            L"R{a1;i1} = I3{;;a2,i1} I4{;;a1,a2}",
        },
        {
            L"I1{;;a1,i1} = A{a1;i2} B{i2;i1}",
            L"I1{;;a1,a2} = C{a1;i1} D{i1;a2}",
            L"R{a1;i1} = I1{;;a2,i1} I1{;;a1,a2}",
        },
        {.is_recyclable = recyclable_prefixes({L"I"})}, 2, 2);
  }

  SECTION("conjugated and adjoint uses") {
    auto leaf_objects = [](const Tree& tree) {
      std::vector<ExprPtr> objects;
      (void)to_expr(tree, [&](const ExprPtr& leaf) {
        objects.push_back(leaf);
        return leaf;
      });
      return objects;
    };

    // Conjugated power of a recyclable scalar
    {
      std::vector<Tree> trees = parse_trees({L"X5 = A B"});
      Power power(ex<Variable>(L"X5"), 2);
      power.conjugate();
      trees.push_back(binarize(ResultExpr(
          Variable(L"R"),
          ex<Product>(ExprPtrList{ex<Power>(power), ex<Variable>(L"C")}))));

      opt::recycle_intermediate_names(
          trees, binarizer, {.is_recyclable = recyclable_prefixes({L"X"})});

      REQUIRE(trees.front()->as_variable().label() == L"X1");
      const std::vector<ExprPtr> leaves = leaf_objects(trees.back());
      auto it = std::ranges::find_if(
          leaves, [](const ExprPtr& leaf) { return leaf->is<Power>(); });
      REQUIRE(it != leaves.end());
      CHECK((*it)->as<Power>().conjugated());
      CHECK((*it)->as<Power>().base()->as<Variable>().label() == L"X1");
    }

    // Adjoint of a recyclable tensor
    {
      std::vector<Tree> trees = parse_trees({L"X5{a1;i1} = A{a1;i2} B{i2;i1}"});
      Tensor adjoint(L"X5", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
                     BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
      adjoint.adjoint();
      REQUIRE(is_adjoint_label(adjoint.label()));
      trees.push_back(binarize(ResultExpr(
          Tensor(L"R", bra{L"i_1"}, ket{L"i_2"}, Symmetry::Nonsymm,
                 BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm),
          ex<Product>(ExprPtrList{ex<Tensor>(adjoint),
                                  deserialize(L"C{a1;i2}", nonsymm_opts)}))));

      opt::recycle_intermediate_names(
          trees, binarizer, {.is_recyclable = recyclable_prefixes({L"X"})});

      REQUIRE(trees.front()->as_tensor().label() == L"X1");
      const std::vector<ExprPtr> leaves = leaf_objects(trees.back());
      auto it = std::ranges::find_if(leaves, [](const ExprPtr& leaf) {
        return leaf->is<Tensor>() &&
               is_adjoint_label(leaf->as<Tensor>().label());
      });
      REQUIRE(it != leaves.end());
      CHECK(strip_adjoint_label((*it)->as<Tensor>().label()) == L"X1");
      CHECK((*it)->as<Tensor>().bra_rank() == 1);
      CHECK(std::ranges::begin((*it)->as<Tensor>().bra())->label() == L"i_1");
    }
  }

  SECTION("extra compatibility") {
    const std::vector<std::wstring> inputs = {L"X1 = A B", L"R1 = X1 C",
                                              L"X2 = D E", L"R2 = X2 F"};
    check(inputs, inputs,
          {.is_recyclable = recyclable_prefixes({L"X"}),
           .extra_compatible = [](const Tree&, const Tree&) { return false; }},
          2, 2);
  }

  SECTION("repeated regions") {
    const std::vector<std::wstring> inputs = {L"X1 = A B", L"R1 = X1 C",
                                              L"X2 = D E", L"R2 = X2 F"};

    // X1 is used inside a loop it isn't defined in
    check(inputs, inputs,
          {.is_recyclable = recyclable_prefixes({L"X"}),
           .repeated_regions = {{1, 3}}},
          2, 2);

    // Definition and uses inside the same loop
    check(inputs, {L"X1 = A B", L"R1 = X1 C", L"X1 = D E", L"R2 = X1 F"},
          {.is_recyclable = recyclable_prefixes({L"X"}),
           .repeated_regions = {{0, 4}}},
          2, 1);
  }

  SECTION("reserved names") {
    // X1 is not a recyclable intermediate
    check({L"X5 = A X1", L"R = X5 B"}, {L"X2 = A X1", L"R = X2 B"},
          {.is_recyclable = recyclable_prefixes({L"X5"})}, 1, 1);
  }

  SECTION("after common subexpression elimination") {
    std::vector<Tree> trees = parse_trees(
        {L"R1 = (A B) C", L"R2 = D (A B)", L"R3 = (E F) G", L"R4 = H (E F)"});

    std::vector<std::size_t> cse_positions =
        opt::eliminate_common_subexpressions(trees, binarizer);
    REQUIRE(cse_positions.size() == 2);

    opt::NameRecyclingResult stats = opt::recycle_intermediate_names(
        trees, binarizer,
        {.is_recyclable =
             [&](std::size_t pos, const Tree&) {
               return std::ranges::find(cse_positions, pos) !=
                      cse_positions.end();
             },
         .starts_new_value = [](std::size_t, const Tree&) { return true; },
         .label_gen =
             [](std::wstring_view, std::size_t counter) {
               return L"CSE" + std::to_wstring(counter);
             }});

    REQUIRE(to_results(trees) ==
            parse_results({L"CSE1 = A B", L"R1 = CSE1 C", L"R2 = D CSE1",
                           L"CSE1 = E F", L"R3 = CSE1 G", L"R4 = H CSE1"}));
    REQUIRE(stats.num_values == 2);
    REQUIRE(stats.num_names == 1);
  }
}
