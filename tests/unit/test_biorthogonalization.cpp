#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/biorthogonalization.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <SeQuant/domain/mbpt/convention.hpp>

#include <string>
#include <vector>

TEST_CASE("biorthogonalization", "[Biorthogonalization]") {
  using namespace sequant;

  auto resetter = set_scoped_default_context(
      Context(mbpt::make_mr_spaces(), Vacuum::SingleProduct));

  SECTION("plain ExprPtr") {
    const std::vector<std::wstring> inputs = {
        L"t{a1;i1}", L"S{i1,i2;a1,a2}:S g{a1,a2;i1,i2}",
        L"S{a1,a2,a3;i1,i2,i3}:S t{i1,i2,i3;a1,a2,a3}"};
    const std::vector<std::wstring> expected_outputs = {
        L"1/2 t{a1;i1}",
        L"S{i1,i2;a1,a2}:S 1/6 (2 g{a1,a2;i1,i2} + g{a2,a1;i1,i2})",
        // TODO: this is not an ideal test case as this simply sums to zero
        L"S{a1,a2,a3;i1,i2,i3}:S 1/20 ("
        "-7 t{i_1,i_2,i_3;a_2,a_3,a_1}:N-C-S "
        "- 7 t{i_1,i_2,i_3;a_3,a_1,a_2}:N-C-S "
        "- t{i_1,i_2,i_3;a_2,a_1,a_3}:N-C-S "
        "- t{i_1,i_2,i_3;a_3,a_2,a_1}:N-C-S "
        "- t{i_1,i_2,i_3;a_1,a_3,a_2}:N-C-S "
        "+ 17t{i_1,i_2,i_3;a_1,a_2,a_3}:N-C-S)"};

    REQUIRE(inputs.size() == expected_outputs.size());

    for (std::size_t i = 0; i < inputs.size(); ++i) {
      CAPTURE(i);

      ExprPtr input_expr = parse_expr(inputs.at(i));

      auto externals = external_indices(input_expr);

      ExprPtr actual = biorthogonal_transform(input_expr, externals);

      REQUIRE_THAT(actual, EquivalentTo(expected_outputs.at(i)));
    }
  }

  SECTION("ResultExpr") {
    std::vector<std::vector<std::wstring>> inputs = {
        {L"R{a1;i1} = t{a1;i1}"},
        {L"R{i1,i2;a1,a2} = S{i1,i2;a1,a2}:S g{a1,a2;i1,i2}"}};
    std::vector<std::vector<std::wstring>> expected_outputs = {
        {L"R{a1;i1} = 1/2 t{a1;i1}"},
        {L"R{a1,a2;i1,i2} = S{i1,i2;a1,a2}:S 1/6 (2 g{a1,a2;i1,i2} + "
         L"g{a2,a1;i1,i2})"}};

    REQUIRE(inputs.size() == expected_outputs.size());

    for (std::size_t i = 0; i < inputs.size(); ++i) {
      CAPTURE(i);

      REQUIRE(inputs.at(i).size() == expected_outputs.at(i).size());

      container::svector<ResultExpr> expressions;
      container::svector<ResultExpr> expected;
      for (std::size_t k = 0; k < inputs.at(i).size(); ++k) {
        ResultExpr parsed = parse_result_expr(inputs.at(i).at(k));
        expressions.push_back(parsed);

        expected.push_back(parse_result_expr(expected_outputs.at(i).at(k)));
      }

      biorthogonal_transform(expressions);

      REQUIRE(expected.size() == expressions.size());

      for (std::size_t k = 0; k < expected.size(); ++k) {
        CAPTURE(k);
        REQUIRE_THAT(expected.at(k), EquivalentTo(expressions.at(k)));
      }
    }
  }
}
