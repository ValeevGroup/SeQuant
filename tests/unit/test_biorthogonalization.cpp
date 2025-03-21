#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/biorthogonalization.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <string>
#include <vector>

TEST_CASE("biorthogonalization", "[Biorthogonalization]") {
  using namespace sequant;

  SECTION("plain ExprPtr") {
    const std::vector<std::wstring> inputs = {
        L"t{a1;i1}", L"S{i1,i2;a1,a2}:S g{a1,a2;i1,i2}"};
    const std::vector<std::wstring> expected_outputs = {
        L"1/2 t{a1;i1}",
        L"S{i1,i2;a1,a2}:S 1/6 (2 g{a1,a2;i1,i2} + g{a2,a1;i1,i2})",
    };

    REQUIRE(inputs.size() == expected_outputs.size());

    for (std::size_t i = 0; i < inputs.size(); ++i) {
      CAPTURE(i);

      ExprPtr input_expr = parse_expr(inputs.at(i));

      std::optional<ExprPtr> symmetrizer = pop_tensor(input_expr, L"S");

      auto externals = external_indices(input_expr);

      ExprPtr actual = biorthogonal_transform(input_expr, externals);

      if (symmetrizer.has_value()) {
        actual = symmetrizer.value() * actual;
        simplify(actual);
      }

      REQUIRE_THAT(actual, EquivalentTo(expected_outputs.at(i)));
    }
  }
}
