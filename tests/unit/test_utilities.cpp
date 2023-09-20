#include "catch.hpp"

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <string>
#include <string_view>
#include <utility>
#include <vector>

sequant::Tensor parse_tensor(std::wstring_view str) {
  return sequant::parse_expr(str)->as<sequant::Tensor>();
}

TEST_CASE("TEST GET_UNCONCTRACTED_INDICES", "[utilities]") {
  using namespace sequant;

  SECTION("dot_product") {
    std::vector<std::pair<std::wstring, std::wstring>> inputs = {
        {L"t{}", L"t{}"},
        {L"t{i1}", L"t{;i1}"},
        {L"t{;i1}", L"t{i1}"},
        {L"t{i1;a1}", L"t{a1;i1}"},
        {L"t{i1;a1;x1}", L"t{a1;i1;x1}"},
        {L"t{;i1;x1}", L"t{i1;;x1}"},
        {L"t{i1;;x1}", L"t{;i1;x1}"},
        {L"t{;;x1}", L"t{;;x1}"},
    };

    for (auto [left, right] : inputs) {
      auto [bra, ket, aux] =
          get_uncontracted_indices(parse_tensor(left), parse_tensor(right));

      REQUIRE(bra.size() == 0);
      REQUIRE(ket.size() == 0);
      REQUIRE(aux.size() == 0);
    }
  }

  SECTION("partial_contraction") {
    auto [bra, ket, aux] = get_uncontracted_indices<std::vector<Index>>(
        parse_tensor(L"t{i1,i2;a1,a2;x1,x2}"), parse_tensor(L"t{a1;i2;x2}")
	);

    std::vector<Index> expectedBra = {Index(L"i_1")};
    std::vector<Index> expectedKet = {Index(L"a_2")};
    std::vector<Index> expectedAux = {Index(L"x_1")};

    REQUIRE_THAT(bra, Catch::UnorderedEquals(expectedBra));
    REQUIRE_THAT(ket, Catch::UnorderedEquals(expectedKet));
    REQUIRE_THAT(aux, Catch::UnorderedEquals(expectedAux));
  }
}
