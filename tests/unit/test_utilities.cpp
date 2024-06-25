#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <codecvt>
#include <iostream>
#include <locale>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace Catch {

// Note: For some reason this template specialization is never used. It works
// for custom types but not for sequant::Index.
template <>
struct StringMaker<sequant::Index> {
  static std::string convert(const sequant::Index &idx) {
    using convert_type = std::codecvt_utf8<wchar_t>;
    std::wstring_convert<convert_type, wchar_t> converter;

    return converter.to_bytes(sequant::to_latex(idx));
  }
};

}  // namespace Catch

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
        parse_tensor(L"t{i1,i2;a1,a2;x1,x2}"), parse_tensor(L"t{a1;i2;x2}"));

    std::vector<Index> expectedBra = {Index(L"i_1")};
    std::vector<Index> expectedKet = {Index(L"a_2")};
    std::vector<Index> expectedAux = {Index(L"x_1")};

    REQUIRE_THAT(bra, Catch::Matchers::UnorderedEquals(expectedBra));
    REQUIRE_THAT(ket, Catch::Matchers::UnorderedEquals(expectedKet));
    REQUIRE_THAT(aux, Catch::Matchers::UnorderedEquals(expectedAux));
  }
}

TEST_CASE("get_unique_indices", "[utilities]") {
  using namespace sequant;
  using namespace Catch::Matchers;

  SECTION("Constant") {
    auto const expression = parse_expr(L"5");

    auto const indices = get_unique_indices(expression);

    REQUIRE(indices.bra.empty());
    REQUIRE(indices.ket.empty());
    REQUIRE(indices == get_unique_indices(expression->as<Constant>()));
  }
  SECTION("Tensor") {
    auto expression = parse_expr(L"t{i1;a1,a2;x1}");

    auto indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra, UnorderedEquals(std::vector<Index>{{L"i_1"}}));
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"a_1", L"a_2"}}));
    REQUIRE_THAT(indices.aux, UnorderedEquals(std::vector<Index>{{L"x_1"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Tensor>()));

    expression = parse_expr(L"t{i1,i2;a1,a2}");

    indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra,
                 UnorderedEquals(std::vector<Index>{{L"i_1"}, {L"i_2"}}));
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"a_1", L"a_2"}}));
    REQUIRE(indices.aux.size() == 0);
    REQUIRE(indices == get_unique_indices(expression->as<Tensor>()));

    expression = parse_expr(L"t{i1,i2;a1,i1}");

    indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra, UnorderedEquals(std::vector<Index>{{L"i_2"}}));
    REQUIRE_THAT(indices.ket, UnorderedEquals(std::vector<Index>{{L"a_1"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Tensor>()));
  }
  SECTION("Product") {
    auto expression = parse_expr(L"t{i1;a1,a2} p{a2;i2;x1}");

    auto indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra, UnorderedEquals(std::vector<Index>{{L"i_1"}}));
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"a_1", L"i_2"}}));
    REQUIRE_THAT(indices.aux, UnorderedEquals(std::vector<Index>{{L"x_1"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Product>()));

    expression = parse_expr(L"1/8 g{a3,a4;i3,i4;x1} t{a1,a4;i1,i4;x1}");

    indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra,
                 UnorderedEquals(std::vector<Index>{{L"a_3"}, {L"a_1"}}));
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"i_3", L"i_1"}}));
    REQUIRE(indices.aux.size() == 0);
    REQUIRE(indices == get_unique_indices(expression->as<Product>()));
  }
  SECTION("Sum") {
    auto expression = parse_expr(L"t{i1;a2;x1} + g{i1;a2;x1}");

    auto indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra, UnorderedEquals(std::vector<Index>{{L"i_1"}}));
    REQUIRE_THAT(indices.ket, UnorderedEquals(std::vector<Index>{{L"a_2"}}));
    REQUIRE_THAT(indices.aux, UnorderedEquals(std::vector<Index>{{L"x_1"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Sum>()));

    expression = parse_expr(L"t{i1;a2} t{i1;a1} + t{i1;a1} g{i1;a2}");

    indices = get_unique_indices(expression);

    REQUIRE(indices.bra.empty());
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"a_1"}, {L"a_2"}}));
    REQUIRE(indices.aux.empty());
    REQUIRE(indices == get_unique_indices(expression->as<Sum>()));
  }
}
