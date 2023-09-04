#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <codecvt>
#include <iostream>
#include <locale>
#include <string>

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
    auto expression = parse_expr(L"t{i1;a1,a2}");

    auto indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra, UnorderedEquals(std::vector<Index>{{L"i_1"}}));
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"a_1", L"a_2"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Tensor>()));

    expression = parse_expr(L"t{i1,i2;a1,a2}");

    indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra,
                 UnorderedEquals(std::vector<Index>{{L"i_1"}, {L"i_2"}}));
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"a_1", L"a_2"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Tensor>()));

    expression = parse_expr(L"t{i1,i2;a1,i1}");

    indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra, UnorderedEquals(std::vector<Index>{{L"i_2"}}));
    REQUIRE_THAT(indices.ket, UnorderedEquals(std::vector<Index>{{L"a_1"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Tensor>()));
  }
  SECTION("Product") {
    auto expression = parse_expr(L"t{i1;a1,a2} p{a2;i2}");

    auto indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra, UnorderedEquals(std::vector<Index>{{L"i_1"}}));
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"a_1", L"i_2"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Product>()));

    expression = parse_expr(L"1/8 g{a3,a4;i3,i4} t{a1,a4;i1,i4}");

    indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra,
                 UnorderedEquals(std::vector<Index>{{L"a_3"}, {L"a_1"}}));
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"i_3", L"i_1"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Product>()));
  }
  SECTION("Sum") {
    auto expression = parse_expr(L"t{i1;a2} + g{i1;a2}");

    auto indices = get_unique_indices(expression);

    REQUIRE_THAT(indices.bra, UnorderedEquals(std::vector<Index>{{L"i_1"}}));
    REQUIRE_THAT(indices.ket, UnorderedEquals(std::vector<Index>{{L"a_2"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Sum>()));

    expression = parse_expr(L"t{i1;a2} t{i1;a1} + t{i1;a1} g{i1;a2}");

    indices = get_unique_indices(expression);

    REQUIRE(indices.bra.empty());
    REQUIRE_THAT(indices.ket,
                 UnorderedEquals(std::vector<Index>{{L"a_1"}, {L"a_2"}}));
    REQUIRE(indices == get_unique_indices(expression->as<Sum>()));
  }
}
