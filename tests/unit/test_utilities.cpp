#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/singleton.hpp>
#include <SeQuant/core/utility/strong.hpp>

#include <codecvt>
#include <iostream>
#include <locale>
#include <string>
#include <string_view>
#include <thread>
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
        parse_tensor(L"t{i1,i2;a1,a2;x1,x2}"), parse_tensor(L"t{a1;i2;x2}"));

    std::vector<Index> expectedBra = {Index(L"i_1")};
    std::vector<Index> expectedKet = {Index(L"a_2")};
    std::vector<Index> expectedAux = {Index(L"x_1")};

    REQUIRE_THAT(bra, Catch::Matchers::UnorderedEquals(expectedBra));
    REQUIRE_THAT(ket, Catch::Matchers::UnorderedEquals(expectedKet));
    REQUIRE_THAT(aux, Catch::Matchers::UnorderedEquals(expectedAux));
  }
}

TEST_CASE("TEST GET_UNIQUE_INDICES", "[utilities]") {
  using namespace sequant;
  using namespace Catch::Matchers;

  SECTION("Constant") {
    auto const expression = parse_expr(L"5");

    auto const indices = get_unique_indices(expression);

    REQUIRE(indices.bra.empty());
    REQUIRE(indices.ket.empty());
    REQUIRE(indices == get_unique_indices(expression->as<Constant>()));
  }
  SECTION("Variable") {
    auto const expression = ex<Variable>(L"x");

    auto const indices = get_unique_indices(expression);

    REQUIRE(indices.bra.empty());
    REQUIRE(indices.ket.empty());
    REQUIRE(indices == get_unique_indices(expression->as<Variable>()));
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

namespace sequant::singleton {

enum Tag { EnableDefaultCtor, DisableDefaultCtor };
template <Tag T>
class S : public sequant::Singleton<S<T>> {
 public:
  int s() const { return s_; }

 private:
  friend class sequant::Singleton<S>;

  template <Tag U = T, typename = std::enable_if_t<U == EnableDefaultCtor>>
  S() {}

  S(int s) : s_(s) {}

  int s_ = 0;
};

}  // namespace sequant::singleton

TEST_CASE("Singleton", "[utilities]") {
  using namespace sequant;
  using namespace sequant::singleton;
  // default-constructible Singleton
  {
    std::vector<std::thread> threads;
    for (int c = 0; c != 5; ++c)
      threads.emplace_back([]() {
        CHECK(Singleton<S<EnableDefaultCtor>>::instance().s() == 0);
      });
    for (auto&& thr : threads) thr.join();
    CHECK_THROWS_AS(Singleton<S<EnableDefaultCtor>>::set_instance(1),
                    std::logic_error);
    CHECK(Singleton<S<EnableDefaultCtor>>::instance().s() == 0);
  }
  // non-default-constructible Singleton
  {
    {
      std::vector<std::thread> threads;
      for (int c = 0; c != 5; ++c)
        threads.emplace_back([]() {
          CHECK_THROWS_AS(Singleton<S<DisableDefaultCtor>>::instance().s(),
                          std::logic_error);
        });
      CHECK(Singleton<S<DisableDefaultCtor>>::instance_ptr() == nullptr);
      for (auto&& thr : threads) thr.join();
    }
    CHECK_NOTHROW(Singleton<S<DisableDefaultCtor>>::set_instance(1));
    {
      std::vector<std::thread> threads;
      for (int c = 0; c != 5; ++c)
        threads.emplace_back([]() {
          CHECK(Singleton<S<DisableDefaultCtor>>::instance().s() == 1);
        });
      for (auto&& thr : threads) thr.join();
    }
    CHECK(Singleton<S<DisableDefaultCtor>>::instance().s() == 1);
  }
}

TEST_CASE("StrongType", "[utilities]") {
  using namespace sequant::detail;

  struct A : strong_type_base<int, A> {
    using strong_type_base::strong_type_base;
  };

  struct B : strong_type_base<int, B> {
    using strong_type_base::strong_type_base;
  };

  struct C : strong_type_base<double, C> {
    using strong_type_base::strong_type_base;
  };

  A a{1};
  B b{2};
  C c0, c1;

  CHECK(a.value() == 1);
  CHECK(int(a) == 1);
  CHECK(b.value() == 2);
  CHECK(c0.value() == double(c1));
  CHECK(double(c0) == std::move(c1).value());

  struct nondefault_constructible_int {
    nondefault_constructible_int() = delete;
    nondefault_constructible_int(int i) : value(i) {}
    int value;
  };

  struct D : strong_type_base<nondefault_constructible_int, C> {
    using strong_type_base::strong_type_base;
  };

  // "D d;" does not compile, but this does
  D d(1);
}
