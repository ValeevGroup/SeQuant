#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <codecvt>
#include <iostream>
#include <locale>
#include <string>
#include <thread>

namespace Catch {

// Note: For some reason this template specialization is never used. It works
// for custom types but not for sequant::Index.
template <>
struct StringMaker<sequant::Index> {
  static std::string convert(const sequant::Index& idx) {
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
  SECTION("Variable") {
    auto const expression = ex<Variable>(L"x");

    auto const indices = get_unique_indices(expression);

    REQUIRE(indices.bra.empty());
    REQUIRE(indices.ket.empty());
    REQUIRE(indices == get_unique_indices(expression->as<Variable>()));
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
