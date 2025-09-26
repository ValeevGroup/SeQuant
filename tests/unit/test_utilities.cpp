#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/singleton.hpp>
#include <SeQuant/core/utility/strong.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <iostream>
#include <locale>
#include <ranges>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

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

sequant::Tensor parse_tensor(std::wstring_view str) {
  return sequant::parse_expr(str)->as<sequant::Tensor>();
}

TEST_CASE("utilities", "[utilities]") {
  using namespace sequant;

  SECTION("get_uncontracted_indices") {
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

  SECTION("get_unique_indices") {
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

  SECTION("Singleton") {
    using namespace sequant::singleton;

    constexpr auto nthreads = 5;

    // default-constructible Singleton
    {
      std::vector<std::thread> threads;
      std::vector<int> thread_results(nthreads, -1);
      for (int t = 0; t != nthreads; ++t) {
        threads.emplace_back([&result = thread_results[t]]() {
          result = Singleton<S<EnableDefaultCtor>>::instance().s();
        });
      }
      for (auto&& thr : threads) thr.join();
      for (auto result : thread_results) CHECK(result == 0);
      CHECK_THROWS_AS(Singleton<S<EnableDefaultCtor>>::set_instance(1),
                      std::logic_error);
      CHECK(Singleton<S<EnableDefaultCtor>>::instance().s() == 0);
    }
    // non-default-constructible Singleton
    {
      {
        std::vector<std::thread> threads;
        std::vector<int> thread_results(nthreads, -1);
        for (int t = 0; t != nthreads; ++t) {
          threads.emplace_back([&result = thread_results[t]]() {
            try {
              Singleton<S<DisableDefaultCtor>>::instance().s();
            } catch (std::logic_error&) {
              result = 0;
              return;
            } catch (...) {
              result = 1;
              return;
            }
            result = 2;
          });
        }
        for (auto&& thr : threads) thr.join();
        CHECK(Singleton<S<DisableDefaultCtor>>::instance_ptr() == nullptr);
        for (auto result : thread_results) CHECK(result == 0);
        CHECK_NOTHROW(Singleton<S<DisableDefaultCtor>>::set_instance(1));
        {
          std::vector<std::thread> threads;
          std::vector<int> thread_results(nthreads, -1);
          for (int t = 0; t != nthreads; ++t) {
            threads.emplace_back([&result = thread_results[t]]() {
              result = Singleton<S<DisableDefaultCtor>>::instance().s();
            });
          }
          for (auto&& thr : threads) thr.join();
          for (auto result : thread_results) CHECK(result == 1);
          CHECK(Singleton<S<DisableDefaultCtor>>::instance().s() == 1);
        }
      }
    }
  }

  SECTION("StrongType") {
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

  SECTION("transform_expr") {
    {
      ExprPtr expr =
          parse_expr(L"- g{a1,i2;a2,i1} t{a2;i2} + 2 g{a1,i2;i1,a2} t{a2;i2}");
      container::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                             {Index{L"i_2"}, Index{L"i_1"}}};
      auto transformed_result = transform_expr(expr, idxmap);
      REQUIRE(transformed_result->is<Sum>());
      REQUIRE(transformed_result->size() == 2);
      REQUIRE_THAT(
          transformed_result,
          EquivalentTo(
              "- g{a1,i1;a2,i2} t{a2;i1} + 2 g{a1,i1;i2,a2} t{a2;i1}"));
    }
    {
      ExprPtr expr =
          parse_expr(L"- g{i2,a1;i1,a2} + 2 g{i2,a1;a2,i1} t{a2;i2}");
      container::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                             {Index{L"i_2"}, Index{L"i_1"}}};
      auto transformed_result = transform_expr(expr, idxmap);
      REQUIRE_THAT(
          transformed_result,
          EquivalentTo("- g{i1,a1;i2,a2} + 2 g{i1,a1;a2,i2} t{a2;i1}"));
    }
  }

  SECTION("TensorBlockComparator") {
    SECTION("TensorBlockEqualComparator") {
      const TensorBlockEqualComparator cmp;

      std::vector<std::tuple<std::wstring, std::wstring, bool>> test_cases = {
          {L"A{a1}", L"B{a2}", false},
          {L"t{a1}", L"t{a2}", true},
          {L"t{;i1}", L"t{;i8}", true},
          {L"t{;;p3}", L"t{;;p4}", true},
          {L"t{a1;i1;p3}", L"t{a1;i1;p4}", true},
          {L"t{a1;}", L"t{;a1}", true},
          {L"t{a1;}", L"t{;;a1}", true},
          {L"t{a1;i1}", L"t{a1;p1}", false},
          {L"t{a1;i1}", L"t{a1;i1}", true},
          {L"t{a4;i3;p2}", L"t{a2;i1;p6}", true},
          {L"I{;;a2,i2}", L"I{;;a1,i1}", true},
      };

      for (const auto& [lhs, rhs, equal] : test_cases) {
        CAPTURE(toUtf8(lhs));
        CAPTURE(toUtf8(rhs));

        const Tensor lhs_tensor = parse_expr(lhs)->as<Tensor>();
        const Tensor rhs_tensor = parse_expr(rhs)->as<Tensor>();
        REQUIRE(cmp(lhs_tensor, rhs_tensor) == equal);
      }
    }

    SECTION("TensorBlockLessThanComparator") {
      const TensorBlockLessThanComparator cmp;
      const TensorBlockEqualComparator equal_cmp;

      REQUIRE(Index(L"i_1").space() < Index(L"p_1").space());

      std::vector<std::tuple<std::wstring, std::wstring, bool>> test_cases = {
          {L"A{a1}", L"B{a2}", true},
          {L"t{a1}", L"t{a2}", false},
          {L"t{;i1}", L"t{;i8}", false},
          {L"t{;;p3}", L"t{;;p4}", false},
          {L"t{a1;i1;p3}", L"t{a1;i1;p4}", false},
          {L"t{a1;i1;p3}", L"t{a1;p4;i2}", true},
          {L"t{a1;}", L"t{;a1}", false},
          {L"t{a1;}", L"t{;;a1}", false},
          {L"t{a1;i1}", L"t{a1;p1}", true},
          {L"t{a1;i1}", L"t{a1;i1}", false},
          {L"t{a4;i3;p2}", L"t{a2;i1;p6}", false},
          {L"I{;;a2,i2}", L"I{;a1;i1}", false},
      };

      for (const auto& [lhs, rhs, less] : test_cases) {
        CAPTURE(toUtf8(lhs));
        CAPTURE(toUtf8(rhs));

        const Tensor lhs_tensor = parse_expr(lhs)->as<Tensor>();
        const Tensor rhs_tensor = parse_expr(rhs)->as<Tensor>();
        REQUIRE(cmp(lhs_tensor, rhs_tensor) == less);

        if (equal_cmp(lhs_tensor, rhs_tensor)) {
          REQUIRE(less == false);
          REQUIRE(cmp(rhs_tensor, lhs_tensor) == less);
        } else {
          REQUIRE(cmp(rhs_tensor, lhs_tensor) == !less);
        }
      }
    }
  }

  SECTION("get_used_indices") {
    for (const auto& [input, expected] :
         std::vector<std::tuple<std::wstring, std::vector<std::wstring>>>{
             {L"Var", {}},
             {L"t{a1;a2}", {L"a_1", L"a_2"}},
             {L"t{a1;a2} f{a2;a1}", {L"a_1", L"a_2"}},
             {L"t{a1;a2} - (Var * (B{p1} T{a1,a2;p1}) + t{a1;a2})",
              {L"a_1", L"a_2", L"p_1"}},
         }) {
      ExprPtr expr = parse_expr(input);
      auto indices =
          expected | std::ranges::views::transform(
                         [](const std::wstring& idx) { return Index(idx); });

      REQUIRE_THAT(get_used_indices(expr),
                   Catch::Matchers::UnorderedRangeEquals(indices));
    }
  }

  SECTION("replace") {
    SECTION("Expr") {
      for (const auto& [input_str, target_str, replacement_str, expected_str] :
           std::vector<std::tuple<std::wstring, std::wstring, std::wstring,
                                  std::wstring>>{
               {L"Var", L"t{a1;i1}", L"Test", L"Var"},
               {L"Var", L"Var", L"Test", L"Test"},
               {L"t{a1;i1} Var", L"Var", L"Test", L"t{a1;i1} Test"},
               {L"t{a1;i1} Var", L"Var", L"K{;;p1,p2}", L"t{a1;i1} K{;;p1,p2}"},
               {L"t{a1;i1} Var", L"t{a1;i1}", L"K{p1} - 1", L"(K{p1} - 1) Var"},
           }) {
        CAPTURE(toUtf8(input_str));
        CAPTURE(toUtf8(target_str));
        CAPTURE(toUtf8(replacement_str));
        CAPTURE(toUtf8(expected_str));

        ExprPtr input = parse_expr(input_str);
        const ExprPtr target = parse_expr(target_str);
        const ExprPtr replacement = parse_expr(replacement_str);

        replace(input, target, replacement);

        REQUIRE_THAT(input, EquivalentTo(expected_str));
      }
    }
    SECTION("ResultExpr") {
      // The big difference to replacing on plain expressions is that the result
      // (indices) are updated as well (if needed)
      for (const auto& [input_str, target_str, replacement_str, expected_str] :
           std::vector<std::tuple<std::wstring, std::wstring, std::wstring,
                                  std::wstring>>{
               {L"R = Var", L"t{a1;i1}", L"Test", L"R = Var"},
               {L"R = Var", L"Var", L"Test", L"R = Test"},
               {L"R{a1;i1} = t{a1;i1} Var", L"Var", L"Test",
                L"R{a1;i1} = t{a1;i1} Test"},
               {L"R{a1;i1} = t{a1;i1} Var", L"Var", L"K{;;p1,p2}",
                L"R{a1;i1;p1,p2} = t{a1;i1} K{;;p1,p2}"},
               {L"R{a1;i1} = t{a1;i1} Var", L"t{a1;i1}", L"K{p1} - 1",
                L"R{p1} = (K{p1} - 1) Var"},
           }) {
        CAPTURE(toUtf8(input_str));
        CAPTURE(toUtf8(target_str));
        CAPTURE(toUtf8(replacement_str));
        CAPTURE(toUtf8(expected_str));

        ResultExpr input = parse_result_expr(input_str);
        const ExprPtr target = parse_expr(target_str);
        const ExprPtr replacement = parse_expr(replacement_str);

        replace(input, target, replacement);

        REQUIRE_THAT(input, EquivalentTo(expected_str));
      }
    }
  }
}
