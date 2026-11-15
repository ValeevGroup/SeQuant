#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>

#include <algorithm>
#include <compare>
#include <cstddef>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

namespace {

sequant::ExprPtr parse(std::wstring_view str) {
  return sequant::deserialize<sequant::ExprPtr>(str);
}

/// @return @p cmp of the reversed operands, i.e. `cmp(b,a)` for `cmp(a,b)`
std::weak_ordering reversed(std::weak_ordering cmp) { return 0 <=> cmp; }

}  // namespace

TEST_CASE("lexicographic_order", "[elements]") {
  using namespace sequant;

  const LexicographicExprCompare cmp(
      LexicographicCompareOptions{.cardinal_tensor_labels = {L"t", L"f"}});
  const auto compare = [&cmp](const ExprPtr& a, const ExprPtr& b) {
    return cmp.compare(*a, *b);
  };

  SECTION("tensors") {
    // cardinal labels first, in the given order, then the others
    CHECK(compare(parse(L"t{i1;a1}"), parse(L"f{i1;a1}")) ==
          std::weak_ordering::less);
    CHECK(compare(parse(L"f{i1;a1}"), parse(L"g{i1;a1}")) ==
          std::weak_ordering::less);
    CHECK(compare(parse(L"g{i1;a1}"), parse(L"h{i1;a1}")) ==
          std::weak_ordering::less);
    // adjoint labels have the precedence of their non-adjoint form
    CHECK(compare(ex<Tensor>(std::wstring(L"f") + adjoint_label, bra{L"i_1"},
                             ket{L"a_1"}),
                  parse(L"g{i1;a1}")) == std::weak_ordering::less);
    // then by rank
    CHECK(compare(parse(L"g{i1;a1}"), parse(L"g{i1,i2;a1,a2}")) ==
          std::weak_ordering::less);
    CHECK(compare(parse(L"g{i1;a1}"), parse(L"g{i1;a1;i2}")) ==
          std::weak_ordering::less);
    // then by indices
    CHECK(compare(parse(L"g{i1;a1}"), parse(L"g{i2;a1}")) ==
          std::weak_ordering::less);
    CHECK(compare(parse(L"g{i1;a2}"), parse(L"g{i2;a1}")) ==
          std::weak_ordering::less);
    CHECK(compare(parse(L"g{i1;a1}"), parse(L"g{i1;a2}")) ==
          std::weak_ordering::less);
    CHECK(compare(parse(L"g{i1;a1}"), parse(L"g{i1;a1}")) ==
          std::weak_ordering::equivalent);
    // the default options use the cardinal labels of the default context
    {
      auto ctx = get_default_context();
      ctx.set_cardinal_tensor_labels({L"g"});
      auto resetter = set_scoped_default_context(ctx);
      CHECK(parse(L"g{i1;a1}")->lexicographic_compare(*parse(L"f{i1;a1}")) ==
            std::weak_ordering::less);
      CHECK(LexicographicExprCompare{}(parse(L"g{i1;a1}"), parse(L"f{i1;a1}")));
    }
  }

  SECTION("scalars") {
    CHECK(compare(ex<Constant>(1), ex<Constant>(2)) ==
          std::weak_ordering::less);
    CHECK(compare(ex<Constant>(rational{1, 2}), ex<Constant>(1)) ==
          std::weak_ordering::less);
    CHECK(compare(ex<Constant>(Complex<rational>{1, 1}),
                  ex<Constant>(Complex<rational>{1, 2})) ==
          std::weak_ordering::less);
    CHECK(compare(ex<Variable>(L"a"), ex<Variable>(L"b")) ==
          std::weak_ordering::less);
    auto conj_a = ex<Variable>(L"a");
    conj_a->adjoint();
    CHECK(compare(ex<Variable>(L"a"), conj_a) == std::weak_ordering::less);
    CHECK(compare(ex<Power>(ex<Variable>(L"a"), 2),
                  ex<Power>(ex<Variable>(L"b"), 1)) ==
          std::weak_ordering::less);
    CHECK(compare(ex<Power>(ex<Variable>(L"a"), 1),
                  ex<Power>(ex<Variable>(L"a"), 2)) ==
          std::weak_ordering::less);
  }

  SECTION("composites") {
    // products: factors first, then scalar
    CHECK(compare(parse(L"t{i1;a1} g{a1;i1}"), parse(L"f{i1;a1} g{a1;i1}")) ==
          std::weak_ordering::less);
    CHECK(compare(parse(L"2 g{i1;a1}"), parse(L"g{i1;a1} g{a1;i1}")) ==
          std::weak_ordering::less);
    CHECK(compare(parse(L"g{i1;a1} g{a1;i1}"), parse(L"2 g{i1;a1} g{a1;i1}")) ==
          std::weak_ordering::less);
    // sums: summand by summand
    CHECK(compare(parse(L"t{i1;a1} + g{i1;a1}"),
                  parse(L"t{i1;a1} + h{i1;a1}")) == std::weak_ordering::less);
    CHECK(compare(parse(L"t{i1;a1} + g{i1;a1}"),
                  parse(L"t{i1;a1} + g{i1;a1} + h{i1;a1}")) ==
          std::weak_ordering::less);
    // unlike types: by type rank
    CHECK(compare(parse(L"g{i1;a1}"), parse(L"2 f{i1;a1}")) ==
          std::weak_ordering::less);
    CHECK(compare(ex<Constant>(100), ex<Variable>(L"a")) ==
          std::weak_ordering::less);
  }

  SECTION("operators") {
    const auto op1 = ex<FNOperator>(cre({L"i_1"}), ann({L"a_1"}));
    const auto op2 = ex<FNOperator>(cre({L"i_1"}), ann({L"a_2"}));
    const auto op3 = ex<FNOperator>(cre({L"i_1", L"i_2"}), ann({L"a_1"}));
    CHECK(compare(op1, op2) == std::weak_ordering::less);
    CHECK(compare(op1, op3) == std::weak_ordering::less);
    CHECK(compare(op1, op1->clone()) == std::weak_ordering::equivalent);
    const auto seq1 = ex<FNOperatorSeq>(
        FNOperatorSeq{op1.as<FNOperator>(), op2.as<FNOperator>()});
    const auto seq2 = ex<FNOperatorSeq>(
        FNOperatorSeq{op2.as<FNOperator>(), op1.as<FNOperator>()});
    CHECK(compare(seq1, seq2) == std::weak_ordering::less);
  }

  SECTION("tensor powers") {
    const auto D = parse(L"D{i1;a1}");
    const auto Dm1 = pow(D, rational{-1});
    const auto Dhalf = pow(D, rational{1, 2});
    const auto D2 = pow(D, rational{2});
    const auto Em1 = pow(parse(L"E{i1;a1}"), rational{-1});
    const auto tm1 = pow(parse(L"t{i1;a1}"), rational{-1});
    REQUIRE(Dm1->is<TensorPower>());
    // by base (cardinal labels first), then by exponent
    CHECK(compare(Dm1, Dhalf) == std::weak_ordering::less);
    CHECK(compare(Dhalf, D2) == std::weak_ordering::less);
    CHECK(compare(Dm1, D2) == std::weak_ordering::less);
    CHECK(compare(Dm1, Em1) == std::weak_ordering::less);
    CHECK(compare(tm1, Dm1) == std::weak_ordering::less);
    CHECK(compare(Dm1, pow(parse(L"D{i2;a1}"), rational{-1})) ==
          std::weak_ordering::less);
    CHECK(compare(Dm1, Dm1->clone()) == std::weak_ordering::equivalent);
    // a TensorPower and a Tensor are ordered by type, consistently
    CHECK(compare(Dm1, D) != std::weak_ordering::equivalent);
    CHECK(compare(Dm1, D) == reversed(compare(D, Dm1)));
    CHECK(compare(Dm1, D) == compare(Em1, parse(L"t{i1;a1}")));
  }

  SECTION("consistency") {
    const std::vector<ExprPtr> exprs = {
        parse(L"t{i1;a1}"),
        parse(L"f{i1;a1}"),
        parse(L"g{i1;a1}"),
        parse(L"g{i2;a1}"),
        parse(L"g{i1,i2;a1,a2}"),
        parse(L"t{i1;a1} g{a1;i1}"),
        parse(L"2 t{i1;a1} g{a1;i1}"),
        parse(L"t{i1;a1} + g{i1;a1}"),
        parse(L"g{i1;a1} + t{i1;a1}"),
        ex<Constant>(1),
        ex<Constant>(2),
        ex<Variable>(L"a"),
        ex<Power>(ex<Variable>(L"a"), 2),
        ex<FNOperator>(cre({L"i_1"}), ann({L"a_1"})),
        ex<FNOperator>(cre({L"i_1"}), ann({L"a_2"})),
        pow(parse(L"D{i1;a1}"), rational{-1}),
        pow(parse(L"D{i1;a1}"), rational{2}),
        pow(parse(L"E{i1;a1}"), rational{-1}),
    };
    for (const auto& a : exprs) {
      for (const auto& b : exprs) {
        const auto ab = compare(a, b);
        CHECK(ab == reversed(compare(b, a)));
        CHECK((ab == std::weak_ordering::equivalent) == (*a == *b));
        CHECK((ab == std::weak_ordering::equivalent) ==
              (compare(a->clone(), b) == std::weak_ordering::equivalent));
      }
    }
  }
}

TEST_CASE("sort_expr", "[algorithms]") {
  using namespace sequant;

  const LexicographicExprCompare cmp(
      LexicographicCompareOptions{.cardinal_tensor_labels = {L"t", L"f"}});

  SECTION("sum") {
    const std::vector<ExprPtr> summands = {
        parse(L"t{i1;a1}"), parse(L"f{i1;a1}"),
        parse(L"Z{i1;a1}"), parse(L"g{i1,i2;a1,a2} t{a1,a2;i1,i2}"),
        ex<Constant>(3),    ex<Variable>(L"w"),
    };
    // tensors before products before constants before variables; within the
    // product, t precedes g
    const std::wstring expected =
        ex<Sum>(ExprPtrList{summands[0], summands[1], summands[2],
                            parse(L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}"),
                            summands[4], summands[5]})
            ->to_latex();

    std::vector<std::size_t> perm(summands.size());
    std::iota(perm.begin(), perm.end(), 0);
    do {
      Sum::summands_type permuted;
      for (auto i : perm) permuted.push_back(summands[i]->clone());
      ExprPtr sum = ex<Sum>(std::move(permuted), Sum::move_only_tag{});
      const Expr* const sum_ptr = sum.get();
      sort_expr(*sum, cmp);
      // sorted in place
      REQUIRE(sum.get() == sum_ptr);
      REQUIRE(sum->to_latex() == expected);
      // the moved constant is still recognized as such
      sum.as<Sum>().append(ex<Constant>(1));
      REQUIRE(sum.as<Sum>().summands().size() == summands.size());
      REQUIRE(sum.as<Sum>().summand(4)->is<Constant>());
      REQUIRE(sum.as<Sum>().summand(4).as<Constant>().value<int>() == 4);
    } while (std::next_permutation(perm.begin(), perm.end()));
  }

  SECTION("constant summand after mutation") {
    ExprPtr sum = parse(L"t{i1;a1} + f{i1;a1}");
    // replace a summand by a constant via the mutable subexpression iterators
    *sum->begin() = ex<Constant>(2);
    sum.as<Sum>().append(ex<Constant>(3));
    REQUIRE(sum.as<Sum>().summands().size() == 2);
    REQUIRE(sum.as<Sum>().summand(0).as<Constant>().value<int>() == 5);
  }

  SECTION("nested") {
    auto expr = parse(L"(h{i1;a1} + g{i1;a1}) * (f{i1;a1} + t{i1;a1})");
    sort_expr(*expr, cmp);
    REQUIRE(
        expr->to_latex() ==
        parse(L"(t{i1;a1} + f{i1;a1}) * (g{i1;a1} + h{i1;a1})")->to_latex());
  }

  SECTION("non-commuting factors") {
    // tensors commute with the operators, and with each other, the operators
    // do not commute with each other
    const auto op1 = ex<FNOperator>(cre({L"p_3"}), ann({L"p_4"}));
    const auto op2 = ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"}));
    REQUIRE(cmp(op2, op1));
    REQUIRE(!op1->commutes_with(*op2));
    ExprPtr prod = ex<Product>(
        ExprPtrList{op1, parse(L"g{p1;p2}"), op2, parse(L"t{p3;p4}")});
    const auto orig = prod->clone();
    sort_expr(*prod, cmp);
    REQUIRE(prod->to_latex() ==
            ex<Product>(
                ExprPtrList{parse(L"t{p3;p4}"), parse(L"g{p1;p2}"), op1, op2})
                ->to_latex());
    REQUIRE(prod->hash_value() != orig->hash_value());
  }

  SECTION("tensor powers") {
    const auto factors = std::vector<ExprPtr>{
        pow(parse(L"E{i1;a1}"), rational{-1}),
        pow(parse(L"D{i1;a1}"), rational{2}),
        pow(parse(L"D{i1;a1}"), rational{-1}), parse(L"t{i1;a1}")};
    std::optional<std::wstring> expected;
    std::vector<std::size_t> perm(factors.size());
    std::iota(perm.begin(), perm.end(), 0);
    do {
      ExprPtr prod = ex<Product>(Product{});
      for (auto i : perm)
        prod.as<Product>().append(1, factors[i]->clone(), Product::Flatten::No);
      sort_expr(*prod, cmp);
      if (!expected) {
        expected = prod->to_latex();
        // the powers of D, the lower exponent first, then that of E (where
        // t, a Tensor, goes depends on the order of the types)
        const auto& sorted = prod.as<Product>().factors();
        REQUIRE(sorted.size() == 4);
        const auto position = [&sorted](const ExprPtr& f) {
          return std::ranges::find_if(
                     sorted, [&f](const ExprPtr& g) { return *g == *f; }) -
                 sorted.begin();
        };
        CHECK(position(factors[2]) + 1 == position(factors[1]));
        CHECK(position(factors[1]) + 1 == position(factors[0]));
      }
      REQUIRE(prod->to_latex() == *expected);
    } while (std::next_permutation(perm.begin(), perm.end()));
  }

  SECTION("custom comparator") {
    // reverse lexicographic order
    auto expr = parse(L"t{i1;a1} + f{i1;a1} + g{i1;a1}");
    sort_expr(*expr,
              [&cmp](const ExprPtr& a, const ExprPtr& b) { return cmp(b, a); });
    REQUIRE(expr->to_latex() ==
            parse(L"g{i1;a1} + f{i1;a1} + t{i1;a1}")->to_latex());
  }

  SECTION("result expr") {
    auto res = deserialize<ResultExpr>(L"R{i1;a1} = g{i1;a1} + t{i1;a1}");
    sort_expr(res, cmp);
    REQUIRE(res.expression()->to_latex() ==
            parse(L"t{i1;a1} + g{i1;a1}")->to_latex());
  }
}
