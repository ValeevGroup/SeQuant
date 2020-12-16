#include <SeQuant/domain/utils/eval_sequence.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>
#include <SeQuant/domain/utils/prod_binarizer.hpp>
#include <vector>

#include "catch.hpp"

TEST_CASE("TEST_PROD_BINARIZER", "[prod_binarizer]") {
  using namespace sequant;
  using utils::eval_expr;
  using utils::parse_expr;
  using eval_seq = utils::eval_sequence<size_t>;
  using utils::prod_binarizer;

  auto make_eval_seq = [](size_t start, size_t n) {
    auto children = ranges::views::iota(start + 1) |
                    ranges::views::take(n - 1) |
                    ranges::to<std::vector<eval_seq>>;

    return eval_seq{start, std::move(children)};
  };

  // helper lambda
  // validates if x is constructible from tspec using parse_expr
  auto validate_tensor = [](const auto& x, std::wstring_view tspec) -> bool {
    return x->to_latex() == parse_expr(tspec, Symmetry::antisymm)->to_latex();
  };
  // helper lambda end

  const auto p1 = parse_expr(
                      L"1/16 "
                      L"* g_(i3, i4)^(a3, a4)"
                      L"* t_(a1, a2)^(i3, i4)"
                      L"* t_(a3,a4)^(i1,i2)",
                      Symmetry::antisymm)
                      ->as<Product>();

  auto s1 = make_eval_seq(0, p1.size());

  auto binarizer = prod_binarizer{p1};

  auto node = utils::binarize_eval_sequence<size_t, eval_expr>(s1, binarizer);

  REQUIRE(validate_tensor(node->data().seq_expr(), L"I_(a1,a2)^(i1,i2)"));

  REQUIRE(
      validate_tensor(node->left()->data().seq_expr(), L"I_(a1,a2)^(a3,a4)"));

  REQUIRE(
      validate_tensor(node->right()->data().seq_expr(), L"t_(a3,a4)^(i1,i2)"));

  REQUIRE(validate_tensor(node->left()->left()->data().seq_expr(),
                          L"g_(i3,i4)^(a3,a4)"));

  REQUIRE(validate_tensor(node->left()->right()->data().seq_expr(),
                          L"t_(a1,a2)^(i3,i4)"));

  REQUIRE(node->right()->leaf());
  REQUIRE(node->left()->left()->leaf());
  REQUIRE(node->left()->right()->leaf());
}
