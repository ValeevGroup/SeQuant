#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/eval_sequence.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>
#include <vector>

#include "catch.hpp"

std::wostream& operator<<(std::wostream& out,
                          const sequant::utils::eval_expr& expr) {
  out << "$" << expr.seq_expr()->to_latex() << "$";
  return out;
}

auto make_eval_seq = [](size_t start, size_t n) {
  auto children =
      ranges::views::iota(start + 1) | ranges::views::take(n - 1) |
      ranges::to<std::vector<sequant::utils::eval_sequence<size_t>>>;

  return sequant::utils::eval_sequence<size_t>{start, std::move(children)};
};

auto latex_label = [](const auto& x) {
  return L'"' + x->data().seq_expr()->to_latex() + L'"';
};

TEST_CASE("TEST_BINARIZE_EXPR", "[binarize_expr]") {
  using namespace sequant;
  using utils::eval_expr;
  using utils::parse_expr;
  using eval_seq = utils::eval_sequence<size_t>;
  using utils::binarize_prod;

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

  auto binarizer = binarize_prod{p1};

  auto snode1 = binarizer(s1);

  REQUIRE(snode1->left()->data().seq_expr()->as<Constant>() ==
          Constant{1.0 / 16});
  REQUIRE(snode1->left()->leaf());

  const auto& node1 = snode1->right();

  REQUIRE(validate_tensor(node1->data().seq_expr(), L"I_(a1,a2)^(i1,i2)"));

  REQUIRE(
      validate_tensor(node1->left()->data().seq_expr(), L"I_(a1,a2)^(a3,a4)"));

  REQUIRE(
      validate_tensor(node1->right()->data().seq_expr(), L"t_(a3,a4)^(i1,i2)"));

  REQUIRE(validate_tensor(node1->left()->left()->data().seq_expr(),
                          L"g_(i3,i4)^(a3,a4)"));

  REQUIRE(validate_tensor(node1->left()->right()->data().seq_expr(),
                          L"t_(a1,a2)^(i3,i4)"));

  REQUIRE(node1->right()->leaf());
  REQUIRE(node1->left()->left()->leaf());
  REQUIRE(node1->left()->right()->leaf());

  const auto& s2 = eval_seq{0, {eval_seq{1, {2}}}};
  // s2 = (0, (1 2))

  auto snode2 = binarizer(s2);
  REQUIRE(snode2->left()->data().seq_expr()->as<Constant>() ==
          Constant{1.0 / 16});

  const auto& node2 = snode2->right();
  REQUIRE(validate_tensor(node2->data().seq_expr(), L"I_(a1,a2)^(i1,i2)"));

  REQUIRE(
      validate_tensor(node2->left()->data().seq_expr(), L"g_(i3,i4)^(a3,a4)"));

  REQUIRE(validate_tensor(node2->right()->data().seq_expr(),
                          L"I_(a1,a2,a3,a4)^(i3,i4,i1,i2)"));

  REQUIRE(validate_tensor(node2->right()->left()->data().seq_expr(),
                          L"t_(a1,a2)^(i3,i4)"));

  REQUIRE(validate_tensor(node2->right()->right()->data().seq_expr(),
                          L"t_(a3,a4)^(i1,i2)"));
}
