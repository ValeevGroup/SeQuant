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

template <typename Os>
Os& operator<<(Os& os, sequant::utils::eval_expr::eval_op op) {
  switch (op) {
    case sequant::utils::eval_expr::eval_op::Id:
      os << "eval: Id";
      break;
    case sequant::utils::eval_expr::eval_op::Sum:
      os << "eval: Sum";
      break;
    case sequant::utils::eval_expr::eval_op::Prod:
      os << "eval: Prod";
      break;
  }
  return os;
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

// validates if x is constructible from tspec using parse_expr
auto validate_tensor = [](const auto& x, std::wstring_view tspec) -> bool {
  return x->to_latex() ==
         sequant::utils::parse_expr(tspec, sequant::Symmetry::antisymm)
             ->to_latex();
};

TEST_CASE("TEST_BINARIZE_EXPR", "[binarize_expr]") {
  using namespace sequant;
  using utils::eval_expr;
  using utils::parse_expr;
  using eval_seq = utils::eval_sequence<size_t>;
  using utils::binarize_flat_prod;

  SECTION("binarize_prod") {
    const auto p1 = parse_expr(
                        L"1/16 "
                        L"* g_{i3, i4}^{a3, a4}"
                        L"* t_{a1, a2}^{i3, i4}"
                        L"* t_{a3,a4}^{i1,i2}",
                        Symmetry::antisymm)
                        ->as<Product>();

    auto s1 = make_eval_seq(0, p1.size());

    auto binarizer = binarize_flat_prod{p1};

    auto snode1 = binarizer(s1);

    REQUIRE(snode1->data().scalar() == Constant{1.0 / 16});

    const auto& node1 = snode1;

    REQUIRE(validate_tensor(node1->data().seq_expr(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(validate_tensor(node1->left()->data().seq_expr(),
                            L"I_{a1,a2}^{a3,a4}"));

    REQUIRE(validate_tensor(node1->right()->data().seq_expr(),
                            L"t_{a3,a4}^{i1,i2}"));

    REQUIRE(validate_tensor(node1->left()->left()->data().seq_expr(),
                            L"g_{i3,i4}^{a3,a4}"));

    REQUIRE(validate_tensor(node1->left()->right()->data().seq_expr(),
                            L"t_{a1,a2}^{i3,i4}"));

    REQUIRE(node1->right()->leaf());
    REQUIRE(node1->left()->left()->leaf());
    REQUIRE(node1->left()->right()->leaf());

    const auto& s2 = eval_seq{0, {eval_seq{1, {2}}}};  // s2 = (0, (1 2))

    auto snode2 = binarizer(s2);

    REQUIRE(snode2->data().scalar() == Constant{1.0 / 16});

    const auto& node2 = snode2;
    REQUIRE(validate_tensor(node2->data().seq_expr(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(validate_tensor(node2->left()->data().seq_expr(),
                            L"g_{i3,i4}^{a3,a4}"));

    REQUIRE(validate_tensor(node2->right()->data().seq_expr(),
                            L"I_{a1,a2,a3,a4}^{i1,i2,i3,i4}"));

    REQUIRE(validate_tensor(node2->right()->left()->data().seq_expr(),
                            L"t_{a1,a2}^{i3,i4}"));

    REQUIRE(validate_tensor(node2->right()->right()->data().seq_expr(),
                            L"t_{a3,a4}^{i1,i2}"));
  }

  SECTION("binarize_evxpr_range") {
    using ranges::views::transform;

    const auto specs = {L"I1_{i1,i2}^{a1,a2}", L"I2_{i1,i2}^{a1,a2}",
                        L"I3_{i1,i2}^{a1,a2}"};

    auto tensor = [](const auto& x) {
      return parse_expr(x, Symmetry::antisymm)->template as<Tensor>();
    };

    auto ev_xpr = [](const auto& x) { return eval_expr{x}; };

    const auto summands =
        specs | transform(tensor) | transform(ev_xpr);  //| transform(eseq);

    auto node = utils::binarize_evxpr_range(summands);

    REQUIRE(validate_tensor(node->data().seq_expr(), L"I_{i1,i2}^{a1,a2}"));

    REQUIRE(validate_tensor(node->right()->data().seq_expr(),
                            L"I3_{i1,i2}^{a1,a2}"));
    REQUIRE(
        validate_tensor(node->left()->data().seq_expr(), L"I_{i1,i2}^{a1,a2}"));

    REQUIRE(validate_tensor(node->left()->left()->data().seq_expr(),
                            L"I1_{i1,i2}^{a1,a2}"));

    REQUIRE(validate_tensor(node->left()->right()->data().seq_expr(),
                            L"I2_{i1,i2}^{a1,a2}"));
  }

  SECTION("binarization of edge cases") {
    auto expr = parse_expr(L"1/4 * g_{i1,i2}^{a1,a2}", Symmetry::antisymm);
    auto root = utils::binarize_flat_prod{expr->as<Product>()}(eval_seq{0});
    std::wcout << "phase = " << root->data().phase().value()
               << "  scalar = " << root->data().scalar().value()
               << "  op = " << root->data().op() << std::endl;
  }

  // SECTION("debinarization") {
  //   auto prod1 = parse_expr(L"1/8 * g_(a1,a2)^(a3,a4) * t_(a3,a4)^(i1,i2)",
  //                           Symmetry::antisymm);
  //   auto prod1_node =
  //       binarize_flat_prod{prod1->as<Product>()}(eval_seq{0, {1}});
  //   auto prod1_denode = utils::debinarize_eval_expr(prod1_node);
  //   std::wcout << "prod1 = $" << prod1->to_latex() << "$\n"
  //              << "de_prod1 = $" << prod1_denode->to_latex() << "$\n";
  // }
}
