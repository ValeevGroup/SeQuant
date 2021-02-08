#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/eval_seq.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>
#include <vector>

#include "catch.hpp"

std::wostream& operator<<(std::wostream& out,
                          const sequant::utils::eval_expr& expr) {
  out << "$" << expr.tensor().to_latex() << "$";
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

auto latex_label = [](const auto& x) {
  return L'"' + x->seq_expr()->to_latex() + L'"';
};

// validates if x is constructible from tspec using parse_expr
auto validate_tensor = [](const auto& x, std::wstring_view tspec) -> bool {
  return x.to_latex() ==
         sequant::utils::parse_expr(tspec, sequant::Symmetry::antisymm)
             ->to_latex();
};

TEST_CASE("TEST_BINARIZE_EXPR", "[binarize_expr]") {
  using namespace sequant;
  using utils::eval_expr;
  using utils::parse_expr;
  using eval_seq = utils::eval_seq<size_t>;
  using utils::binarize_evxpr_range;
  using utils::binarize_flat_prod;

  SECTION("binarize_prod") {
    const auto p1 = parse_expr(
                        L"1/16 "
                        L"* g_{i3, i4}^{a3, a4}"
                        L"* t_{a1, a2}^{i3, i4}"
                        L"* t_{a3,a4}^{i1,i2}",
                        Symmetry::antisymm)
                        ->as<Product>();

    auto s1 = eval_seq{(size_t)0, {size_t{1}, size_t{2}}};

    auto binarizer = binarize_flat_prod{p1};

    auto snode1 = binarizer(s1);

    REQUIRE(snode1->scalar() == Constant{1.0 / 16});

    REQUIRE(snode1 == binarize_flat_prod{p1}());

    const auto& node1 = snode1;

    REQUIRE(validate_tensor(node1->tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(validate_tensor(node1.left()->tensor(), L"I_{a1,a2}^{a3,a4}"));

    REQUIRE(validate_tensor(node1.right()->tensor(), L"t_{a3,a4}^{i1,i2}"));

    REQUIRE(
        validate_tensor(node1.left().left()->tensor(), L"g_{i3,i4}^{a3,a4}"));

    REQUIRE(
        validate_tensor(node1.left().right()->tensor(), L"t_{a1,a2}^{i3,i4}"));

    REQUIRE(node1.right().leaf());
    REQUIRE(node1.left().left().leaf());
    REQUIRE(node1.left().right().leaf());

    const auto& s2 = eval_seq{0, {eval_seq{1, {2}}}};  // s2 = (0, (1 2))

    auto snode2 = binarizer(s2);

    REQUIRE(snode2->scalar() == Constant{1.0 / 16});

    const auto& node2 = snode2;
    REQUIRE(validate_tensor(node2->tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(validate_tensor(node2.left()->tensor(), L"g_{i3,i4}^{a3,a4}"));

    REQUIRE(validate_tensor(node2.right()->tensor(),
                            L"I_{a1,a2,a3,a4}^{i1,i2,i3,i4}"));

    REQUIRE(
        validate_tensor(node2.right().left()->tensor(), L"t_{a1,a2}^{i3,i4}"));

    REQUIRE(
        validate_tensor(node2.right().right()->tensor(), L"t_{a3,a4}^{i1,i2}"));
  }

  SECTION("binarize_evxpr_range") {
    using ranges::views::transform;

    const auto specs = {L"I1_{i1,i2}^{a1,a2}", L"I2_{i1,i2}^{a1,a2}",
                        L"I3_{i1,i2}^{a1,a2}"};

    auto tensor = [](const auto& x) {
      return parse_expr(x, Symmetry::antisymm)->template as<Tensor>();
    };

    auto ev_xpr = [](const auto& x) { return eval_expr{x}; };

    auto const summands = specs | transform(tensor) | transform(ev_xpr);

    auto const node = binarize_evxpr_range(summands);

    REQUIRE(validate_tensor(node->tensor(), L"I_{i1,i2}^{a1,a2}"));

    REQUIRE(validate_tensor(node.right()->tensor(), L"I3_{i1,i2}^{a1,a2}"));
    REQUIRE(validate_tensor(node.left()->tensor(), L"I_{i1,i2}^{a1,a2}"));

    REQUIRE(
        validate_tensor(node.left().left()->tensor(), L"I1_{i1,i2}^{a1,a2}"));

    REQUIRE(
        validate_tensor(node.left().right()->tensor(), L"I2_{i1,i2}^{a1,a2}"));

    auto const prod =
        parse_expr(L"t_{a1,a2}^{i3,i4}*g_{i3,i4}^{i1,i2}", Symmetry::antisymm);

    REQUIRE(binarize_evxpr_range(Constant{1.0 / 16},
                                 *prod | transform([&ev_xpr](auto&& x) {
                                   return ev_xpr(x->template as<Tensor>());
                                 }))
                ->scalar() == Constant{1.0 / 16});
  }

  SECTION("binarization with canonicalization") {
    auto const expr1 =
        parse_expr(L"1/4 * g_{i2,i1}^{a1,a2}", Symmetry::antisymm);

    auto const node1 = binarize_flat_prod{expr1->as<Product>()}();
    //
    // g_{i2, i1}^{a1, a2} =canonized=> -g_{i1, i2}^{a1,a2}
    //

    REQUIRE(node1->scalar() == Constant{-1.0 / 4});
    auto const expr2 = parse_expr(
        L"1/4 * t_{a1,a2}^{i3,i4} * g_{i4,i3}^{i1,i2}", Symmetry::antisymm);
    //
    // g_{i4,i3}^{i1,i2} =canonized=> -g_{i3,i4}^{i1,i2}
    //

    auto const node2 = binarize_flat_prod{expr2->as<Product>()}();

    REQUIRE(node2->scalar() == Constant{-1.0 / 4});
    REQUIRE(node2.left()->scalar() == Constant{1});
    REQUIRE(node2.right()->scalar() == Constant{1});
  }

  SECTION("debinarization") {
      // todo
  }
}
