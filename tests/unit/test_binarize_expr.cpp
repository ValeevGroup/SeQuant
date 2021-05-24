#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/eval_seq.hpp>
#include <SeQuant/domain/utils/parse_expr.hpp>
#include <vector>

#include "catch.hpp"

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
  using utils::binarize_expr;

  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  SECTION("binarize_prod") {
    const auto p1 = parse_expr(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}",
        Symmetry::antisymm);

    auto node1 = binarize_expr(p1);

    REQUIRE(node1->scalar() == Constant{1.0 / 16});

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

    auto node2p = Product{p1->as<Product>().scalar(), {}};
    node2p.append(p1->at(0));
    node2p.append(ex<Product>(Product{p1->at(1), p1->at(2)}));

    auto const node2 = binarize_expr(ex<Product>(node2p));

    REQUIRE(node2->scalar() == Constant{1.0 / 16});

    REQUIRE(validate_tensor(node2->tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(validate_tensor(node2.left()->tensor(), L"g_{i3,i4}^{a3,a4}"));

    REQUIRE(validate_tensor(node2.right()->tensor(),
                            L"I_{a1,a2,a3,a4}^{i1,i2,i3,i4}"));

    REQUIRE(
        validate_tensor(node2.right().left()->tensor(), L"t_{a1,a2}^{i3,i4}"));

    REQUIRE(
        validate_tensor(node2.right().right()->tensor(), L"t_{a3,a4}^{i1,i2}"));
  }

  SECTION("binarization with canonicalization") {
    auto const expr1 =
        parse_expr(L"1/4 * g_{i2,i1}^{a1,a2}", Symmetry::antisymm);

    auto const node1 = binarize_expr(expr1);
    //
    // g_{i2, i1}^{a1, a2} =canonized=> -g_{i1, i2}^{a1,a2}
    //

    REQUIRE(node1->scalar() == Constant{-1.0 / 4});
    auto const expr2 = parse_expr(
        L"1/4 * t_{a1,a2}^{i3,i4} * g_{i4,i3}^{i1,i2}", Symmetry::antisymm);
    //
    // g_{i4,i3}^{i1,i2} =canonized=> -g_{i3,i4}^{i1,i2}
    //

    auto const node2 = binarize_expr(expr2);

    REQUIRE(node2->scalar() == Constant{-1.0 / 4});
    REQUIRE(node2.left()->scalar() == Constant{1});
    REQUIRE(node2.right()->scalar() == Constant{1});
  }

  SECTION("debinarization") {
    using utils::binarize_expr;
    const auto p1 = parse_expr(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}",
        Symmetry::antisymm);
    auto node = binarize_expr(p1);
    // auto lin = linearize(node);
    // std::wcout << "lin = " << lin->to_latex() << std::endl;
  }

  SECTION("binarize_expr") {
    using utils::binarize_expr;
    auto const prod = parse_expr(
                          L"1/4 * g_{i3,i4}^{a3,a4}"
                          L"* t_{a1}^{i3}"
                          L"* t_{a2}^{i4}"
                          L"* t_{a3}^{i1}"
                          L"* t_{a4}^{i2}",
                          Symmetry::antisymm)
                          ->as<Product>();

    auto fac1 = ex<Product>(1.0 / 4.0, ExprPtrList{prod.factor(0)->clone(),
                                                   prod.factor(1)->clone()});
    auto fac2 = ex<Product>(
        ExprPtrList{prod.factor(2)->clone(), prod.factor(3)->clone()});
    auto fac3 = prod.factor(4)->clone();

    auto expr_prod = Product(1.0, ExprPtrList{});
    expr_prod.append(fac1);
    expr_prod.append(fac2);
    expr_prod.append(fac3);

    auto expr = ex<Product>(expr_prod);

    // std::wcout << binarize_expr(expr).digraph<std::wstring>(
    //                   [](auto const& node) {
    //                     return L"\"$" + node->tensor().to_latex() + L"$\"";
    //                   })
    //            << std::endl;

    std::wcout << sequant::utils::linearize_expr(expr)->to_latex() << std::endl;

    // auto node = binarize_expr(p1);

    // std::wcout << "$" << utils::debinarize_eval_expr(node)->to_latex() <<
    // "$\n";

    // std::wcout << node->tensor().to_latex() << "\n";
    // std::wcout << node.tikz<std::wstring>(
    //                   [](auto const& node) {
    //                     return L"$" + node->tensor().to_latex() + L"$";
    //                   },
    //                   [](auto const& node) { return L"draw,circle"; })
    //            << std::endl;
  }
}
