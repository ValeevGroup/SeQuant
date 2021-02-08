#include <SeQuant/domain/factorize/factorize.hpp>
#include <SeQuant/domain/utils/eval_seq.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>

#include "catch.hpp"

TEST_CASE("TEST_FACTORIZE", "[factorize]") {
  using namespace sequant;
  using utils::eval_expr;
  using utils::parse_expr;
  using eval_seq = utils::eval_sequence<size_t>;
  using factorize::detail::product_binarizer;

  SECTION("product_binarizer") {
    const auto p1 = parse_expr(
                        L"1/16 * g_(i3, i4)^(a3, a4) * t_(a1, a2)^(i3, i4) * "
                        L"t_(a3,a4)^(i1,i2)",
                        Symmetry::antisymm)
                        ->as<Product>();

    auto seq1 = eval_seq{0, {1, 2}};
    auto binarizer = product_binarizer{p1};

    auto node =
        utils::binarize_eval_sequence<size_t, eval_expr>(seq1, binarizer);
  }
}
