#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/format.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <format>

int main() {
  using namespace sequant;

  // start-snippet-1
  // Tensor is a leaf Expr: a labeled tensor with a bra and a ket index set
  auto f = ex<Tensor>("f", bra{"i_1"}, ket{"a_1"});
  auto g = ex<Tensor>("g", bra{"i_1", "i_2"}, ket{"a_1", "a_2"});

  // ExprPtr arithmetic builds Product and Sum nodes; both are associative and
  // flatten automatically, so nesting them never grows the tree depth
  auto product = f * g;      // a Product of 2 factors
  auto sum = f * g + g + f;  // a Sum of 3 summands, not a Sum of a Sum and f

  SEQUANT_ASSERT(product.is<Product>());
  SEQUANT_ASSERT(product.as<Product>().factors().size() == 2);
  SEQUANT_ASSERT(sum.is<Sum>());
  SEQUANT_ASSERT(sum.as<Sum>().summands().size() == 3);

  std::cout << sum << std::endl;
  // end-snippet-1

  // start-snippet-2
  // scalar leaves: Constant (a compile-time rational/complex number),
  // Variable (a named runtime scalar), and Power (base^exponent, with base
  // being a Constant or a Variable)
  auto half = ex<Constant>(rational(1, 2));
  auto lambda = ex<Variable>("\\lambda");
  auto lambda_sq = ex<Power>(lambda, 2);

  SEQUANT_ASSERT(lambda_sq.as<Power>().base() == lambda);
  SEQUANT_ASSERT(lambda_sq.as<Power>().exponent() == 2);

  auto scaled = half * g;  // Constant and Tensor combine into one Product
  std::cout << std::format("{}", scaled) << std::endl;
  // end-snippet-2

  // start-snippet-3
  // ResultExpr pairs an expression with an explicit left-hand side, e.g.
  // "R^{a_1}_{i_1} = f^{a_1}_{i_1}"; this is the object handed to optimize()
  // and to the code generators in SeQuant/core/export
  auto result = ResultExpr(Tensor("R", bra{"i_1"}, ket{"a_1"}), f);

  SEQUANT_ASSERT(result.produces_tensor());
  std::cout << std::format("{}", result.expression()) << std::endl;
  // end-snippet-3

  return 0;
}
