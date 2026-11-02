#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/format.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/polynomial_division.hpp>

#include <format>

int main() {
  using namespace sequant;

  // start-snippet-1
  // Tensor is a leaf Expr: a labeled tensor with a bra and a ket index set
  auto f = ex<Tensor>(L"f", bra{L"i_1"}, ket{L"a_1"});
  auto g = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});

  // ExprPtr arithmetic builds Product and Sum nodes; both are associative and
  // flatten automatically, so nesting them never grows the tree depth
  auto product = f * g;      // a Product of 2 factors
  auto sum = f * g + g + f;  // a Sum of 3 summands, not a Sum of a Sum and f

  SEQUANT_ASSERT(product.is<Product>());
  SEQUANT_ASSERT(product.as<Product>().factors().size() == 2);
  SEQUANT_ASSERT(sum.is<Sum>());
  SEQUANT_ASSERT(sum.as<Sum>().summands().size() == 3);

  std::wcout << sum << std::endl;
  // end-snippet-1

  // start-snippet-2
  // scalar leaves: Constant (a compile-time rational/complex number),
  // Variable (a named runtime scalar), and Power (base^exponent, with base
  // being a Constant or a Variable)
  auto half = ex<Constant>(rational(1, 2));
  auto lambda = ex<Variable>(L"\\lambda");
  auto lambda_sq = ex<Power>(lambda, 2);

  SEQUANT_ASSERT(lambda_sq.as<Power>().base() == lambda);
  SEQUANT_ASSERT(lambda_sq.as<Power>().exponent() == 2);

  auto scaled = half * g;  // Constant and Tensor combine into one Product
  std::wcout << std::format(L"{}", scaled) << std::endl;
  // end-snippet-2

  // start-snippet-3
  // ResultExpr pairs an expression with an explicit left-hand side, e.g.
  // "R^{a_1}_{i_1} = f^{a_1}_{i_1}"; this is the object handed to optimize()
  // and to the code generators in SeQuant/core/export
  auto result = ResultExpr(Tensor(L"R", bra{L"i_1"}, ket{L"a_1"}), f);

  SEQUANT_ASSERT(result.produces_tensor());
  std::wcout << std::format(L"{}", result.expression()) << std::endl;
  // end-snippet-3

  // start-snippet-4
  // the intermediate I^{a_1}_{i_1} = f^{a_1}_{a_2} t^{a_2}_{i_1} +
  // g^{a_1}_{i_1}
  auto intermediate =
      ResultExpr(Tensor(L"I", bra{L"a_1"}, ket{L"i_1"}),
                 ex<Tensor>(L"f", bra{L"a_1"}, ket{L"a_2"}) *
                         ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_1"}) +
                     ex<Tensor>(L"g", bra{L"a_1"}, ket{L"i_1"}));

  // an expression containing the intermediate in expanded form (using
  // different dummy indices) plus an unrelated term
  auto h = ex<Tensor>(L"h", bra{L"i_3"}, ket{L"a_3"});
  auto expr = ex<Tensor>(L"f", bra{L"a_3"}, ket{L"a_4"}) *
                  ex<Tensor>(L"t", bra{L"a_4"}, ket{L"i_3"}) * h +
              ex<Tensor>(L"g", bra{L"a_3"}, ket{L"i_3"}) * h +
              ex<Variable>(L"E");

  // polynomial division: expr = divisor * quotient + remainder
  DivisionResult division = divide(expr, intermediate.expression(),
                                   DivisionMode::Split, IndexMatching::Exact);
  SEQUANT_ASSERT(division.groups.size() == 1);
  SEQUANT_ASSERT(!division.exact());  // the remainder is E
  // (f t + g) h + E
  std::wcout << std::format(L"{}", division.to_expr()) << std::endl;

  // factor out the intermediate and replace it by its result tensor
  // I^{a_1}_{i_1} h^{i_1}_{a_1} + E
  auto factored =
      factor_out(expr, intermediate, DivisionMode::Split, IndexMatching::Exact);
  std::wcout << std::format(L"{}", factored) << std::endl;
  // end-snippet-4

  return 0;
}
