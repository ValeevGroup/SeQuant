#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/format.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <format>
#include <vector>

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
  // TensorPower is the element-wise rational power of a tensor leaf; pow()
  // builds it, e.g. the inverse of a symmetric denominator tensor
  auto D = ex<Tensor>(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                      TensorSymmetries{.perm = Symmetry::Symm,
                                       .braket = BraKetSymmetry::Symm,
                                       .column = ColumnSymmetry::Symm});
  auto D_inv = pow(D, rational{-1});

  SEQUANT_ASSERT(D_inv.is<TensorPower>());
  SEQUANT_ASSERT(D_inv.as<AbstractTensor>()._symmetry() == Symmetry::Symm);

  // an even power of an antisymmetric tensor is symmetric
  auto t = ex<Tensor>(L"t", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                      Symmetry::Antisymm);
  SEQUANT_ASSERT(pow(t, rational{2})->as<AbstractTensor>()._symmetry() ==
                 Symmetry::Symm);

  std::wcout << std::format(L"{}", D_inv) << std::endl;
  // end-snippet-4

  // start-snippet-5
  // Compose through ResultExpr heads: each stage is its own equation whose
  // head declares the external indices and symmetries once, and later stages
  // reuse the head via result_as_tensor() instead of respelling it
  const TensorSymmetries antisymm{.perm = Symmetry::Antisymm};
  auto M = ex<Tensor>(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, antisymm);
  auto L = ex<Tensor>(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}, antisymm);

  // ℓ = M ⊙ D^{-1}: a Hadamard (element-wise) product, legal because the
  // head names its externals; the head declares ℓ's antisymmetry
  auto ell = ResultExpr(
      Tensor(L"ℓ", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, antisymm),
      M * pow(D, rational{-1}));

  // E = L · ℓ: an ordinary contraction of the leaf ℓ, taken from the head
  auto E = ResultExpr(Variable(L"E"), L * ex<Tensor>(ell.result_as_tensor()));

  SEQUANT_ASSERT(!E.produces_tensor());
  std::wcout << std::format(L"{}", ell.expression()) << std::endl;
  std::wcout << std::format(L"{}", E.expression()) << std::endl;
  // end-snippet-5

  // start-snippet-6
  // under a head, a summand may omit head indices and is then constant along
  // those modes, e.g. the orbital-energy denominator
  // D{i_1,i_2;a_1,a_2} = ε{;a_1} + ε{;a_2} - ε{i_1;} - ε{i_2;}
  auto denominator =
      ResultExpr(Tensor(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}),
                 ex<Tensor>(L"ε", bra{}, ket{L"a_1"}) +
                     ex<Tensor>(L"ε", bra{}, ket{L"a_2"}) -
                     ex<Tensor>(L"ε", bra{L"i_1"}, ket{}) -
                     ex<Tensor>(L"ε", bra{L"i_2"}, ket{}));
  SEQUANT_ASSERT(is_valid(denominator));

  // binarize() and optimize() lower it first: each summand is multiplied by
  // one reserved all-ones leaf 𝟙 per head index it lacks
  auto lowered = lower_broadcast_sum(denominator);
  [[maybe_unused]] const auto& first_term =
      lowered.expression()->as<Sum>().summand(0);
  SEQUANT_ASSERT(first_term->as<Product>().factors().size() == 4);

  std::wcout << std::format(L"{}", lowered.expression()) << std::endl;
  // end-snippet-6

  // start-snippet-7
  // binarize() of a system inlines each definition into its uses: a leaf
  // names a head if label, ranks and index spaces agree (a leaf that names a
  // head with other symmetries throws); D's head is D, the base of the power
  // in ℓ, defined as the broadcast sum above
  auto D_def = ResultExpr(D->as<Tensor>(), denominator.expression());

  // R{i_2;i_3} = L{a_1,a_2;i_1,i_3} · ℓ{i_1,i_2;a_1,a_2}, ℓ taken from its head
  auto R = ResultExpr(
      Tensor(L"R", bra{L"i_2"}, ket{L"i_3"}),
      ex<Tensor>(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_3"}, antisymm) *
          ex<Tensor>(ell.result_as_tensor()));

  // one tree per equation, in input order; D and ℓ are inlined, R is a root
  [[maybe_unused]] auto sys = binarize(std::vector{D_def, ell, R});
  SEQUANT_ASSERT((sys.consumed == container::svector<bool>{true, true, false}));
  SEQUANT_ASSERT(sys.roots().size() == 1);
  // end-snippet-7

  return 0;
}
