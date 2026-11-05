#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/format.hpp>
#include <SeQuant/core/utility/macros.hpp>

int main() {
  using namespace sequant;

  // start-snippet-1
  // two Products that are mathematically identical, differing only by a
  // relabeling of their (implicitly summed) dummy indices
  auto term_a = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"});
  auto term_b = ex<Tensor>(L"g", bra{L"i_2", L"i_1"}, ket{L"a_2", L"a_1"}) *
                ex<Tensor>(L"t", bra{L"a_2", L"a_1"}, ket{L"i_2", L"i_1"});

  // before canonicalization the two Products are not recognized as equal ...
  SEQUANT_ASSERT(!(term_a == term_b));

  // ... but canonicalize() puts both into the same normal form, using the
  // tensor network machinery (which relies on the bundled bliss graph
  // automorphism library) to consistently relabel dummy indices
  term_a->canonicalize(deduce_named_indices);
  term_b->canonicalize(deduce_named_indices);
  SEQUANT_ASSERT(term_a == term_b);
  // end-snippet-1

  // start-snippet-2
  // simplify() combines canonicalization with cheap algebraic clean-up, so a
  // Sum of the two forms above collapses to a single term with a factor of 2
  auto sum = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}) *
                 ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) +
             ex<Tensor>(L"g", bra{L"i_2", L"i_1"}, ket{L"a_2", L"a_1"}) *
                 ex<Tensor>(L"t", bra{L"a_2", L"a_1"}, ket{L"i_2", L"i_1"});
  simplify(sum, deduce_named_indices);

  SEQUANT_ASSERT(sum.is<Product>());
  std::wcout << sum << std::endl;
  // end-snippet-2

  // start-snippet-3
  // the named (external) indices are never deduced silently: a ResultExpr
  // carries them as the indices of its result ...
  ResultExpr residual(Tensor(L"R", bra{L"a_1"}, ket{L"i_1"}),
                      ex<Tensor>(L"f", bra{L"a_1"}, ket{L"a_2"}) *
                              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_1"}) +
                          ex<Tensor>(L"f", bra{L"a_1"}, ket{L"a_3"}) *
                              ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_1"}));
  simplify(residual);
  // the two terms differ only in the dummy index, so they combine
  SEQUANT_ASSERT(residual.expression().is<Product>());
  // ... while for a plain ExprPtr the caller opts in to deducing them as the
  // indices that occur once in each product (wrong for e.g. tensor
  // hypercontraction, in which an external index may occur several times)
  auto e = ex<Tensor>(L"f", bra{L"a_1"}, ket{L"a_2"}) *
           ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_1"});
  simplify(e, deduce_named_indices);
  // end-snippet-3

  // start-snippet-4
  // canonicalization produces one fixed form, meant to identify equivalent
  // expressions; prettify() returns a copy in a form meant for reading,
  // here with the tensors labeled t first
  const auto pretty = prettify(
      residual,
      PrettifyOptions{
          .presentation = {.cardinal_tensor_labels =
                               container::vector<std::wstring>{L"t"}}});
  std::wcout << to_latex(pretty) << std::endl;
  // end-snippet-4

  return 0;
}
