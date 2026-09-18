#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>

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
  assert(!(term_a == term_b));

  // ... but canonicalize() puts both into the same normal form, using the
  // tensor network machinery (which relies on the bundled bliss graph
  // automorphism library) to consistently relabel dummy indices
  term_a->canonicalize();
  term_b->canonicalize();
  assert(term_a == term_b);
  // end-snippet-1

  // start-snippet-2
  // simplify() combines canonicalization with cheap algebraic clean-up, so a
  // Sum of the two forms above collapses to a single term with a factor of 2
  auto sum = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}) *
                 ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) +
             ex<Tensor>(L"g", bra{L"i_2", L"i_1"}, ket{L"a_2", L"a_1"}) *
                 ex<Tensor>(L"t", bra{L"a_2", L"a_1"}, ket{L"i_2", L"i_1"});
  simplify(sum);

  assert(sum.is<Product>());
  std::wcout << to_latex(sum) << std::endl;
  // end-snippet-2

  return 0;
}
