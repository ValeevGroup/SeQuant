#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/format.hpp>
#include <SeQuant/core/utility/macros.hpp>

int main() {
  using namespace sequant;

  // start-snippet-1
  // two Products that are mathematically identical, differing only by a
  // relabeling of their (implicitly summed) dummy indices
  auto term_a = ex<Tensor>("g", bra{"i_1", "i_2"}, ket{"a_1", "a_2"}) *
                ex<Tensor>("t", bra{"a_1", "a_2"}, ket{"i_1", "i_2"});
  auto term_b = ex<Tensor>("g", bra{"i_2", "i_1"}, ket{"a_2", "a_1"}) *
                ex<Tensor>("t", bra{"a_2", "a_1"}, ket{"i_2", "i_1"});

  // before canonicalization the two Products are not recognized as equal ...
  SEQUANT_ASSERT(!(term_a == term_b));

  // ... but canonicalize() puts both into the same normal form, using the
  // tensor network machinery (which relies on the bundled bliss graph
  // automorphism library) to consistently relabel dummy indices
  term_a->canonicalize();
  term_b->canonicalize();
  SEQUANT_ASSERT(term_a == term_b);
  // end-snippet-1

  // start-snippet-2
  // simplify() combines canonicalization with cheap algebraic clean-up, so a
  // Sum of the two forms above collapses to a single term with a factor of 2
  auto sum = ex<Tensor>("g", bra{"i_1", "i_2"}, ket{"a_1", "a_2"}) *
                 ex<Tensor>("t", bra{"a_1", "a_2"}, ket{"i_1", "i_2"}) +
             ex<Tensor>("g", bra{"i_2", "i_1"}, ket{"a_2", "a_1"}) *
                 ex<Tensor>("t", bra{"a_2", "a_1"}, ket{"i_2", "i_1"});
  simplify(sum);

  SEQUANT_ASSERT(sum.is<Product>());
  std::cout << sum << std::endl;
  // end-snippet-2

  return 0;
}
