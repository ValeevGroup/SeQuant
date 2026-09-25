#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>

#include <iostream>

int main() {
  using namespace sequant;

  // start-snippet-1
  // a small tensor equation: R^{a_1}_{i_1} = f^{a_1}_{i_1} + t^{a_1}_{i_1}
  auto result = ResultExpr(Tensor(L"R", bra{L"i_1"}, ket{L"a_1"}),
                           ex<Tensor>(L"f", bra{L"i_1"}, ket{L"a_1"}) +
                               ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));

  // export_expression() walks a binarized evaluation tree and invokes the
  // matching Generator callback (create/load/compute/unload/...) for every
  // step; TextGenerator is a dependency-free backend that renders those
  // callbacks as human-readable pseudocode -- useful for prototyping a new
  // Generator subclass, or for debugging what a real backend would emit
  TextGenerator<TextGeneratorContext> generator;
  TextGeneratorContext ctx;
  export_expression(to_export_tree(result), generator, ctx);

  std::cout << generator.get_generated_code() << std::endl;
  // end-snippet-1

  return 0;
}
