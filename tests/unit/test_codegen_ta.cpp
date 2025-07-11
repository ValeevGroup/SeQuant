#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/domain/eval/eval_fwd.hpp>
#include <format>

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

namespace sequant {
constexpr std::string_view codegen_label(EvalOp op) noexcept {
  return op == EvalOp::Product ? "*" : "+";
}

std::string codegen_label(meta::can_evaluate auto const& node) {
  if (node->is_scalar()) return node->label();
  return std::format("{}(\"{}\")", to_string(node->as_tensor().label()),
                     node->annot());
}

std::string codegen(meta::can_evaluate auto const& node) {
  if (node.leaf()) return codegen_label(node);
  return std::format("({} {} {})", codegen(node.left()),
                     codegen_label(node->op_type().value()),
                     codegen(node.right()));
}

}  // namespace sequant

using namespace sequant;
namespace vws = ranges::views;
namespace rng = ranges;

TEST_CASE("TA code generation", "[code_gen]") {
  constexpr std::wstring_view expr =
      L"1/2 g{a_1,a_2;i_1,i_2}:N-C-S - 1 f{i_3;i_2}:N-C-S * "
      L"t{a_1,a_2;i_1,i_3}:N-C-S + 1/2 g{a_1,a_2;a_3,a_4}:N-C-S * "
      L"t{a_3,a_4;i_1,i_2}:N-C-S + 2 g{i_3,a_1;a_3,i_1}:N-C-S * "
      L"t{a_2,a_3;i_2,i_3}:N-C-S + 1/2 g{i_3,i_4;i_1,i_2}:N-C-S * "
      L"t{a_1,a_2;i_3,i_4}:N-C-S - 1 g{i_3,a_1;i_2,i_1}:N-C-S * "
      L"t{a_2;i_3}:N-C-S - 1 g{i_3,a_1;i_1,a_3}:N-C-S * "
      L"t{a_2,a_3;i_2,i_3}:N-C-S + f{a_2;a_3}:N-C-S * t{a_1,a_3;i_1,i_2}:N-C-S "
      L"- 1 g{i_3,a_1;a_3,i_1}:N-C-S * t{a_2,a_3;i_3,i_2}:N-C-S - 1 "
      L"g{i_3,a_2;i_1,a_3}:N-C-S * t{a_1,a_3;i_3,i_2}:N-C-S";
  auto nodes =
      parse_expr(expr)->as<Sum>() | vws::transform(binarize<EvalExprTA>);
  for (auto const& n : nodes) std::cout << codegen(n) << std::endl;
}
