//
// Created by Eduard Valeyev on 2019-02-18.
//

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <utility>

namespace sequant {

namespace detail {

std::pair<AbstractTensorPtr, ExprPtr> canonicalize_as_tensor_network(
    const Expr& expr, const CanonicalizeOptions& opts) {
  const auto ctx = get_default_context_snapshot();
  TensorNetwork tn(expr);
  auto byproduct = tn.canonicalize(ctx.cardinal_tensor_labels(), opts);
  SEQUANT_ASSERT(tn.tensors().size() == 1);
  return {tn.tensors()[0], std::move(byproduct)};
}

}  // namespace detail
template <>
const container::svector<std::string>&
NormalOperator<Statistics::FermiDirac>::labels() {
  using namespace std::literals;
  static container::svector<std::string> labels_{"a"s, "ã"s};
  return labels_;
}

template <>
const container::svector<std::string>&
NormalOperator<Statistics::BoseEinstein>::labels() {
  using namespace std::literals;
  static container::svector<std::string> labels_{"b"s, "b̃"s};
  return labels_;
}

}  // namespace sequant
