//
// Created by Eduard Valeyev on 2019-02-18.
//

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>

namespace sequant {
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
