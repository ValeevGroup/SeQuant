//
// Created by Eduard Valeyev on 2019-02-18.
//

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>

#include <limits>

namespace sequant {
namespace detail {
OpIdRegistrar::OpIdRegistrar() {
  auto id = std::numeric_limits<Expr::type_id_type>::max();
  Expr::set_type_id<FNOperator>(id);
  Expr::set_type_id<BNOperator>(--id);
  Expr::set_type_id<FOperator>(--id);
  Expr::set_type_id<BOperator>(--id);
}
}  // namespace detail

template <>
const container::svector<std::wstring>&
NormalOperator<Statistics::FermiDirac>::labels() {
  using namespace std::literals;
  static container::svector<std::wstring> labels_{L"a"s, L"ã"s};
  return labels_;
}

template <>
const container::svector<std::wstring>&
NormalOperator<Statistics::BoseEinstein>::labels() {
  using namespace std::literals;
  static container::svector<std::wstring> labels_{L"b"s, L"b̃"s};
  return labels_;
}

}  // namespace sequant
