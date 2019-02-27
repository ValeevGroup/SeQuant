//
// Created by Eduard Valeyev on 2019-02-18.
//

#include "op.hpp"

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
}  // namespace sequant

