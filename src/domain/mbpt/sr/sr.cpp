//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "sr.hpp"

namespace sequant {
namespace mbpt {
namespace sr {
namespace so {

ExprPtr H1() {
  return Op<OpType::f, 1>();
}
ExprPtr H2() {
  return Op<OpType::g, 2>();
}
ExprPtr H0mp() {
  return H1();
}
ExprPtr H1mp() {
  return H2();
}
ExprPtr H() {
  return H1() + H2();
}
ExprPtr W() {
  return H2();
}

}  // namespace so
}  // namespace sr
}  // namespace mbpt
}  // namespace sequant