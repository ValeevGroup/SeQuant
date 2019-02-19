//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "sr.hpp"

namespace sequant2 {
namespace mbpt {
namespace sr {
namespace so {

ExprPtr H1() {
  return Op<1,OpType::f>();
}
ExprPtr H2() {
  return Op<2,OpType::g>();
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
}  // namespace sequant2