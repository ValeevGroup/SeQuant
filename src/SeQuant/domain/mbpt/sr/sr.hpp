//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_SRCC_HPP
#define SEQUANT_SRCC_HPP

#include <limits>
#include <utility>

#include "SeQuant/domain/mbpt/op.hpp"
#include "SeQuant/core/expr_fwd.hpp"

namespace sequant {
namespace mbpt {
namespace sr {
namespace so {

class make_op {
  std::size_t nbra_, nket_;
  OpType op_;
  bool csv_;

 public:
  make_op(std::size_t nbra, std::size_t nket, OpType op, bool csv);

  ExprPtr operator()(bool complete_unoccupieds = false) const;
};

make_op Op(OpType _Op, std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max());

#include "sr_op.impl.hpp"

ExprPtr H1();
ExprPtr H2();
ExprPtr H0mp();
ExprPtr H1mp();
ExprPtr W();
ExprPtr H();

ExprPtr vac_av(ExprPtr expr, std::initializer_list<std::pair<int,int>> op_connections = {}, bool use_top = true);

namespace csv {

make_op Op(OpType _Op, std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max());

#include "sr_op.impl.hpp"

using sequant::mbpt::sr::so::H;
using sequant::mbpt::sr::so::H0mp;
using sequant::mbpt::sr::so::H1;
using sequant::mbpt::sr::so::H1mp;
using sequant::mbpt::sr::so::H2;
using sequant::mbpt::sr::so::vac_av;

}  // namespace csv

}  // namespace so
}  // namespace sr
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_SRCC_HPP
