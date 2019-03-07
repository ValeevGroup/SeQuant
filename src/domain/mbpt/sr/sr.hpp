//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_SRCC_HPP
#define SEQUANT_SRCC_HPP

#include "../../../SeQuant/op.hpp"
#include "../../../SeQuant/tensor.hpp"
#include "../../../SeQuant/wick.hpp"

namespace sequant {
namespace mbpt {
namespace sr {

enum class OpType { f, g, t, l, A };

enum class OpClass { ex, deex, gen };

inline std::wstring to_wstring(OpType op) {
  switch (op) {
    case OpType::f:
      return L"f";
    case OpType::g:
      return L"g";
    case OpType::t:
      return L"t";
    case OpType::l:
      return L"λ";
    case OpType::A:
      return L"A";
    default:
      throw std::invalid_argument("to_wstring(OpType op): invalid op");
  }
}

inline OpClass to_class(OpType op) {
  switch (op) {
    case OpType::f:
    case OpType::g:
      return OpClass::gen;
    case OpType::t:
      return OpClass::ex;
    case OpType::l:
    case OpType::A: return OpClass::deex;
    default:
      throw std::invalid_argument("to_class(OpType op): invalid op");
  }
}

namespace so {

inline size_t fac(std::size_t n) {
  if (n == 1 || n == 0)
    return 1;
  else
    return n * fac(n - 1);
}

template <std::size_t N, OpType _Op, bool pno>
struct make_op {
  ExprPtr operator()() const {
    assert(N>0);
    size_t n = N;
    OpType op = _Op;
    auto make_idx_vector = [n,op](IndexSpace::Type spacetype) {
      auto space = IndexSpace::instance(spacetype);
      std::vector<Index> result;
      result.reserve(n);
      for (size_t i = 0; i != N; ++i) {
        result.push_back(Index::make_tmp_index(space));
      }
      return result;
    };
    auto make_depidx_vector = [n,op](IndexSpace::Type spacetype, auto&& protoidxs) {
      auto space = IndexSpace::instance(spacetype);
      std::vector<Index> result;
      result.reserve(n);
      for (size_t i = 0; i != N; ++i) {
        result.push_back(Index::make_tmp_index(space, protoidxs, true));
      }
      return result;
    };
    std::vector<Index> braidxs;
    std::vector<Index> ketidxs;
    if (to_class(op) == OpClass::gen) {
//      braidxs = make_idx_vector(IndexSpace::complete);
//      ketidxs = make_idx_vector(IndexSpace::complete);
      braidxs = make_idx_vector(IndexSpace::all);
      ketidxs = make_idx_vector(IndexSpace::all);
    }
    else {
      auto occidxs = make_idx_vector(IndexSpace::active_occupied);
      auto uoccidxs = pno ? make_depidx_vector(IndexSpace::active_unoccupied, occidxs) : make_idx_vector(IndexSpace::active_unoccupied);
      if (to_class(op) == OpClass::ex) {
        braidxs = std::move(uoccidxs);
        ketidxs = std::move(occidxs);
      } else {
        ketidxs = std::move(uoccidxs);
        braidxs = std::move(occidxs);
      }
    }
    return ex<Constant>(1. / (fac(N) * fac(N))) *
    ex<Tensor>(to_wstring(op), braidxs, ketidxs, Symmetry::antisymm) *
    ex<FNOperator>(braidxs, ketidxs, Vacuum::SingleProduct);
  }
};

template <std::size_t N, OpType _Op> static const auto Op = make_op<N, _Op, false>{};

#include "sr_op.impl.hpp"

ExprPtr H1();
ExprPtr H2();
ExprPtr H0mp();
ExprPtr H1mp();
ExprPtr W();
ExprPtr H();

inline ExprPtr vac_av(ExprPtr expr, std::initializer_list<std::pair<int,int>> op_connections = {}) {
//  std::wcout << "vac_av: input = " << to_latex_align(expr, 20, 1) << std::endl;
  auto wick = FWickTheorem{expr}; wick.full_contractions(true).spinfree(false).use_topology(true);
  wick.set_op_connections(op_connections);
  auto result = wick.compute();
  simplify(result);
  std::wcout << "WickTheorem stats: # of contractions attempted = "
             << wick.stats().num_attempted_contractions
             << " # of useful contractions = "
             << wick.stats().num_useful_contractions << std::endl;
  return result;
}

namespace pno {

template <std::size_t N, OpType _Op> static const auto Op = make_op<N, _Op, true>{};

#include "sr_op.impl.hpp"

using sequant::mbpt::sr::so::H;
using sequant::mbpt::sr::so::H0mp;
using sequant::mbpt::sr::so::H1;
using sequant::mbpt::sr::so::H1mp;
using sequant::mbpt::sr::so::H2;
using sequant::mbpt::sr::so::vac_av;

}  // namespace pno

}  // namespace so
}  // namespace sr
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_SRCC_HPP
