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

enum class OpType { f, g, t, l, A, L, R };

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
    case OpType::L:
      return L"L";
    case OpType::R:
      return L"R";
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
    case OpType::R:
      return OpClass::ex;
    case OpType::l:
    case OpType::A:
    case OpType::L: return OpClass::deex;
    default:
      throw std::invalid_argument("to_class(OpType op): invalid op");
  }
}

namespace so {

inline constexpr size_t fac(std::size_t n) {
  if (n == 1 || n == 0)
    return 1;
  else
    return n * fac(n - 1);
}

template <std::size_t Nbra, std::size_t Nket, OpType _Op, bool PNO>
struct make_op {
  ExprPtr operator()() const {
    const auto nbra = Nbra;
    const auto nket = Nket;
    const auto pno = PNO;
    OpType op = _Op;
    auto make_idx_vector = [op](size_t n, IndexSpace::Type spacetype) {
      auto space = IndexSpace::instance(spacetype);
      std::vector<Index> result;
      result.reserve(n);
      for (size_t i = 0; i != n; ++i) {
        result.push_back(Index::make_tmp_index(space));
      }
      return result;
    };
    auto make_depidx_vector = [op](size_t n, IndexSpace::Type spacetype, auto&& protoidxs) {
      auto space = IndexSpace::instance(spacetype);
      std::vector<Index> result;
      result.reserve(n);
      for (size_t i = 0; i != n; ++i) {
        result.push_back(Index::make_tmp_index(space, protoidxs, true));
      }
      return result;
    };
    std::vector<Index> braidxs;
    std::vector<Index> ketidxs;
    if (to_class(op) == OpClass::gen) {
      braidxs = make_idx_vector(nbra, IndexSpace::complete);
      ketidxs = make_idx_vector(nket, IndexSpace::complete);
    }
    else {
      auto make_occidxs = [pno,&make_idx_vector](size_t n) {
        return make_idx_vector(n, IndexSpace::active_occupied);
      };
      auto make_uoccidxs = [pno,&make_idx_vector,&make_depidx_vector](size_t n, auto&& occidxs) {
        return pno ? make_depidx_vector(n, IndexSpace::active_unoccupied, occidxs) : make_idx_vector(n, IndexSpace::active_unoccupied);
      };
      if (to_class(op) == OpClass::ex) {
        ketidxs = make_occidxs(nket);
        braidxs = make_uoccidxs(nbra, ketidxs);
      } else {
        braidxs = make_occidxs(nbra);
        ketidxs = make_uoccidxs(nket, braidxs);
      }
    }
    return ex<Constant>(1. / (fac(Nbra) * fac(Nket))) *
    ex<Tensor>(to_wstring(op), braidxs, ketidxs, Symmetry::antisymm) *
    ex<FNOperator>(braidxs, ketidxs, Vacuum::SingleProduct);
  }
};

template <OpType _Op, std::size_t Nbra, std::size_t Nket = Nbra> static const auto Op = make_op<Nbra, Nket, _Op, false>{};

#include "sr_op.impl.hpp"

ExprPtr H1();
ExprPtr H2();
ExprPtr H0mp();
ExprPtr H1mp();
ExprPtr W();
ExprPtr H();

inline ExprPtr vac_av(ExprPtr expr, std::initializer_list<std::pair<int,int>> op_connections = {}, bool use_top = true) {
  FWickTheorem wick{expr};
  wick.full_contractions(true)
      .spinfree(false)
      .use_topology(use_top)
      .set_op_connections(op_connections);
  auto result = wick.compute();
  simplify(result);
  if (Logger::get_instance().wick_stats) {
    std::wcout << "WickTheorem stats: # of contractions attempted = "
               << wick.stats().num_attempted_contractions
               << " # of useful contractions = "
               << wick.stats().num_useful_contractions << std::endl;
  }
  return result;
}

namespace pno {

template <OpType _Op, std::size_t Nbra, std::size_t Nket = Nbra> static const auto Op = make_op<Nbra, Nket, _Op, true>{};

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
