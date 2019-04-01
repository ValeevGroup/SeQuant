//
// Created by Eduard Valeyev on 2019-02-13.
//

#ifndef SEQUANT_ATTR_HPP
#define SEQUANT_ATTR_HPP

#include <string>

namespace sequant {

enum class IndexSpaceMetric {
  Unit,
  General,
  Invalid
};

/// describes supported symmetries of bra or ket of tensorial objects
enum class Symmetry { symm, antisymm, nonsymm, invalid };

/// describes supported symmetries of tensorial objects w.r.t. bra-ket exchange
enum class BraKetSymmetry { symm, conjugate, nonsymm, invalid };

inline std::wstring to_wolfram(const Symmetry& symmetry) {
  std::wstring result;
  switch (symmetry) {
    case Symmetry::symm:
      result = L"indexSymm[1]";
      break;
    case Symmetry::antisymm:
      result = L"indexSymm[-1]";
      break;
    case Symmetry::nonsymm:
      result = L"indexSymm[0]";
      break;
    default:
      abort();
  }
  return result;
}

enum class BraKetPos { bra, ket };

inline std::wstring to_wolfram(BraKetPos a) {
  using namespace std::literals;
  return L"indexType["s + (a == BraKetPos::bra ? L"bra" : L"ket") + L"]";
}

enum class Statistics { BoseEinstein, FermiDirac };

enum class Action { create, annihilate, invalid };

/// applies (Hermitian) adjoint to @c action
inline Action adjoint(Action action) { return action == Action::create ? Action::annihilate : Action::create; }

inline std::wstring to_wolfram(Action a) {
  using namespace std::literals;
  return L"indexType["s + (a == Action::create ? L"cre" : L"ann") + L"]";
}

enum class Vacuum {
  Physical,
  SingleProduct,
  MultiProduct,
  Invalid
};

inline std::wstring to_string(Vacuum V) {
  switch (V) {
    case Vacuum::Physical:
      return L"PhysicalVacuum";
    case Vacuum::SingleProduct:
      return L"SingleProductVacuum";
    case Vacuum::MultiProduct:
      return L"MultiProductVacuum";
    case Vacuum::Invalid:
      return L"InvalidVacuum";
    default:
      abort();
  }
}

}  // namespace sequant

#endif //SEQUANT_ATTR_HPP
