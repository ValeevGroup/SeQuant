//
// Created by Eduard Valeyev on 2019-02-13.
//

#ifndef SEQUANT2_ATTR_HPP
#define SEQUANT2_ATTR_HPP

#include <string>

namespace sequant2 {

enum class IndexSpaceMetric {
  Unit,
  General,
  Invalid
};

enum class Symmetry { symm, antisymm, nonsymm };

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
  }
  return result;
}

enum class BraKetPos { bra, ket };

inline std::wstring to_wolfram(BraKetPos a) {
  using namespace std::literals;
  return L"indexType["s + (a == BraKetPos::bra ? L"bra" : L"ket") + L"]";
}

enum class Statistics { BoseEinstein, FermiDirac };

enum class Action { create, annihilate };

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

}  // namespace sequant2

#endif //SEQUANT2_ATTR_HPP
