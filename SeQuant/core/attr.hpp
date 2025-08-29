//
// Created by Eduard Valeyev on 2019-02-13.
//

#ifndef SEQUANT_ATTR_HPP
#define SEQUANT_ATTR_HPP

#include <SeQuant/core/utility/macros.hpp>

#include <cassert>
#include <cstdlib>
#include <string>

namespace sequant {

enum class IndexSpaceMetric { Unit, General };

// clang-format off
/// describes supported symmetries of tensorial objects with respect to permutations of columns (in tensor notation)
// clang-format on
enum class ColumnSymmetry { Symm, Nonsymm };

// clang-format off
/// describes supported symmetries of bra or ket of _particle-symmetric_ tensorial objects
/// @note bra or ket can be symmetric or antisymmetric only if the tensor is particle-symmetric, otherwise it does not make sense to permute indices corresponding to distinguishable particles
// clang-format on
enum class Symmetry { Symm, Antisymm, Nonsymm };

/// describes supported symmetries of tensorial objects w.r.t. bra-ket exchange
///
/// @note Currently there is no support for swapping bra index with the ket
///       index for a single particle, only whole bra-ket swaps are considered.
enum class BraKetSymmetry { Symm, Conjugate, Nonsymm };

/// describes type of single-particle basis
enum class SPBasis { Spinor, Spinfree };

inline std::wstring to_wolfram(const Symmetry& symmetry) {
  std::wstring result;
  switch (symmetry) {
    case Symmetry::Symm:
      result = L"indexSymm[1]";
      break;
    case Symmetry::Antisymm:
      result = L"indexSymm[-1]";
      break;
    case Symmetry::Nonsymm:
      result = L"indexSymm[0]";
      break;
  }
  return result;
}

inline std::wstring to_wstring(Symmetry sym) {
  switch (sym) {
    case Symmetry::Symm:
      return L"symmetric";
    case Symmetry::Antisymm:
      return L"antisymmetric";
    case Symmetry::Nonsymm:
      return L"nonsymmetric";
  }

  SEQUANT_UNREACHABLE;
}

enum class BraKetPos {
  Bra,
  Ket,
};

inline std::wstring to_wolfram(BraKetPos a) {
  switch (a) {
    case BraKetPos::Bra:
      return L"indexType[bra]";
    case BraKetPos::Ket:
      return L"indexType[ket]";
  }
}

enum class Statistics {
  FermiDirac,
  BoseEinstein,
  Arbitrary,
};

enum class Action { Create, Annihilate };

/// applies (Hermitian) adjoint to @c action
inline Action adjoint(Action action) {
  return action == Action::Create ? Action::Annihilate : Action::Create;
}

inline std::wstring to_wolfram(Action a) {
  using namespace std::literals;
  return L"indexType["s + (a == Action::Create ? L"cre" : L"ann") + L"]";
}

enum class Vacuum { Physical, SingleProduct, MultiProduct };

inline std::wstring to_string(Vacuum V) {
  switch (V) {
    case Vacuum::Physical:
      return L"PhysicalVacuum";
    case Vacuum::SingleProduct:
      return L"SingleProductVacuum";
    case Vacuum::MultiProduct:
      return L"MultiProductVacuum";
  }

  SEQUANT_UNREACHABLE;
}

/// describes LaTeX typesetting convention for contravariant (bra, annihilation)
/// and covariant (ket, creation) indices
enum class BraKetTypesetting {
  /// contravariants as subscripts
  ContraSub,
  CoSuper = ContraSub,
  BraSub = ContraSub,
  KetSuper = ContraSub,
  /// covariant as subscripts
  CoSub,
  ContraSuper = CoSub,
  BraSuper = CoSub,
  KetSub = CoSub,
};

/// describes typesetting convention for bra/ket tensor slots
enum class BraKetSlotTypesetting {
  /// all indices are typeset naively, with empty slots indicated by
  /// `\testvisiblespace`
  Naive,
  /// tensor package is used to ensure consistent alignment of bra and ket slots
  TensorPackage
};

}  // namespace sequant

#endif  // SEQUANT_ATTR_HPP
