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

enum class IndexSpaceMetric { Unit, General, Invalid };

// clang-format off
/// describes supported symmetries of tensorial objects with respect to permutations of particles (or columns, in tensor notation)
// clang-format on
enum class ParticleSymmetry { symm, nonsymm, invalid };

// clang-format off
/// describes supported symmetries of bra or ket of _particle-symmetric_ tensorial objects
/// @note bra or ket can be symmetric or antisymmetric only if the tensor is particle-symmetric, otherwise it does not make sense to permute indices corresponding to distinguishable particles
// clang-format on
enum class Symmetry { symm, antisymm, nonsymm, invalid };

/// describes supported symmetries of tensorial objects w.r.t. bra-ket exchange
///
/// @note Currently there is no support for swapping bra index with the ket
///       index for a single particle, only whole bra-ket swaps are considered.
enum class BraKetSymmetry { symm, conjugate, nonsymm, invalid };

/// describes type of single-particle basis
enum class SPBasis { spinor, spinfree };

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
    case Symmetry::invalid:
      SEQUANT_ABORT("Invalid symmetry is not allowed here");
  }
  return result;
}

inline std::wstring to_wstring(Symmetry sym) {
  switch (sym) {
    case Symmetry::symm:
      return L"symmetric";
    case Symmetry::antisymm:
      return L"antisymmetric";
    case Symmetry::nonsymm:
      return L"nonsymmetric";
    case Symmetry::invalid:
      return L"invalid";
  }

  SEQUANT_UNREACHABLE;
}

enum class BraKetPos { bra, ket, none };

inline std::wstring to_wolfram(BraKetPos a) {
  switch (a) {
    case BraKetPos::bra:
      return L"indexType[bra]";
    case BraKetPos::ket:
      return L"indexType[ket]";
    case BraKetPos::none:
      return L"indexType[none]";
  }
}

enum class Statistics {
  Null,
  FermiDirac,
  BoseEinstein,
  Arbitrary,
  Invalid = Null
};

enum class Action { create, annihilate, invalid };

/// applies (Hermitian) adjoint to @c action
inline Action adjoint(Action action) {
  return action == Action::create ? Action::annihilate : Action::create;
}

inline std::wstring to_wolfram(Action a) {
  using namespace std::literals;
  return L"indexType["s + (a == Action::create ? L"cre" : L"ann") + L"]";
}

enum class Vacuum { Physical, SingleProduct, MultiProduct, Invalid };

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
