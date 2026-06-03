//
// Created by Eduard Valeyev on 2019-02-13.
//

#ifndef SEQUANT_ATTR_HPP
#define SEQUANT_ATTR_HPP

#include <SeQuant/core/utility/macros.hpp>

#include <cassert>
#include <cstdlib>
#include <ostream>
#include <string>

namespace sequant {

enum class IndexSpaceMetric { Unit, General };

/// describes the scalar field over which the vector spaces representing the
/// bra/ket modes are defined. Physically this is a single global choice for a
/// computation, but it is carried as per-IndexSpace metadata
/// (IndexSpace::field()); a tensor's effective field is resolved from its
/// bra/ket spaces (see sequant::base_field), so a global real/complex switch is
/// expressed by setting the field uniformly on all spaces. It controls how the
/// bra<->ket (Riesz) dual pairing is realized -- linearly via a symmetric
/// bilinear form (`Field::Real`, the adjoint is the transpose) or antilinearly
/// via a Hermitian sesquilinear form (`Field::Complex`, the adjoint is the
/// conjugate-transpose). Together with a tensor's #Hermiticity it determines
/// the tensor's #BraKetSymmetry.
/// @sa Hermiticity, BraKetSymmetry, to_braket_symmetry, IndexSpace::field,
///     sequant::base_field
enum class Field { Real, Complex };

// clang-format off
/// describes supported symmetries of tensorial objects with respect to permutations of columns (in tensor notation), i.e., pairs of {bra[i],ket[i]} slots
// clang-format on
enum class ColumnSymmetry { Symm, Nonsymm };

// clang-format off
/// describes supported symmetries of bra or ket of _particle-symmetric_ tensorial objects
/// @note bra or ket can be symmetric or antisymmetric only if the tensor is particle-symmetric, otherwise it does not make sense to permute indices corresponding to distinguishable particles
// clang-format on
enum class Symmetry { Symm, Antisymm, Nonsymm };

/// describes supported symmetries of tensorial objects w.r.t. bra-ket exchange,
/// i.e., the swap of bra slot bundle with ket slot bundle
///
/// @note Currently there is no support for swapping bra index with the ket
///       index for a single particle, only whole bra-ket swaps are considered.
enum class BraKetSymmetry { Symm, Conjugate, Nonsymm };

/// describes the abstract symmetry of a tensorial object under (Hermitian)
/// adjoint, i.e. whether the abstract tensor equals (`Hermitian`), equals minus
/// (`AntiHermitian`), or is unrelated to (`NonHermitian`) its own adjoint.
///
/// Unlike #BraKetSymmetry this is a field-agnostic property of the abstract
/// tensor: e.g. a 2-electron integral is `Hermitian` whether the computation is
/// real or complex, a cluster amplitude is `NonHermitian` in either. The
/// observable bra<->ket exchange symmetry (#BraKetSymmetry) is the *derived*
/// composition of this trait with the ambient #Field (see to_braket_symmetry):
/// `Hermitian` becomes `Symm` over a real field and `Conjugate` over a complex
/// one; the remaining cases become `Nonsymm`.
/// @note `AntiHermitian` cannot yet be represented in #BraKetSymmetry (which
///       lacks the antisymmetric cases) and therefore currently derives to
///       `Nonsymm`; the trait is recorded for forward compatibility.
enum class Hermiticity { Hermitian, AntiHermitian, NonHermitian };

/// @return the bra<->ket exchange symmetry implied by a tensor's @p hermiticity
///         (a field-agnostic abstract trait) over the scalar @p field
/// @sa Hermiticity, Field, BraKetSymmetry
inline BraKetSymmetry to_braket_symmetry(Hermiticity hermiticity, Field field) {
  switch (hermiticity) {
    case Hermiticity::Hermitian:
      return field == Field::Real ? BraKetSymmetry::Symm
                                  : BraKetSymmetry::Conjugate;
    case Hermiticity::AntiHermitian:
      // BraKetSymmetry has no antisymmetric/anti-conjugate cases yet
      return BraKetSymmetry::Nonsymm;
    case Hermiticity::NonHermitian:
      return BraKetSymmetry::Nonsymm;
  }
  SEQUANT_UNREACHABLE;
}

/// @return the (field-agnostic) #Hermiticity consistent with @p
/// braket_symmetry;
///         used to back-fill the trait when a tensor is constructed by
///         specifying its #BraKetSymmetry directly (legacy/explicit path).
/// @note `Symm` and `Conjugate` both map to `Hermitian` (they differ only in
///       the field, which is not encoded here); `Nonsymm` maps to
///       `NonHermitian` (the `AntiHermitian` preimage is not recoverable).
inline Hermiticity to_hermiticity(BraKetSymmetry braket_symmetry) {
  switch (braket_symmetry) {
    case BraKetSymmetry::Symm:
    case BraKetSymmetry::Conjugate:
      return Hermiticity::Hermitian;
    case BraKetSymmetry::Nonsymm:
      return Hermiticity::NonHermitian;
  }
  SEQUANT_UNREACHABLE;
}

/// describes whether to SEQUANT_ASSERT the vector space semantics of bra/ket
/// slots in tensor networks
/// @sa Context::assert_strict_braket_symmetry
enum class AssertStrictBraKetSymmetry { Yes, No };

/// describes type of single-particle basis
enum class SPBasis { Spinor, Spinfree };

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

/// index slot types
///
/// @note This does not include slot bundles, like braket, etc.
enum class SlotType {
  Bra,
  Ket,
  Aux,
  Proto,
};

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits>& operator<<(
    std::basic_ostream<CharT, Traits>& stream, SlotType origin) {
  switch (origin) {
    case SlotType::Bra:
      stream << "Bra";
      break;
    case SlotType::Ket:
      stream << "Ket";
      break;
    case SlotType::Aux:
      stream << "Aux";
      break;
    case SlotType::Proto:
      stream << "Proto";
      break;
  }
  return stream;
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
