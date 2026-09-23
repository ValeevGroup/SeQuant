//
// Created by Eduard Valeyev on 2019-02-13.
//

#ifndef SEQUANT_ATTR_HPP
#define SEQUANT_ATTR_HPP

#include <SeQuant/core/utility/macros.hpp>

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <ostream>
#include <string>
#include <type_traits>

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

/// behaviour of the operator a c-number tensor represents under complex
/// conjugation K in the position representation: `K I K⁻¹ = +I` (`Even`,
/// a "real operator" in Wigner's sense), `−I` (`Odd`, e.g. the momentum
/// −i∇; ∇ itself is `Even` and anti-Hermitian), or
/// neither (`None`). Basis-agnostic, like #Hermiticity: the elementwise
/// behaviour of the array is derived with the basis #Field, see
/// to_conjugation_symmetry(). Over a complex basis conjugation maps the array
/// to the array in the conjugated basis, which is not an elementwise
/// relation, so only a real basis exposes the parity.
enum class ConjugationParity { Even, Odd, None };

/// relation between the array read with bra and ket exchanged and the array
/// as written: `T{q;p} = +T{p;q}` (`Symm`), `−T{p;q}` (`Antisymm`),
/// `conj(T{p;q})` (`Conjugate`), `−conj(T{p;q})` (`AntiConjugate`), or none.
/// Derived from #Hermiticity, #ConjugationParity and the #Field, see
/// to_braket_symmetry().
/// @note only whole bra<->ket exchanges are described, not the swap of a
///       single particle's bra and ket slot
/// @note The enumerator values enter hashes that decide canonical tie-breaks;
///       new states are appended so the existing values stay fixed.
enum class BraKetSymmetry {
  Symm = 0,
  Conjugate = 1,
  Nonsymm = 2,
  Antisymm = 3,
  AntiConjugate = 4
};

/// relation between the complex-conjugated array and the array as written:
/// `conj(T{p;q}) = +T{p;q}` (`Symm`), `−T{p;q}` (`Antisymm`), or none.
/// Derived from #ConjugationParity and the #Field, see
/// to_conjugation_symmetry().
enum class ConjugationSymmetry { Symm, Antisymm, Nonsymm };

/// describes the abstract symmetry of a tensorial object under (Hermitian)
/// adjoint, i.e. whether the abstract tensor equals (`Hermitian`), equals minus
/// (`AntiHermitian`), or is unrelated to (`NonHermitian`) its own adjoint.
///
/// Unlike #BraKetSymmetry this is a field-agnostic property of the abstract
/// tensor: e.g. a 2-electron integral is `Hermitian` whether the computation is
/// real or complex, a cluster amplitude is `NonHermitian` in either. The
/// observable bra<->ket exchange symmetry (#BraKetSymmetry) is the *derived*
/// composition of this trait with the ambient #Field -- see
/// to_braket_symmetry().
enum class Hermiticity { Hermitian, AntiHermitian, NonHermitian };

/// @return the elementwise conjugation symmetry of an array of @p parity over
///         a basis of the given @p field
constexpr ConjugationSymmetry to_conjugation_symmetry(ConjugationParity parity,
                                                      Field field) noexcept {
  if (field == Field::Complex) return ConjugationSymmetry::Nonsymm;
  switch (parity) {
    case ConjugationParity::Even:
      return ConjugationSymmetry::Symm;
    case ConjugationParity::Odd:
      return ConjugationSymmetry::Antisymm;
    case ConjugationParity::None:
      return ConjugationSymmetry::Nonsymm;
  }
  SEQUANT_UNREACHABLE;
}

/// @return the bra<->ket exchange symmetry implied by the traits @p
///         hermiticity and @p parity over a basis of the given @p field. Over
///         a real basis with known parity the conjugation is a sign, so the
///         exchange relation becomes a plain (anti)symmetry.
constexpr BraKetSymmetry to_braket_symmetry(Hermiticity hermiticity,
                                            ConjugationParity parity,
                                            Field field) noexcept {
  if (hermiticity == Hermiticity::NonHermitian) return BraKetSymmetry::Nonsymm;
  const bool hermitian = hermiticity == Hermiticity::Hermitian;
  if (field == Field::Complex || parity == ConjugationParity::None)
    return hermitian ? BraKetSymmetry::Conjugate
                     : BraKetSymmetry::AntiConjugate;
  // real basis, parity known: T{q;p} = conj(±T{p;q}) = (parity sign)(±T{p;q})
  const bool even = parity == ConjugationParity::Even;
  return (hermitian == even) ? BraKetSymmetry::Symm : BraKetSymmetry::Antisymm;
}

/// @return to_braket_symmetry(@p hermiticity, ConjugationParity::Even, @p
/// field)
constexpr BraKetSymmetry to_braket_symmetry(Hermiticity hermiticity,
                                            Field field) noexcept {
  return to_braket_symmetry(hermiticity, ConjugationParity::Even, field);
}

/// @return the #Hermiticity consistent with an explicitly given @p
///         braket_symmetry and @p parity; back-fills the trait when a tensor is
///         constructed from its #BraKetSymmetry directly. `Symm`/`Antisymm`
///         arise only over a real basis, where the parity decides which
///         hermiticity they came from.
constexpr Hermiticity to_hermiticity(
    BraKetSymmetry braket_symmetry,
    ConjugationParity parity = ConjugationParity::Even) noexcept {
  const bool odd = parity == ConjugationParity::Odd;
  switch (braket_symmetry) {
    case BraKetSymmetry::Symm:
      return odd ? Hermiticity::AntiHermitian : Hermiticity::Hermitian;
    case BraKetSymmetry::Antisymm:
      return odd ? Hermiticity::Hermitian : Hermiticity::AntiHermitian;
    case BraKetSymmetry::Conjugate:
      return Hermiticity::Hermitian;
    case BraKetSymmetry::AntiConjugate:
      return Hermiticity::AntiHermitian;
    case BraKetSymmetry::Nonsymm:
      return Hermiticity::NonHermitian;
  }
  SEQUANT_UNREACHABLE;
}

/// @return the #ConjugationParity implied by an explicitly given @p
///         braket_symmetry over a basis of the given @p field; back-fills the
///         trait when a tensor is constructed from its #BraKetSymmetry
///         directly. A (anti)conjugation pinned over a real basis asserts the
///         adjoint relation and no reality, hence parity `None`; every other
///         case carries no statement about the parity and takes the default,
///         `Even`.
constexpr ConjugationParity to_conjugation_parity(
    BraKetSymmetry braket_symmetry, Field field) noexcept {
  if (field == Field::Real &&
      (braket_symmetry == BraKetSymmetry::Conjugate ||
       braket_symmetry == BraKetSymmetry::AntiConjugate))
    return ConjugationParity::None;
  return ConjugationParity::Even;
}

/// @return the sign `s` in `T{q;p} = s T{p;q}` if the exchange is a plain
///         (anti)symmetry, else nullopt
constexpr std::optional<std::int8_t> braket_swap_sign(
    BraKetSymmetry s) noexcept {
  if (s == BraKetSymmetry::Symm) return std::int8_t{1};
  if (s == BraKetSymmetry::Antisymm) return std::int8_t{-1};
  return std::nullopt;
}
/// @return the sign `s` in `T{q;p} = s conj(T{p;q})` if the exchange is a
///         (anti)conjugation, else nullopt
constexpr std::optional<std::int8_t> braket_conjugate_swap_sign(
    BraKetSymmetry s) noexcept {
  if (s == BraKetSymmetry::Conjugate) return std::int8_t{1};
  if (s == BraKetSymmetry::AntiConjugate) return std::int8_t{-1};
  return std::nullopt;
}
/// @return the sign `s` in `conj(T) = s T` if known, else nullopt
constexpr std::optional<std::int8_t> conjugation_sign(
    ConjugationSymmetry s) noexcept {
  if (s == ConjugationSymmetry::Symm) return std::int8_t{1};
  if (s == ConjugationSymmetry::Antisymm) return std::int8_t{-1};
  return std::nullopt;
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
  Bra = 0b1,
  Ket = 0b10,
  Aux = 0b100,
  Proto = 0b1000,
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

struct SlotTypes {
  std::underlying_type_t<SlotType> active = 0;

  constexpr SlotTypes(SlotType type)
      : active(static_cast<decltype(active)>(type)) {}
  constexpr SlotTypes(decltype(active) val) : active(val) {}

  constexpr bool operator==(const SlotTypes&) const = default;
  constexpr auto operator<=>(const SlotTypes&) const = default;

  constexpr bool operator&(SlotType type) const {
    return active & static_cast<decltype(active)>(type);
  }

  constexpr SlotTypes operator|(SlotType type) const {
    return {active | static_cast<decltype(active)>(type)};
  }
};

constexpr SlotTypes operator|(SlotType lhs, SlotType rhs) {
  using IntType = std::underlying_type_t<SlotType>;
  return {static_cast<IntType>(lhs) | static_cast<IntType>(rhs)};
}

static constexpr const SlotTypes AnySlotType =
    SlotType::Bra | SlotType::Ket | SlotType::Aux | SlotType::Proto;

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
