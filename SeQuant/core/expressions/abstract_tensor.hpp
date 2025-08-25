//
// Created by Eduard Valeyev on 2019-03-22.
//

#ifndef SEQUANT_ABSTRACT_TENSOR_HPP
#define SEQUANT_ABSTRACT_TENSOR_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/index.hpp>

#include <cstdlib>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeinfo>

#include <range/v3/all.hpp>

#include <boost/core/demangle.hpp>

namespace sequant {

/// index slot types
///
/// @note This does not include slot bundles, like braket, etc.
enum class SlotType { Bra, Ket, Aux };

class TensorCanonicalizer;

/// AbstractTensor is a [tensor](https://en.wikipedia.org/wiki/Tensor) over
/// general (i.e., not necessarily commutative)
/// rings. A tensor with \f$ k \geq 0 \f$ ket
/// (see [Dirac notation](https://en.wikipedia.org/wiki/Bra-ket_notation) ) and
/// \f$ b \geq 0 \f$ bra modes
/// describes elements of a tensor product of \f$ b+k \f$ vector spaces.
/// Equivalently it represents a linear map between the tensor product
/// of \f$ k \f$ _primal_ vector spaces to the tensor product of \f$ b \f$
/// _dual_ vector spaces. Tensor modes are 1-to-1 represented by unique
/// [indices](https://en.wikipedia.org/wiki/Abstract_index_notation),
/// represented by Index objects.
///
/// It is also necessary to support modes that are "array-like" in that they
/// do not refer to a vector space or its dual; such modes can represent
/// ordinal indices (e.g. to treat a collection/sequence of tensors as
/// a single tensor) Thus each tensor has zero or more auxiliary (aux) modes.
/// The aux modes are invariant under the transformations of vector spaces and
/// do not contribute to the tensor symmetries.
///
/// Tensors can have the following symmetries:
/// - Tensors can be symmetric or nonsymmetric with respect to the transposition
///   of corresponding {bra,ket} modes in bra/ket mode ranges. This
///   symmetry is used to treat particles as indistinguishable or
///   distinguishable in physical simulation. See more below on how such mode
///   _bundles_ are handled.
/// - Tensors can be symmetric, antisymmetric, and nonsymmetric
///   with respect to the transposition of modes within the bra or ket sets.
///   This symmetry is used to model the
///   distinguishable and indistinguishable (bosonic and fermionic) degrees
///   of freedom in many-body quantum mechanics context. More complicated
///   symmetries are not yet supported.
/// - Tensors can be symmetric, conjugate, or nonsymmetric with respect to
///   swap of bra with ket. This symmetry corresponds to time reversal in
///   physical simulation.
///
/// The supporting rings are not assumed to be scalars, hence tensor
/// product supporting the concept of tensor is not necessarily commutative.
///
/// For some applications the pairing of bra and ket modes is important, hence
/// it is necessary to support not only tensors where first bra mode is paired
/// with the first ket, but it is then necessary to support bra and ket
/// modes that are unpaired. This introduces notion of bra/ket _slots_.
/// A slot is either occupied by a non-null Index (hence corresponds to a bra
/// or ket mode), or empty (i.e., occupied by a null Index).
/// Tensors with non-symmetric bras and kets can have empty slots in arbitrary
/// positions in bra and ket slot bundles. Such tensors symmetric with
/// respect to reordering of braket bundles (`_particle_symmetry()`)
/// assume the canonical order of slots, defined as follows:
/// - paired braket slot bundles (both bra and ket slots are nonempty) appear
///    first
/// - unpaired bra slots appear next (with ket slots empty, or without
///   altogether is there are no unpaired ket slots)
/// - unpaired ket slots appear next (without matching bra slots).
///
/// Empty aux slots are not permitted.
///
/// \note This interface class defines a Tensor _concept_. All Tensor objects
/// must fulfill the is_tensor trait (see below). To adapt an existing class
/// intrusively derive it from AbstractTensor and implement all member
/// functions. This allows to implement heterogeneous containers of objects
/// that meet the Tensor concept.
class AbstractTensor {
  inline auto missing_instantiation_for(const char* fn_name) const {
    std::ostringstream oss;
    oss << "AbstractTensor::" << fn_name << " not implemented in class "
        << boost::core::demangle(typeid(*this).name());
    return std::runtime_error(oss.str());
  }

 public:
  virtual ~AbstractTensor() = default;

  /// clones this object
  //// @throw missing_instantiation_for if not overloaded
  virtual AbstractTensor* _clone() const {
    throw missing_instantiation_for("_clone");
  }

  using const_any_view_rand =
      ranges::any_view<const Index&, ranges::category::random_access>;
  using const_any_view_randsz =
      ranges::any_view<const Index&, ranges::category::random_access |
                                         ranges::category::sized>;
  using any_view_rand =
      ranges::any_view<Index&, ranges::category::random_access>;
  using any_view_randsz =
      ranges::any_view<Index&, ranges::category::random_access |
                                   ranges::category::sized>;

  /// accessor for bra slots
  /// @return view of a contiguous range of bra Index objects; empty slots
  /// may be empty (i.e., occupied by null indices).
  virtual const_any_view_randsz _bra() const {
    throw missing_instantiation_for("_bra");
  }
  /// accesses ket slots
  /// @return view of a contiguous range of ket Index objects; empty slots
  /// may be empty (i.e., occupied by null indices).
  virtual const_any_view_randsz _ket() const {
    throw missing_instantiation_for("_ket");
  }
  /// accesses aux (non-vector-space) indices
  /// @return view of a contiguous range of aux Index objects
  virtual const_any_view_randsz _aux() const {
    throw missing_instantiation_for("_aux");
  }
  /// accesses bra and ket slots
  /// @return concatenated (hence, non contiguous) view of bra and ket slots
  virtual const_any_view_rand _braket() const {
    throw missing_instantiation_for("_braket");
  }
  /// accesses bra, ket, and aux indices
  /// @return view of a not necessarily contiguous range of Index objects
  virtual const_any_view_rand _braketaux() const {
    throw missing_instantiation_for("_braketaux");
  }
  /// accesses all slots
  /// @return view of a not necessarily contiguous range of Index objects
  virtual const_any_view_rand _slots() const {
    throw missing_instantiation_for("_slots");
  }
  /// @return the number of bra slots (some may be empty, hence this is the
  /// gross rank)
  virtual std::size_t _bra_rank() const {
    throw missing_instantiation_for("_bra_rank");
  }
  /// @return the number of nonempty bra slots (i.e., nonnull bra indices)
  virtual std::size_t _bra_net_rank() const {
    throw missing_instantiation_for("_bra_net_rank");
  }
  /// @return the number of bra slots (some may be empty, hence this is the
  /// gross rank)
  virtual std::size_t _ket_rank() const {
    throw missing_instantiation_for("_ket_rank");
  }
  /// @return the number of nonempty ket slots (i.e., nonnull bra indices)
  virtual std::size_t _ket_net_rank() const {
    throw missing_instantiation_for("_ket_net_rank");
  }
  /// @return the number of aux slots
  virtual std::size_t _aux_rank() const {
    throw missing_instantiation_for("_aux_rank");
  }
  /// @return the number of slots
  virtual std::size_t _num_slots() const {
    throw missing_instantiation_for("_num_slots");
  }
  /// @return the number of indices
  virtual std::size_t _num_indices() const {
    throw missing_instantiation_for("_num_indices");
  }
  /// @return the permutational symmetry of the vector space indices of
  /// the tensor
  virtual Symmetry _symmetry() const {
    throw missing_instantiation_for("_symmetry");
  }
  /// @return the symmetry of tensor under exchange of vectors space (bra) and
  /// its dual (ket)
  virtual BraKetSymmetry _braket_symmetry() const {
    throw missing_instantiation_for("_braket_symmetry");
  }
  /// @return the symmetry of tensor under exchange of matching {bra,ket} slot
  /// pairs
  /// @note slots are left-aligned
  virtual ParticleSymmetry _particle_symmetry() const {
    throw missing_instantiation_for("_particle_symmetry");
  }
  virtual std::size_t _color() const {
    throw missing_instantiation_for("_color");
  }
  /// @return true if this is a c-number; if false,
  /// @sa Expr::is_cnumber()
  virtual bool _is_cnumber() const {
    throw missing_instantiation_for("_is_cnumber");
  }
  virtual std::wstring_view _label() const {
    throw missing_instantiation_for("_label");
  }
  virtual std::wstring _to_latex() const {
    throw missing_instantiation_for("_to_latex");
  }

  virtual bool operator<(const AbstractTensor&) const {
    throw missing_instantiation_for("operator<");
  }

  /// hashes the tensor's identity , i.e., this includes hashes of the index
  /// labels (ordinals)
  /// @return hash of the tensor's identity
  /// @warning this hash value is not invariant with respect to tensor
  /// symmetries, i.e. permuting bra indices will not leave the hash invariant
  /// even if the tensor is
  virtual std::size_t _hash_value() const {
    throw missing_instantiation_for("_hash_value");
  }

  virtual bool _transform_indices(const container::map<Index, Index>&) {
    throw missing_instantiation_for("_transform_indices");
  }
  virtual void _reset_tags() { throw missing_instantiation_for("_reset_tags"); }

  /// permutes bra slots according to @p perm
  /// @param perm from-permutation, i.e. Index in slot `permutation[i]` will be
  /// in slot `i`
  virtual void _permute_bra(std::span<const std::size_t> perm) {
    permute_impl(_bra_mutable(), perm);
  }
  /// permutes ket slots according to @p perm
  /// @param perm from-permutation, i.e. Index in slot `permutation[i]` will be
  /// in slot `i`
  virtual void _permute_ket(std::span<const std::size_t> perm) {
    permute_impl(_ket_mutable(), perm);
  }
  /// permutes aux slots according to @p perm
  /// @param perm from-permutation, i.e. Index in slot `permutation[i]` will be
  /// in slot `i`
  virtual void _permute_aux(std::span<const std::size_t> perm) {
    permute_impl(_aux_mutable(), perm);
  }
  /// permutes braket slot groups according to @p perm
  /// @param perm from-permutation, i.e. Index pair in slot `permutation[i]`
  /// will end up in slot `i`
  virtual void _permute_braket(std::span<std::size_t> perm) {
    permute_braket_impl(_bra_mutable(), _ket_mutable(), perm);
  }
  /// swaps bra and ket slots
  virtual void _swap_bra_ket() {
    throw missing_instantiation_for("_swap_bra_ket");
  }

 private:
  /// @return mutable view of bra
  /// @warning this is used for mutable access, flush memoized state before
  /// returning!
  virtual any_view_randsz _bra_mutable() {
    throw missing_instantiation_for("_bra_mutable");
  }
  /// @return mutable view to ket
  /// @warning this is used for mutable access, flush memoized state before
  /// returning!
  virtual any_view_randsz _ket_mutable() {
    throw missing_instantiation_for("_ket_mutable");
  }
  /// @return mutable view of aux indices
  /// @warning this is used for mutable access, flush memoized state before
  /// returning!
  virtual any_view_randsz _aux_mutable() {
    throw missing_instantiation_for("_aux_mutable");
  }

  friend class TensorCanonicalizer;

  static void permute_impl(AbstractTensor::any_view_randsz indices,
                           std::span<const std::size_t> perm) {
    const auto n = indices.size();
    assert(indices.size() == perm.size());
    container::svector<Index> sorted_indices(n);
    for (std::size_t i = 0; i != n; ++i) {
      sorted_indices[i] = std::move(indices[perm[i]]);
    }
    for (std::size_t i = 0; i != n; ++i) {
      indices[i] = std::move(sorted_indices[i]);
    }
  }

  static void permute_braket_impl(AbstractTensor::any_view_randsz bra_indices,
                                  any_view_randsz ket_indices,
                                  std::span<std::size_t> perm_from) {
    const auto n = std::max(bra_indices.size(), ket_indices.size());
    assert(n == perm_from.size());

    // N.B. braket slot bundles are kept in canonical order, see AbstractTensor
    // class dox: paired bundles first, then bra (paired with null ket) and
    // ket (unpaired).

    // 1. count each type of bundles
    std::size_t npaired = 0;
    std::size_t nunpaired_bra = 0;
    for (const auto& [bra_idx, ket_idx] :
         ranges::views::zip(bra_indices, ket_indices)) {
      if (bra_idx.nonnull()) {
        if (ket_idx.nonnull())
          ++npaired;
        else
          ++nunpaired_bra;
      }
    }
    // corner case: there are only unpaired bra slots, no unpaired ket slots
    if (bra_indices.size() > ket_indices.size()) {
      assert(ket_indices.size() == npaired);
      nunpaired_bra += bra_indices.size() - ket_indices.size();
    }

    // adjust perm to ensure that different types of bundles do not mix
    // this ends up just a simple stable sort according to the bundle type
    // i.e. bring elements of perm_from with values [0,npaired) to the front,
    // etc.
    ranges::stable_sort(perm_from, [&](const auto& i, const auto j) {
      enum { paired = 0, bra_unpaired = 1, ket_unpaired = 2 };
      auto to_type = [&](const auto& i) {
        if (i < npaired)
          return paired;
        else if (i < npaired + nunpaired_bra)
          return bra_unpaired;
        else
          return ket_unpaired;
      };
      return to_type(i) < to_type(j);
    });

    container::svector<Index> sorted_indices(n);
    for (std::size_t i = 0; i != bra_indices.size(); ++i) {
      assert(perm_from[i] < bra_indices.size());
      sorted_indices[i] = std::move(bra_indices[perm_from[i]]);
    }
    for (std::size_t i = 0; i != bra_indices.size(); ++i) {
      bra_indices[i] = std::move(sorted_indices[i]);
    }
    for (std::size_t i = 0; i != ket_indices.size(); ++i) {
      assert(perm_from[i] < ket_indices.size());
      sorted_indices[i] = std::move(ket_indices[perm_from[i]]);
    }
    for (std::size_t i = 0; i != ket_indices.size(); ++i) {
      ket_indices[i] = std::move(sorted_indices[i]);
    }
  }
};

/// @name customization points to support generic algorithms on AbstractTensor
/// objects.
/// @{
inline auto braket(const AbstractTensor& t) { return t._braket(); }
inline auto braketaux(const AbstractTensor& t) { return t._slots(); }
/// @return a view of all slots
inline auto slots(const AbstractTensor& t) { return t._slots(); }
inline auto bra_rank(const AbstractTensor& t) { return t._bra_rank(); }
inline auto bra_net_rank(const AbstractTensor& t) { return t._bra_net_rank(); }
inline auto ket_rank(const AbstractTensor& t) { return t._ket_rank(); }
inline auto ket_net_rank(const AbstractTensor& t) { return t._ket_net_rank(); }
inline auto aux_rank(const AbstractTensor& t) { return t._aux_rank(); }
inline auto num_indices(const AbstractTensor& t) { return t._num_indices(); }
inline auto num_slots(const AbstractTensor& t) { return t._num_slots(); }
inline auto symmetry(const AbstractTensor& t) { return t._symmetry(); }
inline auto braket_symmetry(const AbstractTensor& t) {
  return t._braket_symmetry();
}
inline auto particle_symmetry(const AbstractTensor& t) {
  return t._particle_symmetry();
}
inline auto color(const AbstractTensor& t) { return t._color(); }
inline auto is_cnumber(const AbstractTensor& t) { return t._is_cnumber(); }
inline auto label(const AbstractTensor& t) { return t._label(); }
inline auto to_latex(const AbstractTensor& t) { return t._to_latex(); }

/// produces LaTeX representation typeset using <a
/// href="https://ctan.org/pkg/tensor?lang=en">tensor package</a>

/// @param[in] core_label the tensor core label
/// @param[in] bra the tensor bra indices
/// @param[in] ket the tensor ket indices
/// @param[in] aux the tensor aux indices
/// @param[in] bkt the typesetting convention for bra/ket (super or subscript)
/// @param[in] left_align if true, will typeset `bra[0]` above `ket[0]`,
/// `bra[1]` above `ket[1]`, etc.; else will typeset `bra[-1]` above `ket[-1]`,
/// etc.
inline std::wstring to_latex_tensor(
    const std::wstring& core_label, AbstractTensor::const_any_view_randsz bra,
    AbstractTensor::const_any_view_randsz ket,
    AbstractTensor::const_any_view_randsz aux,
    const BraKetTypesetting bkt = BraKetTypesetting::ContraSub,
    bool left_align = true) {
  std::wstring result = L"{\\tensor*{" + core_label + L"}{";

  const auto braket_rank_max = std::max(bra.size(), ket.size());
  const auto num_paired = std::min(bra.size(), ket.size());
  const auto num_unpaired = braket_rank_max - num_paired;
  const auto unpaired_type =
      bra.size() > ket.size() ? SlotType::Bra : SlotType::Ket;

  std::size_t col = 0;

  // loop over left-aligned unpaired slots, if left_align==false
  if (left_align == false) {
    auto* unpaired_indices = unpaired_type == SlotType::Bra
                                 ? &bra[0]
                                 : &ket[0];  // bra/ket are contiguous
    for (; col != num_unpaired; ++col) {
      result += L"*^";
      if ((bkt == BraKetTypesetting::BraSuper &&
           unpaired_type == SlotType::Bra) ||
          (bkt == BraKetTypesetting::KetSuper &&
           unpaired_type == SlotType::Ket)) {
        result += to_latex(unpaired_indices[col]);
      } else
        result += L"{}";
      result += L"_";
      if ((bkt == BraKetTypesetting::BraSub &&
           unpaired_type == SlotType::Bra) ||
          (bkt == BraKetTypesetting::KetSub &&
           unpaired_type == SlotType::Ket)) {
        result += to_latex(unpaired_indices[col]);
      } else
        result += L"{}";
    }
  }

  // loop over paired indices
  auto paired_bra =
      unpaired_type == SlotType::Bra && left_align == false ? col : 0;
  auto paired_ket =
      unpaired_type == SlotType::Ket && left_align == false ? col : 0;
  const auto paired_fence = col + num_paired;
  for (; col != paired_fence; ++col, ++paired_bra, ++paired_ket) {
    result += L"*^";
    result += (bkt == BraKetTypesetting::BraSuper) ? to_latex(bra[paired_bra])
                                                   : to_latex(ket[paired_ket]);
    result += L"_";
    result += (bkt == BraKetTypesetting::BraSub) ? to_latex(bra[paired_bra])
                                                 : to_latex(ket[paired_ket]);
  }

  // loop over right-aligned unpaired slots, if left_align==true
  if (left_align == true) {
    auto* unpaired_indices = unpaired_type == SlotType::Bra
                                 ? &bra[0]
                                 : &ket[0];  // bra/ket are contiguous
    for (; col != braket_rank_max; ++col) {
      result += L"*^";
      if ((bkt == BraKetTypesetting::BraSuper &&
           unpaired_type == SlotType::Bra) ||
          (bkt == BraKetTypesetting::KetSuper &&
           unpaired_type == SlotType::Ket)) {
        result += to_latex(unpaired_indices[col]);
      } else
        result += L"{}";
      result += L"_";
      if ((bkt == BraKetTypesetting::BraSub &&
           unpaired_type == SlotType::Bra) ||
          (bkt == BraKetTypesetting::KetSub &&
           unpaired_type == SlotType::Ket)) {
        result += to_latex(unpaired_indices[col]);
      } else
        result += L"{}";
    }
  }
  result += L"}";

  // aux
  if (aux.size() != 0) {
    result += L"[";
    const auto aux_rank = aux.size();
    for (std::size_t i = 0; i < aux_rank; ++i) {
      result += sequant::to_latex(aux[i]);
      if (i + 1 < aux_rank) {
        result += L",";
      }
    }
    result += L"]";
  }

  result += L"}";
  return result;
}

/// Type trait for checking whether a given class fulfills the Tensor interface
/// requirements Object @c t of a type that meets the concept must satisfy the
/// following:
///         - @c braket(t) and
///         @c indices(t) are valid expressions and evaluate to a range of Index
///         objects;
///         - @c bra_rank(t), @c ket_rank(t) and @c aux_rank(t) are valid
///         expression and return sizes of the @c bra(t), @c ket(t) and
///         @c aux(t) ranges, respectively;
///         - @c symmetry(t) is a valid expression and evaluates to a Symmetry
///         object that describes the symmetry of bra/ket of a
///         _particle-symmetric_ @c t ;
///         - @c braket_symmetry(t) is a valid expression and evaluates to a
///         BraKetSymmetry object that describes the bra-ket symmetry of @c t ;
///         - @c particle_symmetry(t) is a valid expression and evaluates to a
///         ParticleSymmetry object that describes the symmetry of @c t with
///         respect to permutations of particles;
///         - @c color(t) is a valid expression and returns whether a
///         nonnegative integer that identifies the type of a tensor; tensors
///         with different colors can be reordered in a Product at will
///         - @c is_cnumber(t) is a valid expression and returns whether t
///         commutes with other tensor of same color (tensors of different
///         colors are, for now, always assumed to commute)
///         - @c label(t) is a valid expression and its return is convertible to
///         a std::wstring;
///         - @c to_latex(t) is a valid expression and its return is convertible
///         to a std::wstring.
template <typename T>
struct is_tensor
    : std::bool_constant<
          std::is_invocable_v<decltype(braket), T> &&
          std::is_invocable_v<decltype(braketaux), T> &&
          std::is_invocable_v<decltype(bra_rank), T> &&
          std::is_invocable_v<decltype(ket_rank), T> &&
          std::is_invocable_v<decltype(aux_rank), T> &&
          std::is_invocable_v<decltype(symmetry), T> &&
          std::is_invocable_v<decltype(braket_symmetry), T> &&
          std::is_invocable_v<decltype(particle_symmetry), T> &&
          std::is_invocable_v<decltype(color), T> &&
          std::is_invocable_v<decltype(is_cnumber), T> &&
          std::is_invocable_v<decltype(label), T> &&
          std::is_invocable_v<
              decltype(static_cast<std::wstring (*)(const T&)>(to_latex)), T>> {
};
template <typename T>
constexpr bool is_tensor_v = is_tensor<T>::value;
static_assert(is_tensor_v<AbstractTensor>,
              "The AbstractTensor class does not fulfill the requirements of "
              "the Tensor interface");

/// @tparam IndexMap a {source Index -> target Index} map type; if it is not
///         @c container::map<Index,Index>
///         will need to make a copy.
/// @param[in,out] t an AbstractTensor object whose indices will be transformed
/// @param[in] index_map a const reference to an IndexMap object that specifies
/// the transformation
/// @return false if no indices were transformed, true otherwise
/// @pre indices are not tagged, or (if want to protect them from replacement)
/// tagged with (int)0
/// @post transformed indices are tagged with (int)0
template <typename IndexMap = container::map<Index, Index>>
inline bool transform_indices(AbstractTensor& t, const IndexMap& index_map) {
  if constexpr (std::is_same_v<IndexMap, container::map<Index, Index>>) {
    return t._transform_indices(index_map);
  } else {
    container::map<Index, Index> index_map_copy;
    ranges::copy(index_map, index_map_copy);
    return t._transform_indices(index_map_copy);
  }
}

/// Removes tags from tensor indices
/// @param[in,out] t an AbstractTensor object whose indices will be untagged
inline void reset_tags(AbstractTensor& t) { t._reset_tags(); }

/// permutes bra slots of @p t according to @p perm
/// @param t reference to an AbstractTensor object
/// @param perm from-permutation, i.e. Index in input slot `permutation[i]` will
/// be in slot `i`
inline void permute_bra(AbstractTensor& t, std::span<const std::size_t> perm) {
  return t._permute_bra(perm);
}

/// permutes ket slots of @p t according to @p perm
/// @param t reference to an AbstractTensor object
/// @param perm from-permutation, i.e. Index in input slot `permutation[i]` will
/// be in slot `i`
inline void permute_ket(AbstractTensor& t, std::span<const std::size_t> perm) {
  return t._permute_ket(perm);
}

/// permutes braket slot groups of @p t according to @p perm
/// @param t reference to an AbstractTensor object
/// @param perm from-permutation, i.e. Index pair in input slot `permutation[i]`
/// will be in slot `i`
inline void permute_braket(AbstractTensor& t, std::span<std::size_t> perm) {
  return t._permute_braket(perm);
}

// defined in AbstractTensor
// inline bool operator<(const AbstractTensor& first, const AbstractTensor&
// second) {
//  return first.operator<(second);
//}

///@}

using AbstractTensorPtr = std::shared_ptr<AbstractTensor>;

}  // namespace sequant

#endif  // SEQUANT_ABSTRACT_TENSOR_HPP
