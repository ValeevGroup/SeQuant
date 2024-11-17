//
// Created by Eduard Valeyev on 2019-03-22.
//

#ifndef SEQUANT_ABSTRACT_TENSOR_HPP
#define SEQUANT_ABSTRACT_TENSOR_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <memory>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeinfo>
#include <utility>

#include <range/v3/all.hpp>

#include <boost/core/demangle.hpp>

namespace sequant {

class TensorCanonicalizer;

/// AbstractTensor is a [tensor](https://en.wikipedia.org/wiki/Tensor) over
/// general (i.e., not necessarily commutative)
/// rings. A tensor with \f$ k \geq 0 \f$ contravariant
/// (ket, [Dirac notation](https://en.wikipedia.org/wiki/Bra-ket_notation) ) and
/// \f$ b \geq 0 \f$ covariant (bra) modes
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
///   of corresponding (first, second, etc.) modes in bra/ket mode ranges. This
///   symmetry is used to treat particle as indistinguishable or distinguishable
///   in many-particle quantum mechanics context.
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
/// Lastly, the supporting rings are not assumed to be scalars, hence tensor
/// product supporting the concept of tensor is not necessarily commutative.
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

  /// accessor bra (covariant) indices
  /// @return view of a contiguous range of Index objects
  virtual const_any_view_randsz _bra() const {
    throw missing_instantiation_for("_bra");
  }
  /// accesses ket (contravariant) indices
  /// @return view of a contiguous range of Index objects
  virtual const_any_view_randsz _ket() const {
    throw missing_instantiation_for("_ket");
  }
  /// accesses aux (invariant) indices
  /// @return view of a contiguous range of Index objects
  virtual const_any_view_randsz _aux() const {
    throw missing_instantiation_for("_aux");
  }
  /// accesses bra and ket indices
  /// view of a not necessarily contiguous range of Index objects
  virtual const_any_view_rand _braket() const {
    throw missing_instantiation_for("_braket");
  }
  /// accesses bra, ket, and aux indices
  /// @return view of a not necessarily contiguous range of Index objects
  virtual const_any_view_rand _indices() const {
    throw missing_instantiation_for("_indices");
  }
  /// @return the number of bra indices
  virtual std::size_t _bra_rank() const {
    throw missing_instantiation_for("_bra_rank");
  }
  /// @return the number of ket indices
  virtual std::size_t _ket_rank() const {
    throw missing_instantiation_for("_ket_rank");
  }
  /// @return the number of aux indices
  virtual std::size_t _aux_rank() const {
    throw missing_instantiation_for("_aux_rank");
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
  /// @return the symmetry of tensor under exchange of bra and ket
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

  virtual bool operator<(const AbstractTensor& other) const {
    throw missing_instantiation_for("operator<");
  }

  virtual bool _transform_indices(
      const container::map<Index, Index>& index_map) {
    throw missing_instantiation_for("_transform_indices");
  }
  virtual void _reset_tags() { throw missing_instantiation_for("_reset_tags"); }

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
};

/// @name customization points to support generic algorithms on AbstractTensor
/// objects.
/// @{
inline auto bra(const AbstractTensor& t) { return t._bra(); }
inline auto ket(const AbstractTensor& t) { return t._ket(); }
inline auto aux(const AbstractTensor& t) { return t._aux(); }
inline auto braket(const AbstractTensor& t) { return t._braket(); }
inline auto indices(const AbstractTensor& t) { return t._indices(); }
inline auto bra_rank(const AbstractTensor& t) { return t._bra_rank(); }
inline auto ket_rank(const AbstractTensor& t) { return t._ket_rank(); }
inline auto aux_rank(const AbstractTensor& t) { return t._aux_rank(); }
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

/// Type trait for checking whether a given class fulfills the Tensor interface
/// requirements Object @c t of a type that meets the concept must satisfy the
/// following:
///         - @c bra(t) , @c ket(t), @c aux(t), @c braket(t) and
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
          std::is_invocable_v<decltype(bra), T> &&
          std::is_invocable_v<decltype(ket), T> &&
          std::is_invocable_v<decltype(aux), T> &&
          std::is_invocable_v<decltype(braket), T> &&
          std::is_invocable_v<decltype(indices), T> &&
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

// defined in AbstractTensor
// inline bool operator<(const AbstractTensor& first, const AbstractTensor&
// second) {
//  return first.operator<(second);
//}

///@}

using AbstractTensorPtr = std::shared_ptr<AbstractTensor>;

}  // namespace sequant

#endif  // SEQUANT_ABSTRACT_TENSOR_HPP
