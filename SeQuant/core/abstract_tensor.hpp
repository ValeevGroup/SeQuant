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

/// This interface class defines a Tensor concept. Object @c t of a type that
/// meets the concept must satisfy the following:
///         - @c bra(t) , @c ket(t) , and @c braket(t) are valid expressions and
///         evaluate to a range of Index objects;
///         - @c bra_rank(t) and @c ket_rank(t) are valid expression and return
///         sizes of the @c bra(t) and @c ket(t) ranges, respectively;
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
/// To adapt an existing class intrusively derive it from AbstractTensor and
/// implement all member functions. This allows to implememnt heterogeneous
/// containers of objects that meet the Tensor concept.
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

  /// view of a contiguous range of Index objects
  virtual const_any_view_randsz _bra() const {
    throw missing_instantiation_for("_bra");
  }
  /// view of a contiguous range of Index objects
  virtual const_any_view_randsz _ket() const {
    throw missing_instantiation_for("_ket");
  }
  /// view of a not necessarily contiguous range of Index objects
  virtual const_any_view_rand _braket() const {
    throw missing_instantiation_for("_braket");
  }
  virtual std::size_t _bra_rank() const {
    throw missing_instantiation_for("_bra_rank");
  }
  virtual std::size_t _ket_rank() const {
    throw missing_instantiation_for("_ket_rank");
  }
  virtual Symmetry _symmetry() const {
    throw missing_instantiation_for("_symmetry");
  }
  virtual BraKetSymmetry _braket_symmetry() const {
    throw missing_instantiation_for("_braket_symmetry");
  }
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

  friend class TensorCanonicalizer;
};

/// @name customization points to support generic algorithms on AbstractTensor
/// objects.
/// @{
inline auto bra(const AbstractTensor& t) { return t._bra(); }
inline auto ket(const AbstractTensor& t) { return t._ket(); }
inline auto braket(const AbstractTensor& t) { return t._braket(); }
inline auto bra_rank(const AbstractTensor& t) { return t._bra_rank(); }
inline auto ket_rank(const AbstractTensor& t) { return t._ket_rank(); }
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
/// @tparam IndexMap a {source Index -> target Index} map type; if it is not @c
/// container::map<Index,Index>
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
