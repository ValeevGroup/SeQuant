//
// Created by Eduard Valeyev on 2019-03-22.
//

#ifndef SEQUANT_ABSTRACT_TENSOR_HPP
#define SEQUANT_ABSTRACT_TENSOR_HPP

#include <range/v3/all.hpp>

#include <boost/core/demangle.hpp>

#include "algorithm.hpp"
#include "expr.hpp"
#include "index.hpp"

#include <thread>

namespace sequant {

class TensorCanonicalizer;

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
  virtual std::wstring _label() const {
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

/// @brief Base class for Tensor canonicalizers
/// To make custom canonicalizer make a derived class and register an instance
/// of that class with TensorCanonicalizer::register_instance
class TensorCanonicalizer {
 public:
  virtual ~TensorCanonicalizer();

  /// @return ptr to the TensorCanonicalizer object, if any, that had been
  /// previously registered via TensorCanonicalizer::register_instance()
  /// with @c label , or to the default canonicalizer, if any
  static std::shared_ptr<TensorCanonicalizer> instance_ptr(
      std::wstring_view label = L"");

  /// @return ptr to the TensorCanonicalizer object, if any, that had been
  /// previously registered via TensorCanonicalizer::register_instance()
  /// with @c label
  /// @sa instance_ptr
  static std::shared_ptr<TensorCanonicalizer> nondefault_instance_ptr(
      std::wstring_view label);

  /// @return a TensorCanonicalizer previously registered via
  /// TensorCanonicalizer::register_instance() with @c label or to the default
  /// canonicalizer
  /// @throw std::runtime_error if no canonicalizer has been registered
  static std::shared_ptr<TensorCanonicalizer> instance(
      std::wstring_view label = L"");

  /// registers @c canonicalizer to be applied to Tensor objects with label
  /// @c label ; leave the label empty if @c canonicalizer is to apply to Tensor
  /// objects with any label
  /// @note if a canonicalizer registered with label @c label exists, it is
  /// replaced
  static void register_instance(
      std::shared_ptr<TensorCanonicalizer> canonicalizer,
      std::wstring_view label = L"");

  /// tries to register @c canonicalizer to be applied to Tensor objects
  /// with label @c label ; leave the label empty if @c canonicalizer is to
  /// apply to Tensor objects with any label
  /// @return false if there is already a canonicalizer registered with @c label
  /// @sa regiter_instance
  static bool try_register_instance(
      std::shared_ptr<TensorCanonicalizer> canonicalizer,
      std::wstring_view label = L"");

  /// deregisters canonicalizer (if any) registered previously
  /// to be applied to tensors with label @c label
  static void deregister_instance(std::wstring_view label = L"");

  /// @return a list of Tensor labels with lexicographic preference (in order)
  static const auto& cardinal_tensor_labels() {
    return cardinal_tensor_labels_accessor();
  }

  /// @param cardinal_tensor_labels a list of Tensor labels with lexicographic
  /// preference (in order)
  static void set_cardinal_tensor_labels(
      const container::vector<std::wstring>& labels) {
    cardinal_tensor_labels_accessor() = labels;
  }

  /// @return a side effect of canonicalization (e.g. phase), or nullptr if none
  /// @internal what should be returned if canonicalization requires
  /// complex conjugation? Special ExprPtr type (e.g. ConjOp)? Or the actual
  /// return of the canonicalization?
  // TODO generalize for complex tensors
  virtual ExprPtr apply(AbstractTensor&) = 0;

 protected:
  inline auto bra_range(AbstractTensor& t) { return t._bra_mutable(); }
  inline auto ket_range(AbstractTensor& t) { return t._ket_mutable(); }

 private:
  static std::pair<
      container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>*,
      std::unique_lock<std::recursive_mutex>>
  instance_map_accessor();  // map* + locked recursive mutex
  static container::vector<std::wstring>& cardinal_tensor_labels_accessor();
};

/// @brief null Tensor canonicalizer does nothing
class NullTensorCanonicalizer : public TensorCanonicalizer {
 public:
  virtual ~NullTensorCanonicalizer() = default;

  ExprPtr apply(AbstractTensor&) override;
};

class DefaultTensorCanonicalizer : public TensorCanonicalizer {
 public:
  DefaultTensorCanonicalizer() = default;

  /// @tparam IndexContainer a Container of Index objects such that @c
  /// IndexContainer::value_type is convertible to Index (e.g. this can be
  /// std::vector or std::set , but not std::map)
  /// @param external_indices container of external Index objects
  /// @warning @c external_indices is assumed to be immutable during the
  /// lifetime of this object
  template <typename IndexContainer>
  DefaultTensorCanonicalizer(IndexContainer&& external_indices) {
    ranges::for_each(external_indices, [this](const Index& idx) {
      this->external_indices_.emplace(idx.label(), idx);
    });
  }
  virtual ~DefaultTensorCanonicalizer() = default;

  /// Implements TensorCanonicalizer::apply
  /// @note Canonicalizes @c t by sorting its bra (if @c
  /// t.symmetry()==Symmetry::nonsymm ) or its bra and ket (if @c
  /// t.symmetry()!=Symmetry::nonsymm ),
  ///       with the external indices appearing "before" (smaller particle
  ///       indices) than the internal indices
  ExprPtr apply(AbstractTensor& t) override;

  /// Core of DefaultTensorCanonicalizer::apply, only does the canonicalization,
  /// i.e. no tagging/untagging
  template <typename Compare>
  ExprPtr apply(AbstractTensor& t, const Compare& comp) {
    // std::wcout << "abstract tensor: " << to_latex(t) << "\n";
    auto s = symmetry(t);
    auto is_antisymm = (s == Symmetry::antisymm);
    const auto _bra_rank = bra_rank(t);
    const auto _ket_rank = ket_rank(t);
    const auto _rank = std::min(_bra_rank, _ket_rank);

    // nothing to do for rank-1 tensors
    if (_bra_rank == 1 && _ket_rank == 1) return nullptr;

    using ranges::begin;
    using ranges::end;
    using ranges::views::counted;
    using ranges::views::take;
    using ranges::views::zip;

    bool even = true;
    switch (s) {
      case Symmetry::antisymm:
      case Symmetry::symm: {
        auto _bra = bra_range(t);
        auto _ket = ket_range(t);
        //      std::wcout << "canonicalizing " << to_latex(t);
        IndexSwapper::thread_instance().reset();
        // std::{stable_}sort does not necessarily use swap! so must implement
        // sort outselves .. thankfully ranks will be low so can stick with
        // bubble
        bubble_sort(begin(_bra), end(_bra), comp);
        bubble_sort(begin(_ket), end(_ket), comp);
        if (is_antisymm)
          even = IndexSwapper::thread_instance().even_num_of_swaps();
        //      std::wcout << " is " << (even ? "even" : "odd") << " and
        //      produces " << to_latex(t) << std::endl;
      } break;

      case Symmetry::nonsymm: {
        // sort particles with bra and ket functions first,
        // then the particles with either bra or ket index
        auto _bra = bra_range(t);
        auto _ket = ket_range(t);
        auto _zip_braket = zip(take(_bra, _rank), take(_ket, _rank));
        bubble_sort(begin(_zip_braket), end(_zip_braket), comp);
        if (_bra_rank > _rank) {
          auto size_of_rest = _bra_rank - _rank;
          auto rest_of = counted(begin(_bra) + _rank, size_of_rest);
          bubble_sort(begin(rest_of), end(rest_of), comp);
        } else if (_ket_rank > _rank) {
          auto size_of_rest = _ket_rank - _rank;
          auto rest_of = counted(begin(_ket) + _rank, size_of_rest);
          bubble_sort(begin(rest_of), end(rest_of), comp);
        }
      } break;

      default:
        abort();
    }

    ExprPtr result =
        is_antisymm ? (even == false ? ex<Constant>(-1) : nullptr) : nullptr;
    return result;
  }

 private:
  container::map<std::wstring, Index> external_indices_;
};

}  // namespace sequant

#endif  // SEQUANT_ABSTRACT_TENSOR_HPP
