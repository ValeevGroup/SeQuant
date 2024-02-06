//
// Created by Robert Adam on 2023-09-08
//

#ifndef SEQUANT_CORE_TENSOR_CANONICALIZER_HPP
#define SEQUANT_CORE_TENSOR_CANONICALIZER_HPP

#include "abstract_tensor.hpp"
#include "expr.hpp"

#include <memory>
#include <string_view>

namespace sequant {

/// @brief Base class for Tensor canonicalizers
/// To make custom canonicalizer make a derived class and register an instance
/// of that class with TensorCanonicalizer::register_instance
class TensorCanonicalizer {
 public:
	 using index_comparer_t = std::function<bool(const Index&, const Index&)>;
	 using index_pair_t = std::pair<const Index, const Index>;
	 using index_pair_comparer_t = std::function<bool(const index_pair_t &, const index_pair_t)>;

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
  /// @note canonicalization compared indices returned by index_comparer
  // TODO generalize for complex tensors
  virtual ExprPtr apply(AbstractTensor&) const = 0;

  /// @return reference to the object used to compare Index objects
  static const index_comparer_t& index_comparer();

  /// @param comparer the compare object to be used by this
  static void index_comparer(index_comparer_t comparer);

  /// @return reference to the object used to compare Index objects
  static const index_pair_comparer_t& index_pair_comparer();

  /// @param comparer the compare object to be used by this
  static void index_pair_comparer( index_pair_comparer_t comparer);

 protected:
  inline auto bra_range(AbstractTensor& t) const { return t._bra_mutable(); }
  inline auto ket_range(AbstractTensor& t) const { return t._ket_mutable(); }
  inline auto auxiliary_range(AbstractTensor& t) const {
    return t._auxiliary_mutable();
  }

  /// the object used to compare indices
  static index_comparer_t index_comparer_;
  /// the object used to compare pairs of indices
  static index_pair_comparer_t index_pair_comparer_;

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

  ExprPtr apply(AbstractTensor&) const override;
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
  ExprPtr apply(AbstractTensor& t) const override;

  /// Core of DefaultTensorCanonicalizer::apply, only does the canonicalization,
  /// i.e. no tagging/untagging
  template <typename IndexComp, typename IndexPairComp>
  ExprPtr apply(AbstractTensor& t, const IndexComp &idxcmp, const IndexPairComp &paircmp) const {
    // std::wcout << "abstract tensor: " << to_latex(t) << "\n";
    auto s = symmetry(t);
    auto is_antisymm = (s == Symmetry::antisymm);
    const auto _bra_rank = bra_rank(t);
    const auto _ket_rank = ket_rank(t);
    const auto _aux_rank = auxiliary_rank(t);
    const auto _rank = std::min(_bra_rank, _ket_rank);

    // nothing to do for rank-1 tensors
    if (_bra_rank == 1 && _ket_rank == 1 && _aux_rank == 1) return nullptr;

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
        bubble_sort(begin(_bra), end(_bra), idxcmp);
        bubble_sort(begin(_ket), end(_ket), idxcmp);
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
        bubble_sort(begin(_zip_braket), end(_zip_braket), paircmp);
        if (_bra_rank > _rank) {
          auto size_of_rest = _bra_rank - _rank;
          auto rest_of = counted(begin(_bra) + _rank, size_of_rest);
          bubble_sort(begin(rest_of), end(rest_of), idxcmp);
        } else if (_ket_rank > _rank) {
          auto size_of_rest = _ket_rank - _rank;
          auto rest_of = counted(begin(_ket) + _rank, size_of_rest);
          bubble_sort(begin(rest_of), end(rest_of), idxcmp);
        }
      } break;

      default:
        abort();
    }

	// TODO: Handle auxiliary index symmetries once they are introduced
    // auto _aux = auxiliary_range(t);
    // ranges::sort(_aux, comp);

    ExprPtr result =
        is_antisymm ? (even == false ? ex<Constant>(-1) : nullptr) : nullptr;
    return result;
  }

 private:
  container::map<std::wstring, Index> external_indices_;

 protected:
  void tag_indices(AbstractTensor& t) const;
};

class TensorBlockCanonicalizer : public DefaultTensorCanonicalizer {
 public:
  TensorBlockCanonicalizer() = default;
  ~TensorBlockCanonicalizer() = default;

  template <typename IndexContainer>
  TensorBlockCanonicalizer(const IndexContainer& external_indices)
      : DefaultTensorCanonicalizer(external_indices) {}

  ExprPtr apply(AbstractTensor& t) const override;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_TENSOR_CANONICALIZER_HPP
