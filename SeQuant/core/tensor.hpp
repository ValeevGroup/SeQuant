//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT_TENSOR_HPP
#define SEQUANT_TENSOR_HPP

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/utility/strong.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

namespace sequant {

// strong type wrapper for objects associated with bra
DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE(bra);
// strong type wrapper for objects associated with ket
DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE(ket);
// strong type wrapper for objects associated with aux
DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE(aux);

/// @brief a Tensor is an instance of AbstractTensor over a scalar field, i.e.
/// Tensors have commutative addition and product operations
class Tensor : public Expr, public AbstractTensor, public Labeled {
 private:
  using index_container_type = container::svector<Index>;
  static auto make_indices(IndexList indices) { return indices; }
  static auto make_indices(WstrList index_labels) {
    index_container_type result;
    result.reserve(index_labels.size());
    for (const auto &label : index_labels) {
      result.push_back(Index{label});
    }
    return result;
  }
  static auto make_indices(
      std::initializer_list<const wchar_t *> index_labels) {
    index_container_type result;
    result.reserve(index_labels.size());
    for (const auto &label : index_labels) {
      result.push_back(Index{label});
    }
    return result;
  }
  template <typename IndexRange>
  static index_container_type make_indices(IndexRange &&indices) {
    if constexpr (std::is_same_v<index_container_type,
                                 std::decay_t<IndexRange>>) {
      return std::forward<IndexRange>(indices);
    } else {
      using ranges::begin;
      using ranges::end;
      return index_container_type(begin(indices), end(indices));
    }
  }

  /// @return concatenated view of the bra and ket index ranges
  auto braket() { return ranges::views::concat(bra_, ket_); }

  /// @return concatenated view of bra, ket, and aux index ranges
  auto indices() { return ranges::views::concat(bra_, ket_, aux_); }

  /// asserts that @p label is not reserved
  /// @note Tensor with reserved labels are constructed using friends of Tensor
  /// @param label a Tensor label candidate
  void assert_nonreserved_label(std::wstring_view label) const;

  // utility for dispatching to private ctor
  struct reserved_tag {};

  // list of friends who can make Tensor objects with reserved labels
  friend ExprPtr make_overlap(const Index &bra_index, const Index &ket_index);

  template <typename IndexRange1, typename IndexRange2, typename IndexRange3,
            typename = std::enable_if_t<
                (meta::is_statically_castable_v<
                    meta::range_value_t<IndexRange1>,
                    Index>)&&(meta::
                                  is_statically_castable_v<
                                      meta::range_value_t<IndexRange2>,
                                      Index>)&&(meta::
                                                    is_statically_castable_v<
                                                        meta::range_value_t<
                                                            IndexRange3>,
                                                        Index>)>>
  Tensor(std::wstring_view label, const bra<IndexRange1> &bra_indices,
         const ket<IndexRange2> &ket_indices,
         const aux<IndexRange3> &aux_indices, reserved_tag,
         Symmetry s = Symmetry::nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ParticleSymmetry ps = ParticleSymmetry::symm)
      : label_(label),
        bra_(make_indices(bra_indices)),
        ket_(make_indices(ket_indices)),
        aux_(make_indices(aux_indices)),
        symmetry_(s),
        braket_symmetry_(bks),
        particle_symmetry_(ps) {
    validate_symmetries();
  }

  Tensor(std::wstring_view label, bra<index_container_type> &&bra_indices,
         ket<index_container_type> &&ket_indices,
         aux<index_container_type> &&aux_indices, reserved_tag,
         Symmetry s = Symmetry::nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ParticleSymmetry ps = ParticleSymmetry::symm)
      : label_(label),
        bra_(std::move(bra_indices)),
        ket_(std::move(ket_indices)),
        aux_(std::move(aux_indices)),
        symmetry_(s),
        braket_symmetry_(bks),
        particle_symmetry_(ps) {
    validate_symmetries();
  }

 public:
  /// constructs an uninitialized Tensor
  /// @sa Tensor::operator bool()
  Tensor() = default;
  virtual ~Tensor();

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param s the symmetry of bra or ket
  /// @param bks the symmetry with respect to bra-ket exchange
  /// @param ps the symmetry under exchange of particles
  template <
      typename IndexRange1, typename IndexRange2,
      typename = std::enable_if_t<
          (meta::is_statically_castable_v<
              meta::range_value_t<IndexRange1>,
              Index>)&&(meta::
                            is_statically_castable_v<
                                meta::range_value_t<IndexRange2>, Index>)>>
  Tensor(std::wstring_view label, const bra<IndexRange1> &bra_indices,
         const ket<IndexRange2> &ket_indices, Symmetry s = Symmetry::nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ParticleSymmetry ps = ParticleSymmetry::symm)
      : Tensor(label, bra_indices, ket_indices, sequant::aux{}, reserved_tag{},
               s, bks, ps) {
    assert_nonreserved_label(label_);
  }

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param aux_indices list of aux indices (or objects that can be
  /// converted to indices)
  /// @param s the symmetry of bra or ket
  /// @param bks the symmetry with respect to bra-ket exchange
  /// @param ps the symmetry under exchange of particles
  template <typename IndexRange1, typename IndexRange2, typename IndexRange3,
            typename = std::enable_if_t<
                (meta::is_statically_castable_v<
                    meta::range_value_t<IndexRange1>,
                    Index>)&&(meta::
                                  is_statically_castable_v<
                                      meta::range_value_t<IndexRange2>,
                                      Index>)&&(meta::
                                                    is_statically_castable_v<
                                                        meta::range_value_t<
                                                            IndexRange3>,
                                                        Index>)>>
  Tensor(std::wstring_view label, const bra<IndexRange1> &bra_indices,
         const ket<IndexRange2> &ket_indices,
         const aux<IndexRange3> &aux_indices, Symmetry s = Symmetry::nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ParticleSymmetry ps = ParticleSymmetry::symm)
      : Tensor(label, bra_indices, ket_indices, aux_indices, reserved_tag{}, s,
               bks, ps) {
    assert_nonreserved_label(label_);
  }

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param s the symmetry of bra or ket
  /// @param bks the symmetry with respect to bra-ket exchange
  /// @param ps the symmetry under exchange of particles
  Tensor(std::wstring_view label, bra<index_container_type> &&bra_indices,
         ket<index_container_type> &&ket_indices,
         Symmetry s = Symmetry::nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ParticleSymmetry ps = ParticleSymmetry::symm)
      : Tensor(label, std::move(bra_indices), std::move(ket_indices),
               sequant::aux{}, reserved_tag{}, s, bks, ps) {
    assert_nonreserved_label(label_);
  }

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param aux_indices list of aux indices (or objects that can be
  /// converted to indices)
  /// @param s the symmetry of bra or ket
  /// @param bks the symmetry with respect to bra-ket exchange
  /// @param ps the symmetry under exchange of particles
  Tensor(std::wstring_view label, bra<index_container_type> &&bra_indices,
         ket<index_container_type> &&ket_indices,
         aux<index_container_type> &&aux_indices,
         Symmetry s = Symmetry::nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ParticleSymmetry ps = ParticleSymmetry::symm)
      : Tensor(label, std::move(bra_indices), std::move(ket_indices),
               std::move(aux_indices), reserved_tag{}, s, bks, ps) {
    assert_nonreserved_label(label_);
  }

  /// @return true if the Tensor is initialized
  explicit operator bool() const {
    return !label_.empty() && symmetry_ != Symmetry::invalid &&
           braket_symmetry_ != BraKetSymmetry::invalid &&
           particle_symmetry_ != ParticleSymmetry::invalid;
  }

  /// @return "core" label of the tensor
  std::wstring_view label() const override { return label_; }
  /// @return the bra index range
  const auto &bra() const { return bra_; }
  /// @return the ket index range
  const auto &ket() const { return ket_; }
  /// @return the aux index range
  const auto &aux() const { return aux_; }
  /// @return concatenated view of the bra and ket index ranges
  auto braket() const { return ranges::views::concat(bra_, ket_); }
  /// @return concatenated view of all indices of this tensor (bra, ket and
  /// aux)
  auto indices() const { return ranges::views::concat(bra_, ket_, aux_); }
  /// @return view of the bra+ket index ranges
  /// @note this is to work around broken lookup rules
  auto const_braket() const { return this->braket(); }
  /// @return view of all indices
  /// @note this is to work around broken lookup rules
  auto const_indices() const { return this->indices(); }
  /// Returns the Symmetry object describing the symmetry of the bra and ket of
  /// the Tensor, i.e. what effect swapping indices in positions @c i  and @c j
  /// in <em>either bra or ket</em> has on the elements of the Tensor;
  /// Tensor's are <em>always assumed</em> to be particle-symmetric, i.e.
  /// swapping indices in positions @c i and @c j in <b>both bra and ket</b>;
  /// The allowed values are Symmetry::symm, Symmetry::antisymm, and
  /// Symmetry::nonsymm
  /// @return the Symmetry object describing the symmetry of the bra and ket of
  /// the Tensor.
  Symmetry symmetry() const { return symmetry_; }
  /// @return the BraKetSymmetry object describing the symmetry of the Tensor
  /// under exchange of bra and ket.
  BraKetSymmetry braket_symmetry() const { return braket_symmetry_; }
  /// @return the ParticleSymmetry object describing the symmetry of the Tensor
  /// under exchange of particles (columns).
  ParticleSymmetry particle_symmetry() const { return particle_symmetry_; }

  /// @return number of bra indices
  std::size_t bra_rank() const { return bra_.size(); }
  /// @return number of ket indices
  std::size_t ket_rank() const { return ket_.size(); }
  /// @return number of aux indices
  std::size_t aux_rank() const { return aux_.size(); }
  /// @return number of indices in bra/ket
  /// @throw std::logic_error if bra and ket ranks do not match
  std::size_t rank() const {
    if (bra_rank() != ket_rank()) {
      throw std::logic_error("Tensor::rank(): bra rank != ket rank");
    }
    return bra_rank();
  }

  std::wstring to_latex() const override {
    std::wstring result;
    std::vector<std::wstring> labels = {L"g", L"t", L"λ", L"t¹", L"λ¹"};
    bool add_bar =
        ranges::find(labels, this->label()) != labels.end() && this->rank() > 1;

    result = L"{";
    if ((this->symmetry() == Symmetry::antisymm) && add_bar)
      result += L"\\bar{";
    result += utf_to_latex(this->label());
    if ((this->symmetry() == Symmetry::antisymm) && add_bar) result += L"}";
    result += L"^{";
    for (const auto &i : this->ket()) result += sequant::to_latex(i);
    result += L"}_{";
    for (const auto &i : this->bra()) result += sequant::to_latex(i);
    result += L"}";
    if (!this->aux_.empty()) {
      result += L"(";
      const index_container_type &__aux = this->aux();
      for (std::size_t i = 0; i < aux_rank(); ++i) {
        result += sequant::to_latex(__aux[i]);

        if (i + 1 < aux_rank()) {
          result += L",";
        }
      }
      result += L")";
    }
    result += L"}";
    return result;
  }

  ExprPtr canonicalize() override;

  /// @brief adjoint of a Tensor swaps its bra and ket
  virtual void adjoint() override;

  /// Replaces indices using the index map
  /// @param index_map maps Index to Index
  /// @return true if one or more indices changed
  template <template <typename, typename, typename... Args> class Map,
            typename... Args>
  bool transform_indices(const Map<Index, Index, Args...> &index_map) {
    bool mutated = false;
    ranges::for_each(indices(), [&](auto &idx) {
      if (idx.transform(index_map)) mutated = true;
    });
    if (mutated) this->reset_hash_value();
    return mutated;
  }

  type_id_type type_id() const override { return get_type_id<Tensor>(); };

  ExprPtr clone() const override { return ex<Tensor>(*this); }

  void reset_tags() const {
    ranges::for_each(indices(), [](const auto &idx) { idx.reset_tag(); });
  }

  hash_type bra_hash_value() const {
    if (!hash_value_)  // if hash not computed, or reset, recompute
      memoizing_hash();
    return *bra_hash_value_;
  }

 private:
  std::wstring label_{};
  sequant::bra<index_container_type> bra_{};
  sequant::ket<index_container_type> ket_{};
  sequant::aux<index_container_type> aux_{};
  Symmetry symmetry_ = Symmetry::invalid;
  BraKetSymmetry braket_symmetry_ = BraKetSymmetry::invalid;
  ParticleSymmetry particle_symmetry_ = ParticleSymmetry::invalid;
  mutable std::optional<hash_type>
      bra_hash_value_;  // memoized byproduct of memoizing_hash()
  bool is_adjoint_ = false;

  void validate_symmetries() {
    // (anti)symmetric bra or ket makes sense only for particle-symmetric
    // tensors
    if (symmetry_ == Symmetry::symm || symmetry_ == Symmetry::antisymm)
      assert(particle_symmetry_ == ParticleSymmetry::symm);
  }

  hash_type memoizing_hash() const override {
    using std::begin;
    using std::end;
    auto val = hash::range(begin(bra()), end(bra()));
    bra_hash_value_ = val;
    hash::range(val, begin(ket()), end(ket()));
    hash::range(val, begin(aux()), end(aux()));
    hash::combine(val, label_);
    hash::combine(val, symmetry_);
    hash_value_ = val;
    return *hash_value_;
  }
  void reset_hash_value() const override {
    Expr::reset_hash_value();
    bra_hash_value_.reset();
  }

  bool static_equal(const Expr &that) const override {
    const auto &that_cast = static_cast<const Tensor &>(that);
    if (this->label() == that_cast.label() &&
        this->symmetry() == that_cast.symmetry() &&
        this->bra_rank() == that_cast.bra_rank() &&
        this->ket_rank() == that_cast.ket_rank() &&
        this->aux_rank() == that_cast.aux_rank()) {
      // compare hash values first
      if (this->hash_value() ==
          that.hash_value())  // hash values agree -> do full comparison
        return this->bra() == that_cast.bra() && this->ket() == that_cast.ket();
      else
        return false;
    } else
      return false;
  }

  bool static_less_than(const Expr &that) const override {
    if (this == &that) return false;

    const auto &that_cast = static_cast<const Tensor &>(that);
    if (this->label() != that_cast.label()) {
      return this->label() < that_cast.label();
    }

    if (this->bra_rank() != that_cast.bra_rank()) {
      return this->bra_rank() < that_cast.bra_rank();
    }

    if (this->ket_rank() != that_cast.ket_rank()) {
      return this->ket_rank() < that_cast.ket_rank();
    }

    if (this->aux_rank() != that_cast.aux_rank()) {
      return this->aux_rank() < that_cast.aux_rank();
    }

    //          v1: compare hashes only
    //          return Expr::static_less_than(that);
    //          v2: compare fully
    if (this->bra_hash_value() != that_cast.bra_hash_value()) {
      return std::lexicographical_compare(
          this->bra().begin(), this->bra().end(), that_cast.bra().begin(),
          that_cast.bra().end());
    }

    if (this->ket() != that_cast.ket()) {
      return std::lexicographical_compare(
          this->ket().begin(), this->ket().end(), that_cast.ket().begin(),
          that_cast.ket().end());
    }

    return std::lexicographical_compare(this->aux().begin(), this->aux().end(),
                                        that_cast.aux().begin(),
                                        that_cast.aux().end());
  }

  // these implement the AbstractTensor interface
  AbstractTensor::const_any_view_randsz _bra() const override final {
    return ranges::counted_view<const Index *>(
        bra_.empty() ? nullptr : &(bra_[0]), bra_.size());
  }
  AbstractTensor::const_any_view_randsz _ket() const override final {
    return ranges::counted_view<const Index *>(
        ket_.empty() ? nullptr : &(ket_[0]), ket_.size());
  }
  AbstractTensor::const_any_view_randsz _aux() const override final {
    return ranges::counted_view<const Index *>(
        aux_.empty() ? nullptr : &(aux_[0]), aux_.size());
  }
  AbstractTensor::const_any_view_rand _braket() const override final {
    return braket();
  }
  AbstractTensor::const_any_view_rand _indices() const override final {
    return indices();
  }
  std::size_t _bra_rank() const override final { return bra_rank(); }
  std::size_t _ket_rank() const override final { return ket_rank(); }
  std::size_t _aux_rank() const override final { return aux_rank(); }
  Symmetry _symmetry() const override final { return symmetry_; }
  BraKetSymmetry _braket_symmetry() const override final {
    return braket_symmetry_;
  }
  ParticleSymmetry _particle_symmetry() const override final {
    return particle_symmetry_;
  }
  std::size_t _color() const override final { return 0; }
  bool _is_cnumber() const override final { return true; }
  std::wstring_view _label() const override final { return label_; }
  std::wstring _to_latex() const override final { return to_latex(); }
  bool _transform_indices(
      const container::map<Index, Index> &index_map) override final {
    return transform_indices(index_map);
  }
  void _reset_tags() override final { reset_tags(); }
  bool operator<(const AbstractTensor &other) const override final {
    auto *other_tensor = dynamic_cast<const Tensor *>(&other);
    if (other_tensor) {
      const Expr *other_expr = static_cast<const Expr *>(other_tensor);
      return this->static_less_than(*other_expr);
    } else
      return false;  // TODO do we compare typeid? labels? probably the latter
  }

  AbstractTensor::any_view_randsz _bra_mutable() override final {
    this->reset_hash_value();
    return ranges::counted_view<Index *>(bra_.empty() ? nullptr : &(bra_[0]),
                                         bra_.size());
  }
  AbstractTensor::any_view_randsz _ket_mutable() override final {
    this->reset_hash_value();
    return ranges::counted_view<Index *>(ket_.empty() ? nullptr : &(ket_[0]),
                                         ket_.size());
  }
  AbstractTensor::any_view_randsz _aux_mutable() override final {
    this->reset_hash_value();
    return ranges::counted_view<Index *>(aux_.empty() ? nullptr : &(aux_[0]),
                                         aux_.size());
  }

};  // class Tensor

static_assert(is_tensor_v<Tensor>,
              "The Tensor class does not fulfill the requirements of the "
              "Tensor interface");

using TensorPtr = std::shared_ptr<Tensor>;

/// make_overlap tensor label is reserved since it is used by low-level SeQuant
/// machinery. Users can create make_overlap Tensor using make_overlap()
inline std::wstring overlap_label() { return L"s"; }

inline ExprPtr make_overlap(const Index &bra_index, const Index &ket_index) {
  return ex<Tensor>(Tensor(overlap_label(), bra{bra_index}, ket{ket_index},
                           aux{}, Tensor::reserved_tag{}));
}

}  // namespace sequant

#endif  // SEQUANT_TENSOR_HPP
