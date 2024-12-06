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

/// @brief particle-symmetric Tensor, i.e. permuting
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

  /// @return view of the bra+ket index ranges
  auto braket() { return ranges::views::concat(bra_, ket_); }

  /// asserts that @p label is not reserved
  /// @note Tensor with reserved labels are constructed using friends of Tensor
  /// @param label a Tensor label candidate
  void assert_nonreserved_label(std::wstring_view label) const;

  // utility for dispatching to private ctor
  struct reserved_tag {};

  // list of friends who can make Tensor objects with reserved labels
  friend ExprPtr make_overlap(const Index &bra_index, const Index &ket_index);

  template <
      typename IndexRange1, typename IndexRange2,
      typename = std::enable_if_t<
          (meta::is_statically_castable_v<
              meta::range_value_t<IndexRange1>,
              Index>)&&(meta::
                            is_statically_castable_v<
                                meta::range_value_t<IndexRange2>, Index>)>>
  Tensor(std::wstring_view label, const bra<IndexRange1> &bra_indices,
         const ket<IndexRange2> &ket_indices, reserved_tag,
         Symmetry s = Symmetry::nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ParticleSymmetry ps = ParticleSymmetry::symm)
      : label_(label),
        bra_(make_indices(bra_indices)),
        ket_(make_indices(ket_indices)),
        symmetry_(s),
        braket_symmetry_(bks),
        particle_symmetry_(ps) {
    validate_symmetries();
  }

  Tensor(std::wstring_view label, bra<index_container_type> &&bra_indices,
         ket<index_container_type> &&ket_indices, reserved_tag,
         Symmetry s = Symmetry::nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ParticleSymmetry ps = ParticleSymmetry::symm)
      : label_(label),
        bra_(std::move(bra_indices)),
        ket_(std::move(ket_indices)),
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
      : Tensor(label, bra_indices, ket_indices, reserved_tag{}, s, bks, ps) {
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
               reserved_tag{}, s, bks, ps) {
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
  const auto &bra() const { return bra_; }
  const auto &ket() const { return ket_; }
  /// @return joined view of the bra and ket index ranges
  auto braket() const { return ranges::views::concat(bra_, ket_); }
  /// @return view of the bra+ket index ranges
  /// @note this is to work around broken lookup rules
  auto const_braket() const { return this->braket(); }
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
  /// @return number of indices in bra/ket
  /// @throw std::logic_error if bra and ket ranks do not match
  std::size_t rank() const {
    if (bra_rank() != ket_rank()) {
      throw std::logic_error("Tensor::rank(): bra rank != ket rank");
    }
    return bra_rank();
  }

  std::wstring to_latex() const override {
    // either rank > 1 or sum of bra and ket ranks > 1
    const bool add_bar =
        bra_rank() == ket_rank() ? rank() > 1 : bra_rank() + ket_rank() > 1;

    std::wstring result = L"{";
    if ((this->symmetry() == Symmetry::antisymm) && add_bar)
      result += L"\\bar{";
    result += utf_to_latex(this->label());
    if ((this->symmetry() == Symmetry::antisymm) && add_bar) result += L"}";
    result += L"^{";
    for (const auto &i : this->ket()) result += sequant::to_latex(i);
    result += L"}_{";
    for (const auto &i : this->bra()) result += sequant::to_latex(i);
    result += L"}}";
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
    ranges::for_each(braket(), [&](auto &idx) {
      if (idx.transform(index_map)) mutated = true;
    });
    if (mutated) this->reset_hash_value();
    return mutated;
  }

  type_id_type type_id() const override { return get_type_id<Tensor>(); };

  ExprPtr clone() const override { return ex<Tensor>(*this); }

  void reset_tags() const {
    ranges::for_each(braket(), [](const auto &idx) { idx.reset_tag(); });
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
        this->ket_rank() == that_cast.ket_rank()) {
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
    const auto &that_cast = static_cast<const Tensor &>(that);
    if (this == &that) return false;
    if (this->label() == that_cast.label()) {
      if (this->bra_rank() == that_cast.bra_rank()) {
        if (this->ket_rank() == that_cast.ket_rank()) {
          //          v1: compare hashes only
          //          return Expr::static_less_than(that);
          //          v2: compare fully
          if (this->bra_hash_value() == that_cast.bra_hash_value()) {
            return std::lexicographical_compare(
                this->ket().begin(), this->ket().end(), that_cast.ket().begin(),
                that_cast.ket().end());
          } else {
            return std::lexicographical_compare(
                this->bra().begin(), this->bra().end(), that_cast.bra().begin(),
                that_cast.bra().end());
          }
        } else {
          return this->ket_rank() < that_cast.ket_rank();
        }
      } else {
        return this->bra_rank() < that_cast.bra_rank();
      }
    } else {
      return this->label() < that_cast.label();
    }
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
  AbstractTensor::const_any_view_rand _braket() const override final {
    return braket();
  }
  std::size_t _bra_rank() const override final { return bra_rank(); }
  std::size_t _ket_rank() const override final { return ket_rank(); }
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

};  // class Tensor

using TensorPtr = std::shared_ptr<Tensor>;

/// make_overlap tensor label is reserved since it is used by low-level SeQuant
/// machinery. Users can create make_overlap Tensor using make_overlap()
inline std::wstring overlap_label() { return L"s"; }

inline ExprPtr make_overlap(const Index &bra_index, const Index &ket_index) {
  return ex<Tensor>(Tensor(overlap_label(), bra{bra_index}, ket{ket_index},
                           Tensor::reserved_tag{}));
}

}  // namespace sequant

#endif  // SEQUANT_TENSOR_HPP
