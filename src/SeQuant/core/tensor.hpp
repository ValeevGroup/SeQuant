//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT_TENSOR_HPP
#define SEQUANT_TENSOR_HPP

#include <memory>

#include "abstract_tensor.hpp"
#include "algorithm.hpp"
#include "attr.hpp"
#include "expr.hpp"
#include "index.hpp"
#include "sequant.hpp"

namespace sequant {

/// @brief particle-symmetric Tensor, i.e. permuting
class Tensor : public Expr, public AbstractTensor {
 private:
  using index_container_type = container::svector<Index, 4>;
  auto make_indices(IndexList indices) { return indices; }
  auto make_indices(WstrList index_labels) {
    index_container_type result;
    result.reserve(index_labels.size());
    for (const auto &label: index_labels) {
      result.push_back(Index{label});
    }
    return result;
  }
  auto make_indices(std::initializer_list<const wchar_t*> index_labels) {
    index_container_type result;
    result.reserve(index_labels.size());
    for (const auto &label: index_labels) {
      result.push_back(Index{label});
    }
    return result;
  }
  template <typename IndexContainer>
  index_container_type make_indices(IndexContainer &&indices) {
    if constexpr (std::is_same_v<index_container_type,
                                 std::decay_t<IndexContainer>>) {
      return std::forward<IndexContainer>(indices);
    } else {
      return index_container_type(std::begin(indices), std::end(indices));
    }
  }

  /// @return view of the bra+ket index ranges
  auto braket() { return ranges::views::concat(bra_, ket_); }

 public:
  Tensor() = default;
  virtual ~Tensor();

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param symmetry the symmetry of bra or ket
  /// @param braket_symmetry the symmetry with respect to bra-ket exchange
  template <typename IndexContainer,
            typename = std::enable_if_t<
                !meta::is_initializer_list_v<std::decay_t<IndexContainer>>>>
  Tensor(std::wstring_view label, IndexContainer &&bra_indices,
         IndexContainer &&ket_indices, Symmetry s = Symmetry::nonsymm, BraKetSymmetry bks = get_default_context().braket_symmetry())
      : label_(label),
        bra_(make_indices(bra_indices)),
        ket_(make_indices(ket_indices)),
        symmetry_(s),
        braket_symmetry_(bks) {}

  /// @tparam I1 any type convertible to Index)
  /// @tparam I2 any type convertible to Index
  /// @note  I1 and I2 default to Index to allow empty lists
  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted to indices)
  /// @param symmetry the symmetry of bra or ket
  /// @param braket_symmetry the symmetry with respect to bra-ket exchange
  template <typename I1 = Index, typename I2 = Index>
  Tensor(std::wstring_view label,
         std::initializer_list<I1> bra_indices,
         std::initializer_list<I2> ket_indices,
         Symmetry s = Symmetry::nonsymm, BraKetSymmetry bks = get_default_context().braket_symmetry())
      : label_(label), bra_(make_indices(bra_indices)), ket_(make_indices(ket_indices)), symmetry_(s), braket_symmetry_(bks) {}

  std::wstring_view label() const { return label_; }
  const auto& bra() const { return bra_; }
  const auto& ket() const { return ket_; }
  /// @return joined view of the bra and ket index ranges
  auto braket() const { return ranges::views::concat(bra_, ket_); }
  /// @return view of the bra+ket index ranges
  /// @note this is to work around broken lookup rules
  auto const_braket() const { return this->braket(); }
  /// Returns the Symmetry object describing the symmetry of the bra and ket of the Tensor, i.e. what effect swapping indices in positions @c i  and @c j in <b>either bra or ket</em> has on the elements of the Tensor;
  /// Tensor's are <em>always assumed</em> to be particle-symmetric, i.e. swapping indices in positions @c i and @c j in <b>both bra and ket</b>;
  /// The allowed values are Symmetry::symm, Symmetry::antisymm, and Symmetry::nonsymm
  /// @return the Symmetry object describing the symmetry of the bra and ket of the Tensor.
  Symmetry symmetry() const { return symmetry_; }
  /// @return the BraKetSymmetry object describing the symmetry of the Tensor under exchange of bra and ket.
  BraKetSymmetry braket_symmetry() const { return braket_symmetry_; }

  /// @return number of bra indices
  auto bra_rank() const { return bra_.size(); }
  /// @return number of ket indices
  auto ket_rank() const { return ket_.size(); }
  /// @return number of indices in bra/ket
  /// @throw std::logic_error if bra and ket ranks do not match
  auto rank() const {
    if (bra_rank() != ket_rank()) {
      throw std::logic_error("Tensor::rank(): bra rank != ket rank");
    }
    return bra_rank();
  }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    result += this->label();
    result += L"^{";
    for (const auto &i : this->ket())
      result += sequant::to_latex(i);
    result += L"}_{";
    for (const auto &i : this->bra())
      result += sequant::to_latex(i);
    result += L"}}";
    return result;
  }

  std::wstring to_wolfram() const override {
    std::wstring result;
    result = L"SQM[OHead[\"\\!\\(\\*OverscriptBox[\\(";
    result += this->label();
    result += L"\\), \\(_\\)]\\)\",";
    result += sequant::to_wolfram(this->symmetry());
    result += L"],";
    for (const auto &i : this->ket()) {
      result += i.to_wolfram(BraKetPos::ket) + L",";
    }
    for (const auto &i : this->bra()) {
      result += i.to_wolfram(BraKetPos::bra) + L",";
    }
    result = result.erase(result.size() - 1);
    result += L"]";
    return result;
  }

  std::shared_ptr<Expr> canonicalize() override;

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
    if (mutated)
      this->reset_hash_value();
    return mutated;
  }

  type_id_type type_id() const override {
    return get_type_id<Tensor>();
  };

  std::shared_ptr<Expr> clone() const override {
    return ex<Tensor>(*this);
  }

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
  index_container_type bra_{};
  index_container_type ket_{};
  Symmetry symmetry_ = Symmetry::invalid;
  BraKetSymmetry braket_symmetry_ = BraKetSymmetry::invalid;
  mutable std::optional<hash_type>
      bra_hash_value_;  // memoized byproduct of memoizing_hash()

  hash_type memoizing_hash() const override {
    using std::begin;
    using std::end;
    auto val = boost::hash_range(begin(bra()), end(bra()));
    bra_hash_value_ = val;
    boost::hash_range(val, begin(ket()), end(ket()));
    boost::hash_combine(val, label_);
    boost::hash_combine(val, symmetry_);
    hash_value_ = val;
    return *hash_value_;
  }
  void reset_hash_value() const override {
    Expr::reset_hash_value();
    bra_hash_value_.reset();
  }

  bool static_equal(const Expr &that) const override {
    const auto& that_cast = static_cast<const Tensor&>(that);
    if (this->label() == that_cast.label() && this->symmetry() == that_cast.symmetry() && this->bra_rank() == that_cast.bra_rank() && this->ket_rank() == that_cast.ket_rank()) {
      // compare hash values first
      if (this->hash_value() == that.hash_value()) // hash values agree -> do full comparison
        return this->bra() == that_cast.bra() && this->ket() == that_cast.ket();
      else
        return false;
    } else return false;
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
    }
    else {
      return this->label() < that_cast.label();
    }
  }

  // these implement the AbstractTensor interface
  AbstractTensor::const_any_view_randsz _bra() const override final {
    return ranges::counted_view<const Index*>(bra_.empty() ? nullptr : &(bra_[0]), bra_.size());
  }
  AbstractTensor::const_any_view_randsz _ket() const override final {
    return ranges::counted_view<const Index*>(ket_.empty() ? nullptr : &(ket_[0]), ket_.size());
  }
  AbstractTensor::const_any_view_rand _braket() const override final {
    return braket();
  }
  std::size_t _bra_rank() const override final {
    return bra_rank();
  }
  std::size_t _ket_rank() const override final {
    return ket_rank();
  }
  Symmetry _symmetry() const override final {
    return symmetry_;
  }
  BraKetSymmetry _braket_symmetry() const override final {
    return braket_symmetry_;
  }
  std::size_t _color() const override final {
    return 0;
  }
  bool _is_cnumber() const override final {
    return true;
  }
  std::wstring _label() const override final {
    return label_;
  }
  std::wstring _to_latex() const override final {
    return to_latex();
  }
  bool _transform_indices(const container::map<Index, Index>& index_map) override final {
    return transform_indices(index_map);
  }
  void _reset_tags() override final {
    reset_tags();
  }
  bool operator<(const AbstractTensor& other) const override final {
    auto* other_tensor = dynamic_cast<const Tensor*>(&other);
    if (other_tensor) {
      const Expr* other_expr = static_cast<const Expr*>(other_tensor);
      return this->static_less_than(*other_expr);
    }
    else
      return false; // TODO do we compare typeid? labels? probably the latter
  }

  AbstractTensor::any_view_randsz _bra_mutable() override final {
    this->reset_hash_value();
    return ranges::counted_view<Index*>(bra_.empty() ? nullptr : &(bra_[0]), bra_.size());
  }
  AbstractTensor::any_view_randsz _ket_mutable() override final {
    this->reset_hash_value();
    return ranges::counted_view<Index*>(ket_.empty() ? nullptr : &(ket_[0]), ket_.size());
  }

};  // class Tensor

using TensorPtr = std::shared_ptr<Tensor>;

inline std::shared_ptr<Expr> overlap(const Index& bra_index, const Index& ket_index) {
  return ex<Tensor>(L"S", IndexList{bra_index}, IndexList{ket_index});
}

}  // namespace sequant

#endif //SEQUANT_TENSOR_HPP
