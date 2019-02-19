//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_TENSOR_HPP
#define SEQUANT2_TENSOR_HPP

#include <memory>

#include "algorithm.hpp"
#include "attr.hpp"
#include "expr.hpp"
#include "index.hpp"

namespace sequant2 {

class TensorCanonicalizer;

/// @brief particle-symmetric Tensor, i.e. permuting
class Tensor : public Expr {
 private:
  auto make_indices(IndexList indices) {
    return indices;
  }
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

  /// @return view of the bra+ket index ranges
  auto braket() { return ranges::view::concat(bra_, ket_); }

 public:
  Tensor() = default;
  virtual ~Tensor();

  /// @tparam I1 any type convertible to Index)
  /// @tparam I2 any type convertible to Index
  /// @note  I1 and I2 default to Index to allow empty lists
  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted to indices)
  /// @param symmetry the tensor symmetry
  template <typename I1 = Index, typename I2 = Index>
  Tensor(std::wstring_view label,
         std::initializer_list<I1> bra_indices,
         std::initializer_list<I2> ket_indices,
         Symmetry s = Symmetry::nonsymm)
      : label_(label), bra_(make_indices(bra_indices)), ket_(make_indices(ket_indices)), symmetry_(s) {}

  std::wstring_view label() const { return label_; }
  const auto& bra() const { return bra_; }
  const auto& ket() const { return ket_; }
  /// @return view of the bra+ket index ranges
  auto braket() const { return ranges::view::concat(bra_, ket_); }
  Symmetry symmetry() const { return symmetry_; }

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
      result += sequant2::to_latex(i);
    result += L"}_{";
    for (const auto &i : this->bra())
      result += sequant2::to_latex(i);
    result += L"}}";
    return result;
  }

  std::wstring to_wolfram() const override {
    std::wstring result;
    result = L"SQM[OHead[\"\\!\\(\\*OverscriptBox[\\(";
    result += this->label();
    result += L"\\), \\(_\\)]\\)\",";
    result += sequant2::to_wolfram(this->symmetry());
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

  /// Replaced indices using the index map
  /// @return true if one or more indices changed
  template <template <typename, typename, typename... Args> class Map,
            typename... Args>
  bool transform_indices(const Map<Index, Index, Args...> &index_map) {
    bool mutated = false;
    ranges::for_each(braket(), [index_map, &mutated](auto &idx) {
      if (idx.transform(index_map))
        mutated = true;
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
    ranges::for_each(braket(), [](const auto &idx) {
      idx.reset_tag();
    });
  }

 private:
  std::wstring label_{};
  using index_container_type = container::svector<Index>;
  index_container_type bra_{};
  index_container_type ket_{};
  Symmetry symmetry_ = Symmetry::nonsymm;

  hash_type memoizing_hash() const override {
    using std::begin;
    using std::end;
    auto val = boost::hash_range(begin(braket()), end(braket()));
    boost::hash_combine(val, label_);
    boost::hash_combine(val, symmetry_);
    hash_value_ = val;
    return *hash_value_;
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
    if (this->label() == that_cast.label()) {
      if (this->bra_rank() == that_cast.bra_rank()) {
        if (this->ket_rank() == that_cast.ket_rank()) {
          return Expr::static_less_than(that);
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

  friend class TensorCanonicalizer;
};  // class Tensor

using TensorPtr = std::shared_ptr<Tensor>;

/// @brief Base class for Tensor canonicalizers
/// To make custom canonicalizer make a derived class and register an instance of that class with TensorCanonicalizer::register_instance
class TensorCanonicalizer {
 public:
  virtual ~TensorCanonicalizer();
  /// returns a TensorCanonicalizer previously registered via TensorCanonicalizer::register_instance()
  /// with @c label
  static std::shared_ptr<TensorCanonicalizer> instance(std::wstring_view label = L"");
  /// registers @c canonicalizer to be applied to Tensor objects with label @c label ; leave the label
  /// empty if @c canonicalizer is to apply to Tensor with any label)
  static void register_instance(
      std::shared_ptr<TensorCanonicalizer> canonicalizer,
      std::wstring_view label = L"");

  /// @return a list of Tensor labels with lexicographic preference (in order)
  static const auto &cardinal_tensor_labels() {
    return cardinal_tensor_labels_accessor();
  }
  /// @param cardinal_tensor_labels a list of Tensor labels with lexicographic
  /// preference (in order)
  static void set_cardinal_tensor_labels(
      const container::vector<std::wstring> &labels) {
    cardinal_tensor_labels_accessor() = labels;
  }

  auto &bra(Tensor &t) { return t.bra_; };
  auto &ket(Tensor &t) { return t.ket_; };
  auto braket(Tensor &t) { return ranges::view::concat(bra(t), ket(t)); }

  virtual std::shared_ptr<Expr> apply(Tensor &) = 0;

 private:
  static container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>
      &instance_map_accessor();
  static container::vector<std::wstring> &cardinal_tensor_labels_accessor();
};

class DefaultTensorCanonicalizer : public TensorCanonicalizer {
 public:
  DefaultTensorCanonicalizer() = default;

  /// @tparam IndexContainer a Container of Index objects such that @c IndexContainer::value_type is convertible to Index (e.g. this can be std::vector or std::set , but not std::map)
  /// @param external_indices container of external Index objects
  /// @warning @c external_indices is assumed to be immutable during the lifetime of this object
  template<typename IndexContainer>
  DefaultTensorCanonicalizer(IndexContainer &&external_indices) {
    ranges::for_each(external_indices, [this](const Index &idx) {
      this->external_indices_.emplace(std::make_pair(std::wstring(idx.label()), idx));
    });
  }
  virtual ~DefaultTensorCanonicalizer() = default;

  /// Implements TensorCanonicalizer::apply
  /// @note Canonicalizes @c t by sorting its bra (if @c t.symmetry()==Symmetry::nonsymm ) or its bra and ket (if @c t.symmetry()!=Symmetry::nonsymm ),
  ///       with the external indices
  std::shared_ptr<Expr> apply(Tensor &t) override;

  template<typename Compare>
  std::shared_ptr<Expr> apply(Tensor &t, const Compare &comp) {
    auto symmetry = t.symmetry();
    auto is_antisymm = symmetry == Symmetry::antisymm;

    // can only handle anisymmetric case so far
#ifndef NDEBUG
    if (t.bra_rank() > 1 || t.ket_rank() > 1)
      assert(is_antisymm);
#endif

    bool even = true;
    switch (symmetry) {
      case Symmetry::antisymm: {
        auto &_bra = this->bra(t);
        auto &_ket = this->ket(t);
        using std::begin;
        using std::end;
//      std::wcout << "canonicalizing " << to_latex(t);
        IndexSwapper::thread_instance().reset();
        // std::{stable_}sort does not necessarily use swap! so must implement
        // sort outselves .. thankfully ranks will be low so can stick with
        // bubble
        bubble_sort(begin(_bra), end(_bra), comp);
        bubble_sort(begin(_ket), end(_ket), comp);
        even = IndexSwapper::thread_instance().even_num_of_swaps();
//      std::wcout << " is " << (even ? "even" : "odd") << " and produces " << to_latex(t) << std::endl;
      }
        break;

      case Symmetry::symm: {

      }
        break;

      case Symmetry::nonsymm: {

      }
        break;

      default:abort();
    }

    std::shared_ptr<Expr> result = is_antisymm ? (even == false ? ex<Constant>(-1) : nullptr) : nullptr;
    return result;
  }

 private:
  container::map<std::wstring, Index> external_indices_;
};

inline std::shared_ptr<Expr> overlap(const Index& bra_index, const Index& ket_index) {
  return ex<Tensor>(L"S", IndexList{bra_index}, IndexList{ket_index});
}

}  // namespace sequant2

#endif //SEQUANT2_TENSOR_HPP
