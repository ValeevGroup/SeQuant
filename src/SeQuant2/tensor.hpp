//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_TENSOR_HPP
#define SEQUANT2_TENSOR_HPP

#include <memory>

#include "index.hpp"
#include "expr.hpp"

namespace sequant2 {

enum class Symmetry { symm, antisymm, nonsymm };

class TensorCanonicalizer;

/// @brief particle-symmetric Tensor, i.e. permuting
class Tensor : public Expr {
 private:
  auto make_indices(WstrList index_labels) {
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

  Tensor(std::wstring_view label,
         IndexList bra_indices,
         IndexList ket_indices,
         Symmetry s = Symmetry::nonsymm)
      : label_(label), bra_(bra_indices), ket_(ket_indices), symmetry_(s) {}

  Tensor(std::wstring_view label,
         WstrList bra_index_labels,
         WstrList ket_index_labels,
         Symmetry s = Symmetry::nonsymm)
      : label_(label), bra_(make_indices(bra_index_labels)), ket_(make_indices(ket_index_labels)), symmetry_(s) {}

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

  std::shared_ptr<Expr> canonicalize() override;

  Tensor &transform_indices(const std::map<Index, Index> &index_map) {
    bool mutated = false;
    ranges::for_each(braket(), [index_map, &mutated](auto &idx) {
      auto it = index_map.find(idx);
      if (it != index_map.end()) {
        idx = it->second;
        mutated = true;
      }
    });
    if (mutated)
      this->reset_hash_value();
    return *this;
  }

  type_id_type type_id() const override {
    return get_type_id<Tensor>();
  };

  std::shared_ptr<Expr> clone() const override {
    return make<Tensor>(*this);
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
    if (this->bra_rank() == that_cast.bra_rank()) {
      if (this->ket_rank() == that_cast.ket_rank()) {
        if (this->label() == that_cast.label()) {
          return Expr::static_less_than(that);
        } else {
          return this->label() < that_cast.label();
        }
      } else {
        return this->ket_rank() < that_cast.ket_rank();
      }
    } else {
      return this->bra_rank() < that_cast.bra_rank();
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
  static void register_instance(std::shared_ptr<TensorCanonicalizer> canonicalizer, std::wstring_view label = L"");

  auto &bra(Tensor &t) { return t.bra_; };
  auto &ket(Tensor &t) { return t.ket_; };
  auto braket(Tensor &t) { return ranges::view::concat(bra(t), ket(t)); }

  virtual std::shared_ptr<Expr> apply(Tensor &) = 0;

 private:
  static std::map<std::wstring, std::shared_ptr<TensorCanonicalizer>> &instance_map_accessor();
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
    assert(is_antisymm);

    //
    auto sort_swappables = [](const auto &begin, const auto &end, const auto &compare) {
      const auto len = end - begin;
      switch (len) {
        case 0: return;

        case 1: return;

        case 2: {
          auto &elem1 = *begin;
          auto &elem2 = *(begin + 1);
          if (compare(elem1, elem2))
            return;
          else {
            swap(elem1, elem2);
            return;
          }
        }
          break;

        case 3: {
          // bubble sort
          // {2,3} -> [2,3]
          {
            auto &elem2 = *(begin + 1);
            auto &elem3 = *(begin + 2);
            const auto lt23 = compare(elem2, elem3);
            if (!lt23)
              swap(elem2, elem3);
          }
          // sort {1,[2,3]} -> [1, 2, 3]
          {
            auto &elem1 = *(begin);
            auto &elem2 = *(begin + 1);
            auto &elem3 = *(begin + 2);
            const auto lt12 = compare(elem1, elem2);
            const auto lt13 = compare(elem1, elem3);
            if (!lt12)
              swap(elem1, elem2);
            if (!lt13)
              swap(elem2, elem3);
          }
          return;
        }
          break;

        default:abort();  // not yet implemented
      }
    };

    bool even = true;
    switch (symmetry) {
      case Symmetry::antisymm: {
        auto &_bra = this->bra(t);
        auto &_ket = this->ket(t);
        using std::begin;
        using std::end;
//      std::wcout << "canonicalizing " << to_latex(t);
        IndexSwapper::thread_instance().reset();
        // std::{stable_}sort does not necessarily use swap! so must implement sort outselves .. thankfully ranks will be low so can stick with bubble
        sort_swappables(begin(_bra), end(_bra), comp);
        sort_swappables(begin(_ket), end(_ket), comp);
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

    std::shared_ptr<Expr> result = is_antisymm ? (even == false ? make<Constant>(-1) : nullptr) : nullptr;
    return result;
  }

 private:
  std::map<std::wstring, Index> external_indices_;
};

inline std::shared_ptr<Expr> overlap(const Index& bra_index, const Index& ket_index) {
  return std::make_shared<Tensor>(L"S", IndexList{bra_index}, IndexList{ket_index});
}

}  // namespace sequant2

#endif //SEQUANT2_TENSOR_HPP
