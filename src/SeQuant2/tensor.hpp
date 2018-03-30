//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_TENSOR_HPP
#define SEQUANT2_TENSOR_HPP

#include <memory>

#include "expr.hpp"

namespace sequant2 {

enum class Symmetry { symm, antisymm, nonsymm };

/// @brief particle-symmetric Tensor, i.e. permuting
class Tensor : public Expr {
 private:
  auto make_indices(std::initializer_list<std::wstring_view> index_labels) {
    index_container_type result;
    result.reserve(index_labels.size());
    for (const auto &label: index_labels) {
      result.push_back(Index{label});
    }
    return result;
  }

 public:
  Tensor() = default;
  virtual ~Tensor() = default;

  Tensor(std::wstring_view label,
         std::initializer_list<Index> bra_indices,
         std::initializer_list<Index> ket_indices,
         Symmetry s = Symmetry::nonsymm)
      : label_(label), bra_(bra_indices), ket_(ket_indices), symmetry_(s) {}

  Tensor(std::wstring_view label,
         std::initializer_list<std::wstring_view> bra_index_labels,
         std::initializer_list<std::wstring_view> ket_index_labels,
         Symmetry s = Symmetry::nonsymm)
      : label_(label), bra_(make_indices(bra_index_labels)), ket_(make_indices(ket_index_labels)) {}

  std::wstring_view label() const { return label_; }
  const auto& bra() const { return bra_; }
  const auto& ket() const { return ket_; }
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

  std::shared_ptr<Expr> canonicalize() override {
    abort();
  }

  type_id_type type_id() const override {
    return get_type_id<Tensor>();
  };

 private:
  std::wstring label_{};
  using index_container_type = container::svector<Index>;
  index_container_type bra_{};
  index_container_type ket_{};
  Symmetry symmetry_ = Symmetry::nonsymm;

  bool static_compare(const Expr& that) const override {
    const auto& that_cast = static_cast<const Tensor&>(that);
    if (this->label() == that_cast.label() && this->symmetry() == that_cast.symmetry() && this->bra_rank() == that_cast.bra_rank() && this->ket_rank() == that_cast.ket_rank()) {
      // compare hash values first
      if (this->hash_value() == that.hash_value()) // hash values agree -> do full comparison
        return this->bra() == that_cast.bra() && this->ket() == that_cast.ket();
      else
        return false;
    } else return false;
  }

};

inline std::shared_ptr<Expr> overlap(const Index& bra_index, const Index& ket_index) {
  return std::make_shared<Tensor>(L"S", std::initializer_list<Index>{bra_index}, std::initializer_list<Index>{ket_index});
}

}  // namespace sequant2

#endif //SEQUANT2_TENSOR_HPP
