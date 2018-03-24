//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_TENSOR_HPP
#define SEQUANT2_TENSOR_HPP

#include "expr.hpp"

namespace sequant2 {

enum class Symmetry { symm, antisymm, nonsymm };

/// @brief particle-symmetric Tensor, i.e. permuting
class Tensor : public Expr {
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
  const std::vector<Index> &bra() const { return bra_; }
  const std::vector<Index> &ket() const { return ket_; }
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

 private:
  std::wstring label_{};
  std::vector<Index> bra_{};
  std::vector<Index> ket_{};
  Symmetry symmetry_ = Symmetry::nonsymm;

  std::vector<Index> make_indices(std::initializer_list<std::wstring_view> index_labels) {
    std::vector<Index> result;
    result.reserve(index_labels.size());
    for (const auto &label: index_labels) {
      result.push_back(Index{label});
    }
    return result;
  }
};

inline std::wstring to_latex(const Tensor &t) {
  std::wstring result;
  result = L"{";
  result += t.label();
  result += L"^{";
  for (const auto &i : t.ket())
    result += to_latex(i);
  result += L"}_{";
  for (const auto &i : t.bra())
    result += to_latex(i);
  result += L"}}";
  return result;
}

}  // namespace sequant2

#endif //SEQUANT2_TENSOR_HPP
