#ifndef SEQUANT_UTILITY_TENSOR_HPP
#define SEQUANT_UTILITY_TENSOR_HPP

#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <functional>
#include <vector>

namespace sequant {

/// Comparator template that compares tensor blocks (slots). This means that it
/// takes only the spaces of indices into account but not their concrete
/// labeling (i.e. i1 and i2 are considered equal).
/// Note that the comparison will NOT take the exact location in bra/ket/aux
/// into consideration explicitly. It compares the sequence of all indices which
/// will lead to t{a1,i1} comparing equal to e.g. t{;a1,i1;}. The assumption
/// here is that the bra/ket/aux is not important at the level where tensor
/// blocks are relevant (e.g. during numerical evaluation).
template <template <class> class Comparator, template <class> class Selector,
          bool fallback>
struct TensorBlockComparator {
  auto operator()(const Tensor &lhs, const Tensor &rhs) const {
    if (lhs.label() != rhs.label()) {
      Comparator<decltype(lhs.label())> cmp;
      return cmp(lhs.label(), rhs.label());
    }

    if (lhs.num_slots() != rhs.num_slots()) {
      Comparator<decltype(lhs.num_slots())> cmp;
      return cmp(lhs.num_slots(), rhs.num_slots());
    }

    if (lhs.num_indices() != rhs.num_indices()) {
      Comparator<decltype(lhs.num_indices())> cmp;
      return cmp(lhs.num_indices(), rhs.num_indices());
    }

    auto &&lhs_indices = lhs.indices();
    auto &&rhs_indices = rhs.indices();

    Comparator<IndexSpace> cmp;
    Selector<IndexSpace> selector;

    for (auto lhs_it = lhs_indices.begin(), rhs_it = rhs_indices.begin();
         lhs_it != lhs_indices.end(); ++lhs_it, ++rhs_it) {
      SEQUANT_ASSERT(rhs_it != rhs_indices.end());

      const IndexSpace &left = lhs_it->space();
      const IndexSpace &right = rhs_it->space();

      if (selector(left, right)) {
        return cmp(left, right);
      }
    }

    return fallback;
  }

  bool operator()(const Expr &lhs, const Expr &rhs) const {
    if (lhs.is<Tensor>() && rhs.is<Tensor>()) {
      return (*this)(lhs.as<Tensor>(), rhs.as<Tensor>());
    }

    Comparator cmp;

    return cmp(lhs, rhs);
  }

  bool operator()(const ExprPtr &lhs, const ExprPtr &rhs) const {
    return (*this)(*lhs, *rhs);
  }
};

/// Compares tensor blocks (slots) for equality
using TensorBlockEqualComparator =
    TensorBlockComparator<std::equal_to, std::not_equal_to, true>;
/// Compares tensor blocks (slots) on a less-than relationship
using TensorBlockLessThanComparator =
    TensorBlockComparator<std::less, std::not_equal_to, false>;

/// Similar to TensorBlockComparator but in case of tensors that compare equal
/// based on their tensor blocks (slots), they are compared on if and where they
/// have specific indices (compared the usual way, including ordinals).
template <template <class> class Comparator, template <class> class Selector,
          bool fallback>
struct IndexSpecificTensorBlockComparator {
  IndexSpecificTensorBlockComparator() = default;
  IndexSpecificTensorBlockComparator(std::vector<Index> indices)
      : indices_(std::move(indices)) {}

  const std::vector<Index> &indices() const { return indices_; }
  void set_indices(std::vector<Index> indices) {
    indices_ = std::move(indices);
  }

  auto operator()(const Tensor &lhs, const Tensor &rhs) const {
    bool equal = TensorBlockEqualComparator{}(lhs, rhs);
    if (!equal) {
      return TensorBlockComparator<Comparator, Selector, fallback>{}(lhs, rhs);
    }

    // Tensor blocks are equal, now check for specific indices
    auto &&lhs_indices = lhs.indices();
    auto &&rhs_indices = rhs.indices();

    Selector<Index> idx_selector;

    SEQUANT_ASSERT(lhs.num_indices() == rhs.num_indices());
    auto lit = lhs_indices.begin();
    auto rit = rhs_indices.begin();
    while (lit != lhs_indices.end()) {
      if (idx_selector(*lit, *rit)) {
        const std::size_t lhs_pos = std::ranges::distance(
            indices_.begin(), std::ranges::find(indices_, *lit));
        const std::size_t rhs_pos = std::ranges::distance(
            indices_.begin(), std::ranges::find(indices_, *rit));

        Selector<std::size_t> pos_select;
        if (pos_select(lhs_pos, rhs_pos)) {
          Comparator<std::size_t> pos_cmp;
          return pos_cmp(lhs_pos, rhs_pos);
        }
      }

      ++lit;
      ++rit;
    }

    return fallback;
  }

  bool operator()(const Expr &lhs, const Expr &rhs) const {
    if (lhs.is<Tensor>() && rhs.is<Tensor>()) {
      return (*this)(lhs.as<Tensor>(), rhs.as<Tensor>());
    }

    Comparator cmp;

    return cmp(lhs, rhs);
  }

  bool operator()(const ExprPtr &lhs, const ExprPtr &rhs) const {
    return (*this)(*lhs, *rhs);
  }

 private:
  std::vector<Index> indices_;
};

/// Compares tensor blocks (slots) and specific indices for equality
using IndexSpecificTensorBlockEqualComparator =
    IndexSpecificTensorBlockComparator<std::equal_to, std::not_equal_to, true>;
/// Compares tensor blocks (slots) and specific indices on a less-than
/// relationship
using IndexSpecificTensorBlockLessThanComparator =
    IndexSpecificTensorBlockComparator<std::less, std::not_equal_to, false>;

}  // namespace sequant

#endif  // SEQUANT_UTILITY_TENSOR_HPP
