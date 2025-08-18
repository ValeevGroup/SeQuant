#ifndef SEQUANT_UTILITY_TENSOR_HPP
#define SEQUANT_UTILITY_TENSOR_HPP

#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <cassert>
#include <functional>

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
      assert(rhs_it != rhs_indices.end());

      const IndexSpace &left = lhs_it->space();
      const IndexSpace &right = rhs_it->space();

      if (selector(left, right)) {
        return cmp(left, right);
      }
    }

    return fallback;
  }
};

/// Compares tensor blocks (slots) for equality
using TensorBlockEqualComparator =
    TensorBlockComparator<std::equal_to, std::not_equal_to, true>;
/// Compares tensor blocks (slots) on a less-than relationship
using TensorBlockLessThanComparator =
    TensorBlockComparator<std::less, std::not_equal_to, false>;

}  // namespace sequant

#endif  // SEQUANT_UTILITY_TENSOR_HPP
