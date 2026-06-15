#include <SeQuant/core/export/tensor_comparator.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <algorithm>
#include <vector>

namespace sequant::detail {

ExportTensorComparator::ExportTensorComparator(
    std::vector<Index> batchedIndices)
    : m_batched(std::move(batchedIndices)) {}

const std::vector<Index> &ExportTensorComparator::batch_indices() const {
  return m_batched;
}
void ExportTensorComparator::set_batch_indices(std::vector<Index> indices) {
  m_batched = std::move(indices);
}

bool ExportTensorComparator::operator()(const Tensor &lhs,
                                        const Tensor &rhs) const {
  TensorBlockLessThanComparator cmp;
  bool blocksLess = cmp(lhs, rhs);
  if (blocksLess || cmp(rhs, lhs)) {
    return blocksLess;
  }

  // Blocks are identical -> check for batched indices
  auto &&lhs_indices = lhs.indices();
  auto &&rhs_indices = rhs.indices();

  SEQUANT_ASSERT(lhs.num_indices() == rhs.num_indices());
  auto lit = lhs_indices.begin();
  auto rit = rhs_indices.begin();
  while (lit != lhs_indices.end()) {
    if (*lit != *rit) {
      const std::size_t lhs_pos = std::ranges::distance(
          m_batched.begin(), std::ranges::find(m_batched, *lit));
      const std::size_t rhs_pos = std::ranges::distance(
          m_batched.begin(), std::ranges::find(m_batched, *rit));

      const bool lhs_is_batched = lhs_pos < m_batched.size();
      const bool rhs_is_batched = rhs_pos < m_batched.size();

      if (lhs_is_batched != rhs_is_batched) {
        return lhs_is_batched;
      } else if (lhs_is_batched && rhs_is_batched) {
        return lhs_pos < rhs_pos;
      }
    }

    ++lit;
    ++rit;
  }

  // Tensors are equal
  return false;
}

}  // namespace sequant::detail
