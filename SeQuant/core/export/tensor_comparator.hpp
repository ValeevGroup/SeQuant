#ifndef SEQUANT_CORE_EXPORT_TENSOR_COMPARATOR_HPP
#define SEQUANT_CORE_EXPORT_TENSOR_COMPARATOR_HPP

#include <SeQuant/core/index.hpp>

#include <vector>

namespace sequant {

class Tensor;

namespace detail {

class ExportTensorComparator {
 public:
  ExportTensorComparator(std::vector<Index> batchedIndices = {});
  bool operator()(const Tensor &lhs, const Tensor &rhs) const;

  const std::vector<Index> &batch_indices() const;
  void set_batch_indices(std::vector<Index> indices);

 protected:
  std::vector<Index> m_batched;
};

}  // namespace detail
}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_TENSOR_COMPARATOR_HPP
