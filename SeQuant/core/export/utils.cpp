#include <SeQuant/core/export/utils.hpp>

namespace sequant {

bool TensorBlockCompare::operator()(const Tensor &lhs,
                                    const Tensor &rhs) const {
  if (lhs.label() != rhs.label()) {
    return lhs.label() < rhs.label();
  }

  const auto &lhs_indices = lhs.const_indices();
  const auto &rhs_indices = rhs.const_indices();

  if (lhs_indices.size() != rhs_indices.size()) {
    return lhs_indices.size() < rhs_indices.size();
  }

  for (std::size_t i = 0; i < lhs_indices.size(); ++i) {
    if (lhs_indices.at(i).space() != rhs_indices.at(i).space()) {
      return lhs_indices.at(i).space() < rhs_indices.at(i).space();
    }
  }

  return false;
}

}  // namespace sequant
