#include <SeQuant/core/export/utils.hpp>

namespace sequant {

bool TensorBlockCompare::operator()(const Tensor &lhs,
                                    const Tensor &rhs) const {
  if (lhs.label() != rhs.label()) {
    return lhs.label() < rhs.label();
  }
  if (lhs.braket().size() != rhs.braket().size()) {
    return lhs.braket().size() < rhs.braket().size();
  }
  auto lhsBraket = lhs.braket();
  auto rhsBraket = rhs.braket();

  for (std::size_t i = 0; i < lhsBraket.size(); ++i) {
    if (lhsBraket.at(i).space() != rhsBraket.at(i).space()) {
      return lhsBraket.at(i).space() < rhsBraket.at(i).space();
    }
  }

  return false;
}

}  // namespace sequant
