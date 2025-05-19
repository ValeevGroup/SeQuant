#include <SeQuant/core/export/compute_selection.hpp>

#include <type_traits>

namespace sequant {

ComputeSelection operator|(ComputeSelection lhs, ComputeSelection rhs) {
  using underlying_type = std::underlying_type_t<ComputeSelection>;

  return static_cast<ComputeSelection>(static_cast<underlying_type>(lhs) |
                                       static_cast<underlying_type>(rhs));
}

ComputeSelection operator&(ComputeSelection lhs, ComputeSelection rhs) {
  using underlying_type = std::underlying_type_t<ComputeSelection>;

  return static_cast<ComputeSelection>(static_cast<underlying_type>(lhs) &
                                       static_cast<underlying_type>(rhs));
}

ComputeSelection &operator|=(ComputeSelection &lhs, ComputeSelection rhs) {
  lhs = lhs | rhs;
  return lhs;
}

ComputeSelection &operator&=(ComputeSelection &lhs, ComputeSelection rhs) {
  lhs = lhs & rhs;
  return lhs;
}

ComputeSelection operator~(ComputeSelection selection) {
  using underlying_type = std::underlying_type_t<ComputeSelection>;

  return static_cast<ComputeSelection>(
      ~static_cast<underlying_type>(selection));
}

}  // namespace sequant
