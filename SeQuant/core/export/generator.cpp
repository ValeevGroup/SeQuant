#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>

#include <type_traits>

namespace sequant {

PrunableScalars operator&(PrunableScalars lhs, PrunableScalars rhs) {
  using UT = std::underlying_type_t<PrunableScalars>;

  return static_cast<PrunableScalars>(static_cast<UT>(lhs) &
                                      static_cast<UT>(rhs));
}

PrunableScalars operator|(PrunableScalars lhs, PrunableScalars rhs) {
  using UT = std::underlying_type_t<PrunableScalars>;

  return static_cast<PrunableScalars>(static_cast<UT>(lhs) |
                                      static_cast<UT>(rhs));
}

}  // namespace sequant
