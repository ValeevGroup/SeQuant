//
// Created by Eduard Valeyev on 10/14/22.
//

#include <SeQuant/version.hpp>

namespace sequant {

const char* revision() noexcept {
  static const char revision[] = SEQUANT_REVISION;
  return revision;
}

}  // namespace sequant
