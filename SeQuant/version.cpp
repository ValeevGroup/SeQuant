//
// Created by Eduard Valeyev on 10/14/22.
//

#include <SeQuant/version.hpp>

namespace sequant {

const char* revision() noexcept {
  static const char revision[] = SEQUANT_GIT_REVISION;
  return revision;
}

const char* git_description() noexcept {
  static const char description[] = SEQUANT_GIT_DESCRIPTION;
  return description;
}

}  // namespace sequant
