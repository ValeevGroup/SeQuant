//
// Created by Eduard Valeyev on 7/30/25.
//

#ifndef SEQUANT_CORE_UTILITY_DEBUG_HPP
#define SEQUANT_CORE_UTILITY_DEBUG_HPP

#include <cstdio>
#include <sstream>

namespace sequant {

/// uses std::wostringstream + std::wprintf to print to stdout even if wcout is
/// captured (e.g., by Catch2)
template <typename... Args>
void wprintf(Args&&... args) {
  std::wostringstream oss;
  (oss << ... << std::forward<Args>(args)) << std::endl;
  std::wprintf(oss.str().c_str());
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_DEBUG_HPP
