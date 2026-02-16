//
// Created by Eduard Valeyev on 7/31/23
//

#include <SeQuant/core/utility/macros.hpp>

#include <cstdlib>
#include <iostream>
#include <sstream>

#define SEQUANT_ASSERT_BEHAVIOR \
  SEQUANT_CONCAT(SEQUANT_ASSERT_, SEQUANT_ASSERT_BEHAVIOR_)

namespace sequant {

bool assert_enabled() {
#ifdef SEQUANT_ASSERT_ENABLED
  return true;
#else
  return false;
#endif
}

int assert_behavior() { return SEQUANT_ASSERT_BEHAVIOR; }

#ifdef SEQUANT_ASSERT_ENABLED
void assert_failed(const std::string &errmsg,
                   const std::source_location location) {
#if SEQUANT_ASSERT_BEHAVIOR == SEQUANT_ASSERT_THROW
  std::ostringstream oss;
  oss << errmsg << " at " << location.file_name() << ":" << location.line()
      << " in function '" << location.function_name() << "'";
  throw sequant::Exception(oss.str());
#elif SEQUANT_ASSERT_BEHAVIOR == SEQUANT_ASSERT_ABORT
  std::cerr << errmsg << " at " << location.file_name() << ":"
            << location.line() << " in function '" << location.function_name()
            << "'" << std::endl;
  std::abort();
#endif
}
#else
void assert_failed(const std::string &, const std::source_location) {}
#endif

[[noreturn]] void abort_msg(const std::string &errmsg,
                            const std::source_location location) {
  std::cerr << errmsg << " at " << location.file_name() << ":"
            << location.line() << " in function '" << location.function_name()
            << "'";
  std::abort();
}

}  // namespace sequant
