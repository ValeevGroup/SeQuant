//
// Created by Eduard Valeyev on 7/31/23
//

#include <SeQuant/core/utility/macros.hpp>

#include <cstdlib>
#include <iostream>
#include <sstream>

namespace sequant {

[[noreturn]] void throw_failure(const std::string &errmsg,
                                const std::source_location location) {
  std::ostringstream oss;
  oss << errmsg << " at " << location.file_name() << ":" << location.line()
      << " in function '" << location.function_name() << "'";
  throw sequant::Exception(oss.str());
}

void assert_failed(const std::string &errmsg,
                   const std::source_location location) {
  if constexpr (assert_behavior() == AssertBehavior::Throw) {
    throw_failure(errmsg, location);
  } else if constexpr (assert_behavior() == AssertBehavior::Abort) {
    std::cerr << errmsg << " at " << location.file_name() << ":"
              << location.line() << " in function '" << location.function_name()
              << "'" << std::endl;
    std::abort();
  } else {
    std::abort();
  }
}

[[noreturn]] void enforce_failed(const std::string &errmsg,
                                 const std::source_location location) {
#ifdef SEQUANT_ASSERT_ENABLED
  assert_failed(errmsg, location);
#else
  throw_failure(errmsg, location);
#endif
}

[[noreturn]] void abort_msg(const std::string &errmsg,
                            const std::source_location location) {
  std::cerr << errmsg << " at " << location.file_name() << ":"
            << location.line() << " in function '" << location.function_name()
            << "'";
  std::abort();
}

}  // namespace sequant
