#include <SeQuant/core/runtime.hpp>

#include <iostream>
#include <locale>

namespace sequant {

void set_locale() {
  // set global C++ locale (and C locale)
  std::locale target_locale{};
  try {  // use en_US.UTF-8, if supported
    target_locale = std::locale{"en_US.UTF-8"};
  } catch (std::exception &) {  // use default if en_US.UTF-8 not available
  }
  std::locale::global(target_locale);
  // set C++ streams locale to target
  std::ios_base::sync_with_stdio(false);
  std::cout.imbue(target_locale);
  std::cerr.imbue(target_locale);
  std::clog.imbue(target_locale);
  std::wcout.imbue(target_locale);
  std::wcerr.imbue(target_locale);
  std::wclog.imbue(target_locale);
  std::ios_base::sync_with_stdio(true);
}

}  // namespace sequant
