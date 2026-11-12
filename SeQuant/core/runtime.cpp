#include <SeQuant/core/runtime.hpp>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <locale>
#include <thread>

namespace sequant {

namespace detail {
int& nthreads_accessor() {
  auto init_nthreads = []() {
    const auto nthreads_cstr = std::getenv("SEQUANT_NUM_THREADS");
    int nthreads = nthreads_cstr ? string_to<unsigned int>(nthreads_cstr)
                                 : (std::thread::hardware_concurrency() > 0
                                        ? std::thread::hardware_concurrency()
                                        : 1);
    return nthreads;
  };
  static int nthreads = init_nthreads();
  return nthreads;
}
}  // namespace detail

void set_locale() {
  // set global C++ locale (and C locale)
  std::locale target_locale{};
  try {  // use C.UTF-8, if supported
    target_locale = std::locale{"C.UTF-8"};
  } catch (std::exception&) {  // use default if C.UTF-8 not available
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
