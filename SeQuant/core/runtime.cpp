#include <SeQuant/core/runtime.hpp>

#include <exception>
#include <iostream>
#include <locale>

namespace sequant {

void set_locale() {
  // set global C++ locale (and C locale)
  std::locale target_locale{};
  try {  // use C.UTF-8, if supported
    target_locale = std::locale{"C.UTF-8"};
  } catch (std::exception &) {  // use default if C.UTF-8 not available
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

namespace {
// see https://en.cppreference.com/w/cpp/locale/numpunct/thousands_sep
class no_thousands_separator : public std::numpunct<char> {
 protected:
  char do_thousands_sep() const override { return '\0'; }
  std::string do_grouping() const override { return ""; }
};

class no_thousands_separator_w : public std::numpunct<wchar_t> {
 protected:
  wchar_t do_thousands_sep() const override { return L'\0'; }
  std::string do_grouping() const override { return ""; }
};
}  // namespace

void disable_thousands_separator() {
  auto current_locale = std::locale();

  std::locale modified_locale(current_locale, new no_thousands_separator);
  modified_locale = std::locale(modified_locale, new no_thousands_separator_w);

  std::locale::global(modified_locale);
}

}  // namespace sequant
