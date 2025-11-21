//
// Created by Eduard Valeyev on 2/15/23.
//

#ifndef SEQUANT_CORE_SCOPE_HPP
#define SEQUANT_CORE_SCOPE_HPP

#include <utility>

#ifdef SEQUANT_CORE_UTILITY_SCOPE_HPP_USE_STD_EXPERIMENTAL_SCOPE
#error \
    "SEQUANT_CORE_UTILITY_SCOPE_HPP_USE_STD_EXPERIMENTAL_SCOPE already defined, please create an issue at github.com/ValeevGroup/SeQuant"
#endif

// availability of std::experimental::scope_exit requires experimental/scope to
// be present ...
#if __has_include(<experimental/scope>)

// ... and , if using libstdc++ >= 13, C++20 or later
#if defined(__GLIBCXX__) && __GLIBCXX__ >= 13 && ___cplusplus >= 202002L
#define SEQUANT_CORE_UTILITY_SCOPE_HPP_USE_STD_EXPERIMENTAL_SCOPE 1
#endif

#endif

#if SEQUANT_CORE_UTILITY_SCOPE_HPP_USE_STD_EXPERIMENTAL_SCOPE

#include <experimental/scope>

namespace sequant::detail {
using std::experimental::scope_exit;
}

#else

namespace sequant::detail {

template <typename EF>
class scope_exit {
 public:
  template <typename Fn>
  explicit scope_exit(Fn&& fn) noexcept
      : fn_(std::forward<Fn>(fn)), released_(false) {}
  scope_exit(scope_exit&& other) noexcept = default;
  scope_exit(const scope_exit&) = delete;

  scope_exit() noexcept {}

  ~scope_exit() noexcept {
    if (!released_) {
      fn_();
      released_ = true;
    }
  }

  scope_exit& operator=(const scope_exit&) = delete;
  scope_exit& operator=(scope_exit&&) = delete;

  void release() noexcept { released_ = true; }

 private:
  EF fn_;
  bool released_ = false;
};

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_UTILITY_SCOPE_HPP_USE_STD_EXPERIMENTAL_SCOPE
#undef SEQUANT_CORE_UTILITY_SCOPE_HPP_USE_STD_EXPERIMENTAL_SCOPE

namespace sequant::detail {
template <typename EF>
[[nodiscard]] scope_exit<EF> make_scope_exit(EF&& ef) {
  return scope_exit<EF>(std::forward<EF>(ef));
}
}  // namespace sequant::detail

#endif  // SEQUANT_CORE_SCOPE_HPP
