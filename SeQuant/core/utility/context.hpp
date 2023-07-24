#ifndef SEQUANT_CORE_UTILITY_CONTEXT_HPP
#define SEQUANT_CORE_UTILITY_CONTEXT_HPP

#include <optional>

/// \name reusable components for manipulation of global contexts

/// @{

namespace sequant::detail {

template <typename Ctx>
inline Ctx& implicit_context_instance() {
  static Ctx instance_;
  return instance_;
}

template <typename Ctx>
const Ctx& get_implicit_context() {
  return implicit_context_instance<Ctx>();
}

template <typename Ctx>
void set_implicit_context(const Ctx& ctx) {
  implicit_context_instance<Ctx>() = ctx;
}

template <typename Ctx>
void reset_implicit_context() {
  implicit_context_instance<Ctx>() = Ctx{};
}

/// used to auto-reset implicit context after leaving scope
template <typename Ctx>
struct ImplicitContextResetter {
  ImplicitContextResetter() = default;
  ImplicitContextResetter(const Ctx& previous_ctx) noexcept
      : previous_ctx_(previous_ctx) {}
  ~ImplicitContextResetter() noexcept {
    if (previous_ctx_) set_implicit_context<Ctx>(*previous_ctx_);
  }

  // ContextResetter is move-only
  ImplicitContextResetter(const ImplicitContextResetter&) = delete;
  ImplicitContextResetter& operator=(const ImplicitContextResetter&) = delete;

 private:
  std::optional<Ctx> previous_ctx_;
};

template <typename Ctx>
ImplicitContextResetter<Ctx> set_scoped_implicit_context(const Ctx& ctx) {
  if (get_implicit_context<Ctx>() != ctx) {
    auto previous_ctx = get_implicit_context<Ctx>();
    set_implicit_context(ctx);
    return previous_ctx;
  } else
    return {};
}

}  // namespace sequant::detail

/// @}

#endif  // SEQUANT_CORE_UTILITY_CONTEXT_HPP
