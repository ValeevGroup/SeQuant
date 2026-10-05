#ifndef SEQUANT_CORE_UTILITY_CONTEXT_HPP
#define SEQUANT_CORE_UTILITY_CONTEXT_HPP

#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

/// \name reusable components for manipulation of global contexts
///
/// An implicit context of type `Ctx` is a process-wide object
/// (implicit_context_instance()) that a thread can shadow by scoped overlays
/// (set_scoped_implicit_context()). Overlays are thread-local and are
/// propagated into the workers of the parallel primitives in
/// SeQuant/core/runtime.hpp, so a scoped context is seen by the code that
/// runs on behalf of its thread but never by other threads.
/// \note a process-wide change (set_implicit_context(),
/// reset_implicit_context()) is not seen by a thread with overlays until its
/// scopes end

/// @{

namespace sequant::detail {

template <typename Ctx>
inline Ctx& implicit_context_instance() {
  static Ctx instance_;
  return instance_;
}

/// a scoped overlay of the implicit context whose type is identified by `key`
struct ImplicitContextOverlay {
  const void* key;
  std::shared_ptr<const void> ctx;
};

/// the overlays of this thread, innermost last
using ImplicitContextOverlays = std::vector<ImplicitContextOverlay>;

inline ImplicitContextOverlays& implicit_context_overlays() {
  thread_local ImplicitContextOverlays overlays;
  return overlays;
}

template <typename Ctx>
const void* implicit_context_key() {
  static const char key = 0;
  return &key;
}

/// @return the innermost overlay of this thread for `Ctx`, or nullptr
/// @note the overlay lives at least as long as it is on this thread's stack
template <typename Ctx>
const Ctx* implicit_context_overlay() {
  const auto& overlays = implicit_context_overlays();
  for (auto it = overlays.rbegin(); it != overlays.rend(); ++it) {
    if (it->key == implicit_context_key<Ctx>())
      return static_cast<const Ctx*>(it->ctx.get());
  }
  return nullptr;
}

/// @return the innermost overlay of this thread for `Ctx`, if any, else the
/// process-wide context
template <typename Ctx>
const Ctx& get_implicit_context() {
  if (const auto* overlay = implicit_context_overlay<Ctx>()) return *overlay;
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

/// pops the overlay pushed by set_scoped_implicit_context() when leaving scope
template <typename Ctx>
struct ImplicitContextResetter {
  ImplicitContextResetter() = default;
  explicit ImplicitContextResetter(const void* overlay) noexcept
      : overlay_(overlay) {}
  /// @note a scope ended out of order trips an assertion; since this is
  /// noexcept, under SEQUANT_ASSERT_BEHAVIOR=THROW that calls std::terminate
  ~ImplicitContextResetter() noexcept {
    if (overlay_) {
      auto& overlays = implicit_context_overlays();
      const auto it = std::find_if(
          overlays.rbegin(), overlays.rend(),
          [this](const auto& o) { return o.ctx.get() == overlay_; });
      SEQUANT_ASSERT(it != overlays.rend() && it == overlays.rbegin() &&
                     "scoped contexts must end in the reverse order of their "
                     "creation, on the thread that created them");
      if (it != overlays.rend()) overlays.erase(std::next(it).base());
    }
  }

  // neither copyable nor movable: a moved-from resetter would still hold
  // overlay_ and pop the live scope; returning by value relies on guaranteed
  // copy elision
  ImplicitContextResetter(const ImplicitContextResetter&) = delete;
  ImplicitContextResetter(ImplicitContextResetter&&) = delete;
  ImplicitContextResetter& operator=(const ImplicitContextResetter&) = delete;
  ImplicitContextResetter& operator=(ImplicitContextResetter&&) = delete;

 private:
  const void* overlay_ = nullptr;
};

/// makes @p ctx this thread's implicit context of type `Ctx` until the
/// returned object is destroyed
template <typename Ctx>
ImplicitContextResetter<Ctx> set_scoped_implicit_context(Ctx ctx) {
  auto overlay = std::make_shared<const Ctx>(std::move(ctx));
  const void* overlay_ptr = overlay.get();
  implicit_context_overlays().push_back(
      {implicit_context_key<Ctx>(), std::move(overlay)});
  return ImplicitContextResetter<Ctx>(overlay_ptr);
}

/// makes @p overlays (captured on another thread with
/// implicit_context_overlays()) this thread's overlays until destroyed
class ImplicitContextOverlaysScope {
 public:
  explicit ImplicitContextOverlaysScope(const ImplicitContextOverlays& overlays)
      : active_(!overlays.empty() || !implicit_context_overlays().empty()) {
    if (active_)
      previous_ = std::exchange(implicit_context_overlays(), overlays);
  }
  ~ImplicitContextOverlaysScope() {
    if (active_) implicit_context_overlays() = std::move(previous_);
  }

  ImplicitContextOverlaysScope(const ImplicitContextOverlaysScope&) = delete;
  ImplicitContextOverlaysScope& operator=(const ImplicitContextOverlaysScope&) =
      delete;

 private:
  bool active_;
  ImplicitContextOverlays previous_;
};

}  // namespace sequant::detail

/// @}

#endif  // SEQUANT_CORE_UTILITY_CONTEXT_HPP
