#ifndef SEQUANT_CORE_UTILITY_CONTEXT_HPP
#define SEQUANT_CORE_UTILITY_CONTEXT_HPP

#include <SeQuant/api.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstring>
#include <iterator>
#include <memory>
#include <typeinfo>
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
/// \warning implicit_context_instance() holds a function-local static; since
/// SeQuant libraries are built with hidden symbol visibility, each shared
/// library has its own copy of it. Hence getting/setting the process-wide
/// context of a given `Ctx` must only be instantiated in (the non-inline
/// functions of) the one library that owns `Ctx`. The thread's overlays
/// (see implicit_context_overlays()) and their keys (see
/// implicit_context_key()) are shared by all libraries, i.e. querying the
/// overlays of a `Ctx` is fine anywhere.

/// @{

namespace sequant::detail {

template <typename Ctx>
inline Ctx& implicit_context_instance() {
  static Ctx instance_;
  return instance_;
}

/// a scoped overlay of the implicit context whose type is identified by `key`
/// (see implicit_context_key())
struct ImplicitContextOverlay {
  const char* key;
  std::shared_ptr<const void> ctx;
};

/// the overlays of this thread, innermost last
using ImplicitContextOverlays = std::vector<ImplicitContextOverlay>;

/// @return the overlays of this thread (shared by all SeQuant libraries)
SEQUANT_SYMB_INTERNAL_API ImplicitContextOverlays& implicit_context_overlays();

/// @return the key that identifies the overlays of `Ctx`
/// @note the key is the (implementation-defined, but unique) name of `Ctx`
/// rather than an address: those would differ between shared libraries
/// @sa implicit_context_key_equal()
template <typename Ctx>
const char* implicit_context_key() {
  return typeid(Ctx).name();
}

/// @return whether @p key1 and @p key2 identify the same `Ctx`
inline bool implicit_context_key_equal(const char* key1, const char* key2) {
  return key1 == key2 || std::strcmp(key1, key2) == 0;
}

/// @return the innermost overlay of this thread for `Ctx`, or nullptr
/// @note the overlay lives at least as long as it is on this thread's stack
template <typename Ctx>
const Ctx* implicit_context_overlay() {
  const auto& overlays = implicit_context_overlays();
  for (auto it = overlays.rbegin(); it != overlays.rend(); ++it) {
    if (implicit_context_key_equal(it->key, implicit_context_key<Ctx>()))
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
