#include "sequant.hpp"

namespace sequant {

bool operator==(const SeQuant& ctx1, const SeQuant& ctx2) {
  return ctx1.vacuum() == ctx2.vacuum() && ctx1.metric() == ctx2.metric() &&
         ctx1.braket_symmetry() == ctx2.braket_symmetry() &&
         ctx1.spbasis() == ctx2.spbasis();
}

bool operator!=(const SeQuant& ctx1, const SeQuant& ctx2) {
  return !(ctx1 == ctx2);
}

namespace detail {

SeQuant& default_context_instance() {
  static SeQuant instance_;
  return instance_;
}

}  // namespace detail

const SeQuant& get_default_context() {
  return detail::default_context_instance();
}

void set_default_context(const SeQuant& ctx) {
  detail::default_context_instance() = ctx;
}

void reset_default_context() { detail::default_context_instance() = SeQuant{}; }

detail::ContextResetter::ContextResetter() noexcept
    : previous_ctx_(get_default_context()) {}

detail::ContextResetter::~ContextResetter() noexcept {
  set_default_context(previous_ctx_);
}

detail::ContextResetter set_scoped_default_context(const SeQuant& ctx) {
  detail::ContextResetter resetter;
  detail::default_context_instance() = ctx;
  return resetter;
}

}  // namespace sequant
