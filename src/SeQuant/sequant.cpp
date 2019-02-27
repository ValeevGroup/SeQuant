#include "sequant.hpp"

namespace sequant {

namespace detail {

SeQuant &default_context_instance() {
  static SeQuant instance_;
  return instance_;
}

}  // anonymous namespace

const SeQuant &get_default_context() {
  return detail::default_context_instance();
}

void set_default_context(const SeQuant &ctx) {
  detail::default_context_instance() = ctx;
}

void reset_default_context() {
  detail::default_context_instance() = SeQuant{};
}

}  // namespace sequant
