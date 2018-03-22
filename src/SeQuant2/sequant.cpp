#include "sequant.hpp"

namespace sequant2 {

namespace {

SeQuant2 &default_context_instance() {
  static SeQuant2 instance_;
  return instance_;
}

}  // anonymous namespace

const SeQuant2 &default_context() {
  return default_context_instance();
}

void set_default_context(const SeQuant2 &ctx) {
  default_context_instance() = ctx;
}

void reset_default_context() {
  default_context_instance() = SeQuant2{};
}

}  // namespace sequant2
