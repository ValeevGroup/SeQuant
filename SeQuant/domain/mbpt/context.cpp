#include <SeQuant/domain/mbpt/context.hpp>

namespace sequant::mbpt {

Context::Context(CSV csv) noexcept : csv_(csv) {}

bool operator==(Context const& left, Context const& right) {
  return left.csv() == right.csv();
}

bool operator!=(Context const& left, Context const& right) {
  return !(left == right);
}

const Context& get_default_mbpt_context() {
  return sequant::detail::get_implicit_context<Context>();
}

void set_default_mbpt_context(const Context& ctx) {
  sequant::detail::set_implicit_context(ctx);
}

void reset_default_mbpt_context() {
  sequant::detail::reset_implicit_context<Context>();
}

[[nodiscard]] sequant::detail::ImplicitContextResetter<Context>
set_scoped_default_mbpt_context(const Context& f) {
  return sequant::detail::set_scoped_implicit_context(f);
}

}  // namespace sequant::mbpt
