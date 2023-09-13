#include "imed.hpp"

namespace sequant {

Imed::Imed(size_t k) : key{k}, flops{}, memory{}, pos{}, expr{} {}

bool Imed::CompKey::operator()(Imed const& left,
                               Imed const& right) const noexcept {
  return left.key < right.key;
}

bool Imed::CompKey::operator()(size_t k, const Imed& right) const noexcept {
  return k < right.key;
}

bool Imed::CompKey::operator()(const Imed& left, size_t k) const noexcept {
  return left.key < k;
}

bool operator<(Imed const& left, Imed const& right) noexcept {
  return left.key < right.key;
}
}  // namespace sequant