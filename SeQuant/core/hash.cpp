#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>

namespace sequant {

std::size_t hash_value(const ExprPtr &expr) { return hash_value(*expr); }

}  // namespace sequant
