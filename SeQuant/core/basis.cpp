#include <SeQuant/core/basis.hpp>

#include <string>

namespace sequant {

std::wstring IndexBasis::instance_suffix() const {
  return basis_instance_ ? L";" + std::to_wstring(*basis_instance_)
                         : std::wstring{};
}

bool includes(const IndexBasis& basis, const IndexBasis& subbasis) {
  return includes(basis.space(), subbasis.space()) &&
         (!basis.has_basis_instance() ||
          basis.basis_instance() == subbasis.basis_instance());
}

bool different_instances(const IndexBasis& b1, const IndexBasis& b2) {
  return b1.has_basis_instance() && b2.has_basis_instance() &&
         b1.basis_instance() != b2.basis_instance();
}

}  // namespace sequant
