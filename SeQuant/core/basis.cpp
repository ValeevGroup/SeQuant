#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <string>

namespace sequant {

IndexBasis::IndexBasis(IndexSpace space, optional_instance basis_instance,
                       std::wstring name, std::optional<std::size_t> extent,
                       IndexSpaceMetric metric, std::optional<Field> field)
    : space_(std::move(space)),
      basis_instance_(basis_instance),
      name_(std::move(name)),
      extent_(extent),
      metric_(metric),
      field_(field) {
  SEQUANT_ASSERT(name_.empty() || basis_instance_);
  SEQUANT_ASSERT(!name_.empty() ||
                 (!extent_ && metric_ == IndexSpaceMetric::Unit && !field_));
}

std::wstring IndexBasis::instance_suffix() const {
  return basis_instance_ ? L";" + std::to_wstring(*basis_instance_)
                         : std::wstring{};
}

bool same_instance(const IndexBasis& b1, const IndexBasis& b2) noexcept {
  return b1.space() == b2.space() && b1.basis_instance() == b2.basis_instance();
}

bool includes(const IndexBasis& basis, const IndexBasis& subbasis) {
  return includes(basis.space(), subbasis.space()) &&
         (!basis.has_basis_instance() ||
          (basis.basis_instance() == subbasis.basis_instance() &&
           basis.name() == subbasis.name()));
}

bool different_instances(const IndexBasis& b1, const IndexBasis& b2) {
  return b1.has_basis_instance() && b2.has_basis_instance() &&
         (b1.basis_instance() != b2.basis_instance() || b1.name() != b2.name());
}

}  // namespace sequant
