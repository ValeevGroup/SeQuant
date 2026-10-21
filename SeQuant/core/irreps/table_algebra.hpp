#ifndef SEQUANT_CORE_IRREPS_TABLE_ALGEBRA_HPP
#define SEQUANT_CORE_IRREPS_TABLE_ALGEBRA_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/irreps/algebra.hpp>

#include <string>
#include <utility>

namespace sequant {

/// IRREP algebra defined by an explicit (symbolic) multiplication table. This
/// can be used to describe arbitrary (also non-abelian) groups, with IRREPs
/// being identified by their name.
class TableAlgebra {
 public:
  using irrep_type = std::string;
  /// Maps pairs of IRREPs to their product. For every pair of IRREPs, only
  /// one of the two orderings has to be present. Products involving the
  /// totally symmetric IRREP are implicit.
  using table_type =
      container::map<std::pair<irrep_type, irrep_type>, IrrepSet<irrep_type>>;
  /// Maps IRREPs to their complex conjugate. IRREPs not listed are taken to
  /// be their own complex conjugate.
  using conjugation_map = container::map<irrep_type, irrep_type>;

  /// @param totally_symmetric Name of the totally symmetric IRREP
  /// @param table The multiplication table
  /// @param conjugates The complex conjugation map
  /// @throws Exception if the table is inconsistent
  TableAlgebra(irrep_type totally_symmetric, table_type table,
               conjugation_map conjugates = {});

  const irrep_type &totally_symmetric() const;

  /// @throws Exception if the table doesn't contain the requested product
  IrrepSet<irrep_type> product(const irrep_type &lhs,
                               const irrep_type &rhs) const;

  irrep_type conjugate(const irrep_type &irrep) const;

 private:
  irrep_type m_totallySymmetric;
  table_type m_table;
  conjugation_map m_conjugates;
};

static_assert(IrrepAlgebra<TableAlgebra>);

}  // namespace sequant

#endif  // SEQUANT_CORE_IRREPS_TABLE_ALGEBRA_HPP
