#include <SeQuant/core/irreps/table_algebra.hpp>
#include <SeQuant/core/utility/exception.hpp>

namespace sequant {

TableAlgebra::TableAlgebra(irrep_type totally_symmetric, table_type table,
                           conjugation_map conjugates)
    : m_totallySymmetric(std::move(totally_symmetric)),
      m_table(std::move(table)),
      m_conjugates(std::move(conjugates)) {
  for (const auto &[factors, result] : m_table) {
    const auto &[lhs, rhs] = factors;

    if (result.is_any()) {
      throw Exception("Product of IRREPs '" + lhs + "' and '" + rhs +
                      "' must not be 'any'");
    }

    if (lhs == m_totallySymmetric || rhs == m_totallySymmetric) {
      const irrep_type &other = lhs == m_totallySymmetric ? rhs : lhs;
      if (result != IrrepSet<irrep_type>{other}) {
        throw Exception("Product of IRREP '" + other +
                        "' with the totally symmetric IRREP '" +
                        m_totallySymmetric + "' must yield '" + other +
                        "' but table says " + to_string(result));
      }
    }

    auto it = m_table.find(std::make_pair(rhs, lhs));
    if (it != m_table.end() && it->second != result) {
      throw Exception("Multiplication table is not commutative for IRREPs '" +
                      lhs + "' and '" + rhs + "': " + to_string(result) +
                      " vs. " + to_string(it->second));
    }
  }
}

const TableAlgebra::irrep_type &TableAlgebra::totally_symmetric() const {
  return m_totallySymmetric;
}

IrrepSet<TableAlgebra::irrep_type> TableAlgebra::product(
    const irrep_type &lhs, const irrep_type &rhs) const {
  if (lhs == m_totallySymmetric) {
    return {rhs};
  }
  if (rhs == m_totallySymmetric) {
    return {lhs};
  }

  auto it = m_table.find(std::make_pair(lhs, rhs));
  if (it == m_table.end()) {
    it = m_table.find(std::make_pair(rhs, lhs));
  }
  if (it == m_table.end()) {
    throw Exception("Multiplication table doesn't contain the product of '" +
                    lhs + "' and '" + rhs + "'");
  }

  return it->second;
}

TableAlgebra::irrep_type TableAlgebra::conjugate(
    const irrep_type &irrep) const {
  auto it = m_conjugates.find(irrep);

  return it == m_conjugates.end() ? irrep : it->second;
}

}  // namespace sequant
