#ifndef SEQUANT_CORE_OPTIMIZE_SUM_HPP
#define SEQUANT_CORE_OPTIMIZE_SUM_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>

namespace sequant::opt {

///
/// \brief Create clusters out of positions of terms in a sum that share common
///        intermediates.
///
/// \param expr A Sum to find clusters in.
/// \return A vector of clusters (vectors of position index of terms
///         in \c expr).
///
container::vector<container::vector<size_t>> clusters(Sum const& expr);

///
/// \brief Reorder summands so that terms having common intermediates appear
///        closer.
///
/// \param sum Expression to reorder.
/// \return Expression with summands re-ordered.
///
Sum reorder(Sum const& sum);

}  // namespace sequant::opt

#endif  // SEQUANT_CORE_OPTIMIZE_SUM_HPP
