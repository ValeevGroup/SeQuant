#ifndef SEQUANT_CORE_OPTIMIZE_SUM_HPP
#define SEQUANT_CORE_OPTIMIZE_SUM_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
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

/// Same as \ref clusters(Sum const&) but uses precomputed eval nodes for the
/// summands, avoiding an internal \c binarize pass. \p nodes[i] must be the
/// binarized form of \c expr.at(i).
container::vector<container::vector<size_t>> clusters(
    Sum const& expr, container::vector<FullBinaryNode<EvalExpr>> const& nodes);

///
/// \brief Reorder summands so that terms having common intermediates appear
///        closer.
///
/// \param sum Expression to reorder.
/// \return Expression with summands re-ordered.
///
Sum reorder(Sum const& sum);

/// Same as \ref reorder(Sum const&) but uses precomputed eval nodes for the
/// summands. \p nodes[i] must be the binarized form of \c sum.at(i).
Sum reorder(Sum const& sum,
            container::vector<FullBinaryNode<EvalExpr>> const& nodes);

}  // namespace sequant::opt

#endif  // SEQUANT_CORE_OPTIMIZE_SUM_HPP
