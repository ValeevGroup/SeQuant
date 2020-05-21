#ifndef SEQUANT_FACTORIZER_HPP
#define SEQUANT_FACTORIZER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>

// IndexSpace type based hashing of tensors for ColorMatrix
#include "SeQuant/domain/evaluate/eval_fwd.hpp"

#include <tuple>

///
/// Find common sub-networks between a pair of tensor networks.
/// @author Bimal Gaudel
/// @version May 2020
///

namespace sequant::factorize {
/// Adjacency matrix of an undirected tensor network.
/// @author Bimal Gaudel
/// @version May 16, 2020
class AdjacencyMatrix {
 public:
  using pos_type = std::size_t;
  using color_mat_type =
      container::svector<container::svector<evaluate::HashType>>;

 private:
  /// Color data matrix of an undirected tensor network.
  color_mat_type colorMatrix_;

  /// Constant reference to the color data matrix.
  const color_mat_type& color_mat() const;

 public:
  /// @param expr is a ExprPtr of Tensors'.
  explicit AdjacencyMatrix(const ExprPtr& expr);

  /// @param tensors is a vector of ExprPtr to Tensors'
  AdjacencyMatrix(const container::svector<ExprPtr>& tensors);

  /// Getter for the number of vertices.
  size_t num_verts() const;

  /// Check if pos1 and pos2 are adjacent to each other.
  bool are_connected(pos_type pos1, pos_type pos2) const;

  /// Get the color between the vertices at positions pos1 and pos2.
  color_mat_type::value_type::value_type color(pos_type pos1,
                                               pos_type pos2) const;

  /// Check if two tensors t1 and t2 are connected.
  /// @return True if at least one of the braket label is common in both
  /// tensors.
  static bool are_connected(const ExprPtr& t1, const ExprPtr& t2);
};

/// Check if a Tensor expresion exists in a potentially nested Product of
/// tensors.
bool tensor_exists(const ExprPtr& expr, const ExprPtr& tnsr);

///
/// Factorize common sub-networks from a pair of tensor networks.
///
/// @param exprA Expression to find common subnetwork of. All of the
/// subexpressions of @exprA must be ExprPtr to sequant Tensors'.
///
/// @param exprB Expression to find common subnetwork of. All of the
/// subexpressions of @exprB must be ExprPtr to sequant Tensors'.
///
/// @return Tuple of two ExprPtr that are factorized forms of exprA and exprB.
///
/// @note Only pair of Tensor Products' are supported now.
std::tuple<ExprPtr, ExprPtr> factorize_pair(const ExprPtr& exprA,
                                            const ExprPtr& exprB);
}  // namespace sequant::factorize

#endif /* #ifndef SEQUANT_FACTORIZER_HPP */
