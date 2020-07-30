#ifndef SEQUANT_FACTORIZER_HPP
#define SEQUANT_FACTORIZER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/tensor_network.hpp>

// IndexSpace type based hashing of tensors for connectivity
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
  explicit AdjacencyMatrix(
      const ExprPtr& expr,
      const TensorNetwork::named_indices_t& external_indices = {});
  /// @param tensors is a vector of ExprPtr to Tensors'
  explicit AdjacencyMatrix(
      const container::svector<ExprPtr>& tensors,
      const TensorNetwork::named_indices_t& external_indices = {});

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

/// Check if a Tensor expression exists in a potentially nested Product of
/// Tensors'.
bool tensor_exists(const ExprPtr& expr, const ExprPtr& tnsr);

/// Get positions of common type of tensors in a pair of Exprs'.
std::tuple<container::set<AdjacencyMatrix::pos_type>,
           container::set<AdjacencyMatrix::pos_type>>
common_tensors(const ExprPtr& expr1, const ExprPtr& expr2);

/// Get the target indices of an expression
TensorNetwork::named_indices_t target_indices(const ExprPtr& expr);

/// Get a map of common sub-expr pair from a pair of AdjacencyMatrix objects.
///
/// @return A tuple-tuple map. eg. {(0, 2): (1, 5), ... }
///         implies the node 0 and 2 from the AdjacencyMatrix mat1 are
///         equivalently connected as the nodes 1 and 5 from mat2 and so on.
container::map<std::tuple<AdjacencyMatrix::pos_type, AdjacencyMatrix::pos_type>,
               std::tuple<AdjacencyMatrix::pos_type, AdjacencyMatrix::pos_type>>
common_pairs(const AdjacencyMatrix& mat1, const AdjacencyMatrix& mat2);

/// Get a tuple of two vectors, where elements of the vectors are the set of
/// pos_types'.
std::tuple<container::svector<container::set<AdjacencyMatrix::pos_type>>,
           container::svector<container::set<AdjacencyMatrix::pos_type>>>
common_nets(const container::map<
            std::tuple<AdjacencyMatrix::pos_type, AdjacencyMatrix::pos_type>,
            std::tuple<AdjacencyMatrix::pos_type, AdjacencyMatrix::pos_type>>&
                pairs);

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
std::tuple<ExprPtr, ExprPtr> factorize_pair(const ExprPtr& expr1,
                                            const ExprPtr& expr2);

///
/// Fuse two factorized Exprs'.
///
/// @param expr1 factorized term to be fused
/// @param expr2 factorized term to be fused
///              expr1 and expr2 are of the form (AB..)(CD..)
///
/// @param symOp ExprPtr to 'S' or 'A' type tensor
/// @return ExprPtr to fused form or nullptr if no fusion is possible.
///
ExprPtr fuse_pair(const ExprPtr& expr1, const ExprPtr& expr2,
                  const ExprPtr& symop, evaluate::Operation op);

/// @param symOp ExprPtr to 'S' or 'A' type tensor
/// @param arithOp evaluate::Operation type
ExprPtr make_intermediate(const ExprPtr& expr1, const ExprPtr& expr2,
                          const ExprPtr& symOp, evaluate::Operation arithOp);

}  // namespace sequant::factorize

#endif /* #ifndef SEQUANT_FACTORIZER_HPP */
