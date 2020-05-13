///
/// A tree data structure of positions of SeQuant Expr objects. Choosing
/// positions to be ordinals makes it easier for eliminating redundant path
/// computations.
///
/// Such a tree for a product of four tensors will have paths that look
/// something like this:
///
/// (0, 1, 2, 3)
/// (0, 1, 3, 2)
/// (0, 2, 1, 3)
/// (0, 2, 3, 1)
/// (0, 3, 1, 2)
/// (0, 3, 2, 1)
/// (1, 2, 0, 3)
/// (1, 2, 3, 0)
/// ((0, 3), (1, 2))
/// (1, 3, 0, 2)
/// (1, 3, 2, 0)
/// ((0, 2), (1, 3))
/// (2, 3, 0, 1)
/// (2, 3, 1, 0)
/// ((0, 1), (2, 3))
///
/// Created by Bimal Gaudel on December 2019
///

#ifndef SEQUANT_EVALUATE_PATH_TREE_HPP
#define SEQUANT_EVALUATE_PATH_TREE_HPP

#include <SeQuant/core/container.hpp>

#include <cstddef>
#include <memory>
#include <string>

namespace sequant::evaluate {

class PathTree;

using PathTreePtr = std::shared_ptr<PathTree>;

///
/// PathTree is a tree data structure where each node has a size_t type data
/// also called a label herein. The label represents the position of a
/// tensor expression in a tensor network.
///
/// Eg.  (T_i * T_j * T_k * T_l * T_m) is a tensor network of tensor
/// contractions, then the positions of the tensors are 0..4 for T_i..T_m.
/// So the corresponding path_tree will be (0, 1, 2, 3).
///
/// ((T_i * T_j * T_k) * T_l * T_m) is yet another tensor network of tensor
/// contractions where
///      Index
///         0 -> child path_tree (0 -> T_i
///                               1 -> T_j
///                               2 -> T_k)
///         1 -> T_l
///         2 -> T_m
///
/// @author Bimal  Gaudel
/// @version December 2019
///
class PathTree {
 public:
  /// Default constructor.
  PathTree() = default;

  /// Construct PathTree from a size_t label.
  ///
  /// @param label The label stored at the root node.
  ///
  explicit PathTree(size_t label);

  /// Destructor is the default destructor.
  virtual ~PathTree() = default;

  /// Getter for the label_ field.
  size_t get_label() const;

  /// Constant reference to the children_ container.
  const sequant::container::svector<PathTreePtr>& get_children() const;

  /// Mutable reference to the children_ container.
  sequant::container::svector<PathTreePtr>& get_children();

  /// Add a children node to this root.
  void add_child(PathTreePtr&);

  /// Remove the last appended child from the children_ container.
  void pop_last_child();

  /// Check if it is the leaf of the tree.
  bool is_leaf() const;

  ///
  /// Get a printable wstring for the tree.
  ///
  /// @return std::wstring for printing.
  ///
  std::wstring print_tree() const;

 private:
  /// The data stored by this node.
  size_t label_{};

  /// The children of this root.
  sequant::container::svector<PathTreePtr> children_;
};

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVALUATE_PATH_TREE_HPP
