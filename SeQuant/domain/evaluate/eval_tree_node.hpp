#ifndef SEQUANT_EVALUATE_EVAL_TREE_NODE_HPP
#define SEQUANT_EVALUATE_EVAL_TREE_NODE_HPP

#include "SeQuant/core/expr_fwd.hpp"

#include "eval_fwd.hpp"

namespace sequant::evaluate {

class EvalTreeNode {
 protected:
  /// Vector of index labels.
  IndexContainer indices_;

  /// Hash value of the node.
  HashType hash_value_{0};

  /// Scalar that should scale this node's result while evaluating.
  ScalarType scalar_{1.0};

  /// Hashing method for the node. All derived classes should implement their
  /// own hashing methods with this name.
  [[nodiscard]] virtual HashType hash_node() const = 0;

 public:
  /// Default destructor made virtual.
  virtual ~EvalTreeNode() = default;

  /// Getter for the index label container.
  /// \return Const reference to the index label container.
  const IndexContainer& indices() const;

  /// Get the hash value of the node.
  [[nodiscard]] HashType hash_value() const;

  /// Update the hash value of the node.
  void update_hash();

  /// Get the scalar of the node.
  ScalarType scalar() const;

  /// Set the scalar of the node.
  void scale(ScalarType);

  /// Get the LaTex code for the node.
  virtual std::wstring to_latex() const = 0;

  /// Check if the node is a leaf node.
  virtual bool is_leaf() const = 0;
};

class EvalTreeInternalNode : public EvalTreeNode {
 private:
  /// Pointer to the left node.
  EvalNodePtr left_{nullptr};

  /// Pointer to the right node.
  EvalNodePtr right_{nullptr};

  /// Operation type of the internal node.
  Operation operation_{Operation::INVALID};

  /// Hashing method for the internal node.
  HashType hash_node() const override;

 public:
  /// Default destructor made virtual.
  virtual ~EvalTreeInternalNode() = default;

  /// Construct internal node.
  EvalTreeInternalNode(const EvalNodePtr&, const EvalNodePtr&, Operation);

  /// Get left node pointer.
  const EvalNodePtr& left() const;

  /// Get right node pointer.
  const EvalNodePtr& right() const;

  /// Get operation type.
  Operation operation() const;

  /// Get the LaTeX code for the node.
  std::wstring to_latex() const override;

  /// Returns false as internal node is not leaf.
  bool is_leaf() const override;
};

class EvalTreeLeafNode : public EvalTreeNode {
 private:
  /// Pointer to the SeQuant Expr of Tensor type.
  ExprPtr expr_{nullptr};
  //
  /// Hashing method for the leaf node.
  HashType hash_node() const override;

  /// Set true if the bra-ket index labels were swapped during object
  /// construction.
  bool swapped_labels_{false};

  /// Set the bra-ket labels of the node based on the field swapped_labels_
  void set_labels();

  /// Determine if a SeQuant Expr of Tensor type requires swapping index labels
  /// swapped during construction of the object.
  /// \return True as soon as first occurence where the ket Index's IndexSpace
  /// attribute is lexicographically smaller than that of the bra Index's at
  /// corresponding position.
  static bool canonize_swap(const ExprPtr&);

 public:
  /// Default destructor made virtual.
  virtual ~EvalTreeLeafNode() = default;

  /// Construct leaf node.
  /// \param tnsr_expr       shared pointer to the SeQuant Expr of Tensor type.
  //
  /// \param canonize_braket If true, as soon as first occurence where the ket
  ///                        Index's IndexSpace attribute is lexicographically
  ///                        smaller than that of the bra Index's at
  ///                        corresponding position in the bra() and ket() of
  ///                        tnsr_expr is encountered, the index labels will be
  ///                        filled from ket() followed by bra(), in reverse
  ///                        otherwise.

  EvalTreeLeafNode(const ExprPtr& tnsr_expr, bool canonize_braket);

  /// Getter of the pointer to the SeQuant Expr associated to this object.
  const ExprPtr& expr() const;

 /// Get the LaTeX code for the node.
 std::wstring to_latex() const override;

  /// Returns true as this object is a leaf node.
  bool is_leaf() const override;

  /// Swap tensor bra-ket index labels.
  /// eg. if was  T(b_1, b_2, k_1, k_2)
  ///     will be T(k_1, k_2, b_1, b_2)
  /// @note don't forget to update_hash()
  void swap_labels();

  /// Return true if the bra-ket labels have been swapped from how they appear
  /// in the tensor Expr.
  bool swapped_labels() const;

  //
};

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVALUATE_EVAL_TREE_NODE_HPP
