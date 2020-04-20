#ifndef SEQUANT_EVALUATE_EVAL_TREE_NODE_HPP
#define SEQUANT_EVALUATE_EVAL_TREE_NODE_HPP

#include "SeQuant/core/expr_fwd.hpp"

#include "eval_fwd.hpp"

namespace sequant::evaluate {

class EvalTreeNode {
 protected:
  IndexContainer indices_;

  HashType hash_value_{0};

  ScalarType scalar_{1.0};

  virtual HashType hash_node() const = 0;

 public:
  const IndexContainer& indices() const;

  HashType hash_value() const;

  ScalarType scalar() const;

  void scale(ScalarType);

  virtual bool is_leaf() const = 0;
};

class EvalTreeInternalNode : public EvalTreeNode {
 private:
  EvalNodePtr left_{nullptr};

  EvalNodePtr right_{nullptr};

  Operation operation_{Operation::INVALID};

  HashType hash_node() const override;

 public:
  EvalTreeInternalNode(const EvalNodePtr&, const EvalNodePtr&, Operation);

  const EvalNodePtr& left() const;

  const EvalNodePtr& right() const;

  Operation operation() const;

  bool is_leaf() const override;
};

class EvalTreeLeafNode : public EvalTreeNode {
 private:
  ExprPtr expr_{nullptr};

  bool swapped_bra_ket_{false};

  HashType hash_node() const override;

  static bool need_bra_ket_swap(const ExprPtr&);

 public:
  EvalTreeLeafNode(const ExprPtr& tnsr_expr, bool canonize_braket);

  const ExprPtr& expr() const;

  void uncanonize_braket();

  bool is_leaf() const override;
};

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVALUATE_EVAL_TREE_NODE_HPP
