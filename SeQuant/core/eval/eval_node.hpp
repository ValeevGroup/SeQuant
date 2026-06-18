//
// Created by Bimal Gaudel on 5/24/21.
//

#ifndef SEQUANT_EVAL_EVAL_NODE_HPP
#define SEQUANT_EVAL_EVAL_NODE_HPP

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/count_if.hpp>

namespace sequant {

template <meta::eval_node Node>
ExprPtr linearize_eval_node(Node const& node) {
  if (node.leaf()) return to_expr(node);

  // Adjoint is unary: round-trip back to the marker-bearing tensor stored in
  // the node's own ExprPtr — that's the symbolic form a Sum/Product parent
  // would see — and ignore the Constant(1) sentinel right child.
  if (node->op_type() == EvalOp::Adjoint) return to_expr(node);

  ExprPtr lres = linearize_eval_node(node.left());
  ExprPtr rres = linearize_eval_node(node.right());

  SEQUANT_ASSERT(lres);
  SEQUANT_ASSERT(rres);

  if (node->op_type() == EvalOp::Sum) return ex<Sum>(ExprPtrList{lres, rres});
  SEQUANT_ASSERT(node->op_type() == EvalOp::Product);
  return ex<Product>(
      Product{1, ExprPtrList{lres, rres}, Product::Flatten::Yes});
}

// implementation details of eval-node cost analysis; prefer sequant::detail
// over an unnamed namespace in a header (see CppCoreGuidelines SF.21)
namespace detail {

enum NodePos { Left = 0, Right, This };

/// Tally indices per IndexSpace appearing in the bra+ket+aux of `t`. The aux
/// slots are included so that auxiliary spaces (e.g. density-fitting/THC) show
/// up in the cost. AsyCost is only used for analysis/reporting, so spaces are
/// taken verbatim from the indices, with no registry lookup.
///
/// Cost is expressed in terms of whatever IndexSpace each index carries. A
/// union-space index (e.g. MR hole `I=i∪u`) contributes its union dimension
/// `|I|=|i|+|u|` as a single symbol: it is not decomposed into base-space (`i`,
/// `u`, ...) block terms.
inline AsyCost::ExponentMap space_counts(Tensor const& t) {
  AsyCost::ExponentMap counts;
  for (auto const& idx : t.const_braketaux_indices()) ++counts[idx.space()];
  return counts;
}

class ContractedIndexCount {
 public:
  using Counts = AsyCost::ExponentMap;

  explicit ContractedIndexCount(meta::eval_node auto const& n) {
    auto const L = NodePos::Left;
    auto const R = NodePos::Right;
    auto const T = NodePos::This;

    SEQUANT_ASSERT(n->is_tensor() && n.left()->is_tensor() &&
                   n.right()->is_tensor());

    // For the result and both operands, record the tensor's rank and collect
    // every index into `distinct`. The set dedupes indices shared across
    // tensors (a contracted index appears in both operands; a batched index in
    // all three), so each is counted exactly once below.
    container::set<Index> distinct;
    for (auto p : {L, R, T}) {
      auto const& t = (p == L ? n.left() : p == R ? n.right() : n)->as_tensor();
      auto const counts = space_counts(t);
      ranks_[p] = 0;
      for (auto const& [_, c] : counts) ranks_[p] += c;
      for (auto const& idx : t.const_braketaux_indices()) distinct.insert(idx);
    }

    // An outer product contracts nothing, so the result keeps all operand
    // indices and the ranks add up exactly.
    is_outerprod_ = ranks_[L] + ranks_[R] == ranks_[T];

    // Cost exponent per space = number of distinct indices touching that space.
    for (auto const& idx : distinct) ++unique_[idx.space()];
  }

  [[nodiscard]] bool is_outerprod() const noexcept { return is_outerprod_; }

  /// Per-space count of distinct indices participating in the contraction.
  /// Each external, contracted, or batched index is counted once. For ordinary
  /// contractions (no batched index) this equals (L[s] + R[s] + T[s]) / 2.
  [[nodiscard]] Counts const& unique_counts() const { return unique_; }

 private:
  Counts unique_;
  std::array<size_t, 3> ranks_{0, 0, 0};
  bool is_outerprod_ = false;
};
}  // namespace detail

///
/// \brief This function object takes an evaluation node and returns the
///        symbolic cost of flops required for evaluation, as an AsyCost object.
///        @see AsyCost.
/// \note
///        - The cost of evaluation of leaf nodes is assumed to be zero.
///        - Cost of evaluation of children nodes are not counted.
///
struct Flops {
  [[nodiscard]] AsyCost operator()(meta::eval_node auto const& n) const {
    using detail::space_counts;
    if (n.leaf()) return AsyCost::zero();
    if (n->op_type() == EvalOp::Product  //
        && n.left()->is_tensor()         //
        && n.right()->is_tensor()) {
      if (n->is_tensor()) {
        auto const idx_count = detail::ContractedIndexCount{n};
        auto c = AsyCost{idx_count.unique_counts()};
        return idx_count.is_outerprod() ? c : 2 * c;
      } else {  // full contraction to scalar
        SEQUANT_ASSERT(n->is_scalar());
        SEQUANT_ASSERT(space_counts(n.left()->as_tensor()) ==
                       space_counts(n.right()->as_tensor()));
        return 2 * AsyCost{space_counts(n.left()->as_tensor())};
      }
    } else if (n->is_tensor()) {
      // scalar times a tensor
      // or a tensor plus a tensor
      return AsyCost{space_counts(n->as_tensor())};
    } else /* scalar (+|*) scalar */
      return AsyCost::zero();
  }
};

///
/// \brief This function object takes an evaluation node and returns the
///        memory storage required to evaluate it in AsyCost form.
/// \note
///       - The memory requirement for non-tensor objects (eg. variables and
///         scalar constants) are taken to be zero.
///       - The memory requirement for evaluation of children nodes is not
///         counted.
///
struct Memory {
  [[nodiscard]] AsyCost operator()(meta::eval_node auto const& n) const {
    AsyCost result;
    auto add_cost = [&result](ExprPtr const& expr) {
      result += expr.is<Tensor>()
                    ? AsyCost{detail::space_counts(expr.as<Tensor>())}
                    : AsyCost::zero();
    };

    // Adjoint is unary — the right child is the Constant(1) sentinel (zero
    // cost). But the backend materializes the permuted/conjugated result into
    // a fresh array, so the bare-leaf operand (left) and the adjoint result
    // (node) are both live at peak — two tensors' worth, not one.
    if (n->op_type() == EvalOp::Adjoint) {
      add_cost(n.left()->expr());
      add_cost(n->expr());
      return result;
    }

    add_cost(n.left()->expr());
    add_cost(n.right()->expr());
    add_cost(n->expr());
    return result;
  }
};

///
/// \brief This function object takes an evaluation node and returns the
///        symbolic cost of flops required for evaluation as an AsyCost object.
///        @see AsyCost. If the cost can be reduced due to symmetry, it is done
///        so.
/// \note The cost of evaluation of leaf nodes is assumed to be zero.
///
struct FlopsWithSymm {
  [[nodiscard]] AsyCost operator()(meta::eval_node auto const& n) const {
    auto cost = Flops{}(n);
    // Adjoint is unary (right is Constant(1) sentinel) and does no
    // symmetry-driven reduction — return Flops's permute-style cost as-is.
    if (n.leaf() || n->op_type() == EvalOp::Adjoint ||
        !(n->is_tensor()            //
          && n.left()->is_tensor()  //
          && n.right()->is_tensor()))
      return cost;

    // confirmed:
    // left, right and this node
    // all have tensor expression
    auto const& t = n->as_tensor();
    auto const tsymm = t.symmetry();
    //

    // ------
    // the rules of cost reduction are taken from
    //   doi:10.1016/j.procs.2012.04.044
    // ------
    if (tsymm == Symmetry::Symm || tsymm == Symmetry::Antisymm) {
      auto const op = n->op_type();
      auto const tbrank = t.bra_rank();
      auto const tkrank = t.ket_rank();
      if (op == EvalOp::Sum)
        cost = tsymm == Symmetry::Symm
                   ? cost / (factorial(tbrank) * factorial(tkrank))
                   : cost / factorial(tbrank);
      else if (op == EvalOp::Product) {
        auto const lsymm = n.left()->as_tensor().symmetry();
        auto const rsymm = n.right()->as_tensor().symmetry();
        cost = (lsymm == rsymm && lsymm == Symmetry::Nonsymm)
                   ? cost / factorial(tbrank)
                   : cost / (factorial(tbrank) * factorial(tkrank));
      } else
        SEQUANT_ASSERT(false &&
                       "Unsupported evaluation operation for asymptotic cost "
                       "computation.");
    }
    return cost;
  }
};

///
/// \brief This function object takes an evaluation node and returns the cost of
///        evaluating the node as an AsyCost object. The cost is computed by
///        summing the cost of evaluation of children nodes and the cost of
///        evaluation of the node itself.
/// \param node The evaluation node whose cost to be evaluated
/// \param cost_fn A function object that takes an evaluation node and returns
///                the symbolic cost of flops required for evaluation as an
///                AsyCost object. @see AsyCost.
/// \return The asymptotic cost of evaluating the given node.
///
template <meta::eval_node Node, typename F = Flops>
  requires requires(F const& fn, Node const& n) {
    { fn(n) } -> std::same_as<AsyCost>;
  }
AsyCost asy_cost(Node const& node, F const& cost_fn = {}) {
  return node.leaf() ? cost_fn(node)
                     : asy_cost(node.left(), cost_fn) +
                           asy_cost(node.right(), cost_fn) + cost_fn(node);
}

///
/// \brief This function object takes an evaluation node and returns the minimum
///        storage required for evaluating the node as an AsyCost object. The
///        minimum storage is the largest amount of storage required for
///        evaluating the node and its children.
/// \return The minimum storage required for evaluating the given node.
///
AsyCost min_storage(meta::eval_node auto const& node) {
  using detail::space_counts;
  auto result = AsyCost::zero();
  auto visitor = [&result](meta::eval_node auto const& n) {
    auto cost = AsyCost::zero();
    if (n.leaf() && n->is_tensor())
      cost = AsyCost{space_counts(n->as_tensor())};
    else if (!n.leaf()) {
      cost +=
          (n.left()->is_tensor() ? AsyCost{space_counts(n.left()->as_tensor())}
                                 : AsyCost::zero());
      cost += (n.right()->is_tensor()
                   ? AsyCost{space_counts(n.right()->as_tensor())}
                   : AsyCost::zero());
      cost += (n->is_tensor() ? AsyCost{space_counts(n->as_tensor())}
                              : AsyCost::zero());
    } else {
      // do nothing
    }
    result = std::max(result, cost);
  };
  node.visit(visitor);
  return result;
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_NODE_HPP
