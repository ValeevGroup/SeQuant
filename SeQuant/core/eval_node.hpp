//
// Created by Bimal Gaudel on 5/24/21.
//

#ifndef SEQUANT_EVAL_NODE_HPP
#define SEQUANT_EVAL_NODE_HPP

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/math.hpp>

namespace sequant {

template <meta::eval_node Node>
ExprPtr linearize_eval_node(Node const& node) {
  if (node.leaf()) return to_expr(node);

  ExprPtr lres = linearize_eval_node(node.left());
  ExprPtr rres = linearize_eval_node(node.right());

  assert(lres);
  assert(rres);

  if (node->op_type() == EvalOp::Sum) return ex<Sum>(ExprPtrList{lres, rres});
  assert(node->op_type() == EvalOp::Product);
  return ex<Product>(
      Product{1, ExprPtrList{lres, rres}, Product::Flatten::Yes});
}

namespace {

enum NodePos { Left = 0, Right, This };

[[maybe_unused]] std::pair<size_t, size_t> occ_virt(Tensor const& t) {
  auto bk_rank = t.bra_net_rank() + t.ket_net_rank();
  auto nocc = ranges::count_if(t.const_braket_indices(), [](Index const& idx) {
    return idx.space() ==
           get_default_context().index_space_registry()->hole_space(
               idx.space().qns());
  });
  auto nvirt = bk_rank - nocc;
  return {nocc, nvirt};
}

class ContractedIndexCount {
 public:
  explicit ContractedIndexCount(meta::eval_node auto const& n) {
    auto const L = NodePos::Left;
    auto const R = NodePos::Right;
    auto const T = NodePos::This;

    assert(n->is_tensor() && n.left()->is_tensor() && n.right()->is_tensor());

    for (auto p : {L, R, T}) {
      auto const& t = (p == L ? n.left() : p == R ? n.right() : n)->as_tensor();
      std::tie(occs_[p], virts_[p]) = occ_virt(t);
      ranks_[p] = occs_[p] + virts_[p];
    }

    // no. of contractions in occupied index space (always a whole number)
    occ_ = (occs_[L] + occs_[R] - occs_[T]) / 2;

    // no. of contractions in virtual index space (always a whole number)
    virt_ = (virts_[L] + virts_[R] - virts_[T]) / 2;

    is_outerprod_ = ranks_[L] + ranks_[R] == ranks_[T];
  }

  [[nodiscard]] size_t occ(NodePos p) const noexcept { return occs_[p]; }

  [[nodiscard]] size_t virt(NodePos p) const noexcept { return virts_[p]; }

  [[nodiscard]] size_t rank(NodePos p) const noexcept { return ranks_[p]; }

  [[nodiscard]] size_t occ() const noexcept { return occ_; }

  [[nodiscard]] size_t virt() const noexcept { return virt_; }

  [[nodiscard]] bool is_outerpod() const noexcept { return is_outerprod_; }

  [[nodiscard]] size_t unique_occs() const noexcept {
    return occ(NodePos::Left) + occ(NodePos::Right) - occ();
  }

  [[nodiscard]] size_t unique_virts() const noexcept {
    return virt(NodePos::Left) + virt(NodePos::Right) - virt();
  }

 private:
  std::array<size_t, 3> occs_{0, 0, 0};
  std::array<size_t, 3> virts_{0, 0, 0};
  std::array<size_t, 3> ranks_{0, 0, 0};
  size_t occ_ = 0;
  size_t virt_ = 0;
  bool is_outerprod_ = false;
};
}  // namespace

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
    if (n.leaf()) return AsyCost::zero();
    if (n->op_type() == EvalOp::Product  //
        && n.left()->is_tensor()         //
        && n.right()->is_tensor()) {
      auto const idx_count = ContractedIndexCount{n};
      auto c = AsyCost{idx_count.unique_occs(), idx_count.unique_virts()};
      return idx_count.is_outerpod() ? c : 2 * c;
    } else if (n->is_tensor()) {
      // scalar times a tensor
      // or a tensor plus a tensor
      return AsyCost{occ_virt(n->as_tensor())};
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
      result += expr.is<Tensor>() ? AsyCost{occ_virt(expr.as<Tensor>())}
                                  : AsyCost::zero();
    };

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
    if (n.leaf() || !(n->is_tensor()            //
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
        assert(false &&
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
  auto result = AsyCost::zero();
  auto visitor = [&result](meta::eval_node auto const& n) {
    auto cost = AsyCost::zero();
    if (n.leaf() && n->is_tensor())
      cost = AsyCost{occ_virt(n->as_tensor())};
    else if (!n.leaf()) {
      cost += (n.left()->is_tensor() ? AsyCost{occ_virt(n.left()->as_tensor())}
                                     : AsyCost::zero());
      cost +=
          (n.right()->is_tensor() ? AsyCost{occ_virt(n.right()->as_tensor())}
                                  : AsyCost::zero());
      cost += (n->is_tensor() ? AsyCost{occ_virt(n->as_tensor())}
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

#endif  // SEQUANT_EVAL_NODE_HPP
