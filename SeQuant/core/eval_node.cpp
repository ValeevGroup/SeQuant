#include "eval_node.hpp"
#include "expr.hpp"

#include <boost/math/special_functions/factorials.hpp>

namespace sequant {

EvalNode to_eval_node(ExprPtr const& expr) {
  using ranges::accumulate;
  using ranges::views::tail;
  using ranges::views::transform;

  assert(!expr->is<Constant>() &&
         "constant type expression"
         "not allowed in eval node");

  if (expr->is<Tensor>()) return EvalNode{EvalExpr{expr->as<Tensor>()}};

  auto subxprs = *expr | ranges::views::transform([](auto const& x) {
    return to_eval_node(x);
  }) | ranges::to_vector;

  assert(expr->is<Sum>() || expr->is<Product>());
  auto const op = expr->is<Sum>() ? EvalOp::Sum : EvalOp::Prod;

  auto bnode = ranges::accumulate(
      ranges::views::tail(subxprs), std::move(*subxprs.begin()),
      [op](auto& lnode, auto& rnode) {
        auto pxpr = EvalExpr{*lnode, *rnode, op};
        if (pxpr.op() == EvalOp::Prod) {
          pxpr *= lnode->scalar();
          pxpr *= rnode->scalar();

          lnode->scale(1.0);
          rnode->scale(1.0);
        }

        return EvalNode(std::move(pxpr), std::move(lnode), std::move(rnode));
      });

  if (expr->is<Product>()) *bnode *= expr->as<Product>().scalar();

  return bnode;
}

EvalNode to_eval_node(EvalNode l, EvalNode r, EvalOp op) {
  return EvalNode{EvalExpr{*l, *r, op}, std::move(l), std::move(r)};
}

EvalNode to_eval_node(std::wstring_view label, ExprPtr const& expr) {
  assert(label == L"A" || label == L"S");
  auto op = label == L"A" ? EvalOp::Antisymm : EvalOp::Symm;

  auto rhs = to_eval_node(expr);
  auto bra = rhs->tensor().bra();
  auto ket = rhs->tensor().ket();
  ranges::sort(bra, Index::LabelCompare{});
  ranges::sort(ket, Index::LabelCompare{});

  auto const tnsr = Tensor{label, bra, ket,
                           label == L"A" ? Symmetry::antisymm : Symmetry::symm};
  auto lhs = EvalNode{EvalExpr{tnsr}};
  auto pdata = EvalExpr{*lhs, *rhs, op};
  return EvalNode{std::move(pdata), std::move(lhs), std::move(rhs)};
}

EvalNode to_eval_node_antisymm(ExprPtr const& expr) {
  return to_eval_node(L"A", expr);
}

EvalNode to_eval_node_symm(ExprPtr const& expr) {
  return to_eval_node(L"S", expr);
}

AsyCost asy_cost_single_node_symm_off(EvalNode const& node) {
  if (node.leaf()) return AsyCost::zero();

  auto bks = ranges::views::concat(node.left()->tensor().const_braket(),
                                   node.right()->tensor().const_braket(),
                                   node->tensor().const_braket());
  auto const uniques =
      bks | ranges::to<container::set<Index, Index::LabelCompare>>;

  size_t const nocc = ranges::count_if(uniques, [](auto&& idx) {
    return idx.space() == IndexSpace::active_occupied;
  });

  size_t const nvirt = uniques.size() - nocc;

  switch (node->op()) {
    case EvalOp::Symm: {
      auto f = static_cast<int>(
          boost::math::factorial<double>(node->tensor().rank()));
      return AsyCost{nocc, nvirt, static_cast<int>(f)};
    }
    case EvalOp::Antisymm: {
      auto f = static_cast<int>(
          boost::math::factorial<double>(node->tensor().rank()));
      return AsyCost{nocc, nvirt, static_cast<int>(f * f)};
    }
    default:
      // for matrix multiplication the flops will be doubled
      // to account for the summation ops in a `dot(row, col)` operation
      return AsyCost{nocc, nvirt, node->op() == EvalOp::Prod ? 2 : 1};
  }
}

AsyCost asy_cost_single_node(const EvalNode& node) {
  auto cost = asy_cost_single_node_symm_off(node);
  auto factorial = [](auto x) {
    return static_cast<int>(boost::math::factorial<double>(x));
  };
  // parent node symmetry
  auto const psym = node->tensor().symmetry();
  // parent node bra symmetry
  auto const pbrank = node->tensor().bra_rank();
  // parent node ket symmetry
  auto const pkrank = node->tensor().ket_rank();

  if (psym == Symmetry::nonsymm || psym == Symmetry::invalid) {
    // do nothing
  } else {
    // ------
    // psym is Symmetry::symm or Symmetry::antisymm
    //
    // the rules of cost reduction are taken from
    //   doi:10.1016/j.procs.2012.04.044
    // ------

    auto const op = node->op();
    if (op == EvalOp::Sum) {
      cost = psym == Symmetry::symm
                 ? cost / (factorial(pbrank) * factorial(pkrank))
                 : cost / factorial(pbrank);
    } else if (op == EvalOp::Prod) {
      auto const lsym = node.left()->tensor().symmetry();
      auto const rsym = node.right()->tensor().symmetry();
      cost = (lsym == rsym && lsym == Symmetry::nonsymm)
                 ? cost / factorial(pbrank)
                 : cost / (factorial(pbrank) * factorial(pkrank));
    } else if (op == EvalOp::Symm) {
      cost = cost / factorial(pbrank);
    } else if (op == EvalOp::Antisymm) {
      cost = cost / (factorial(pbrank) * factorial(pkrank));
    } else {
      assert(
          false &&
          "Unsupported evaluation operation for asymptotic cost computation.");
    }
  }

  return cost;
}

}  // namespace sequant
