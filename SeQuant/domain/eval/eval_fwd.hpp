//
// Created by Bimal Gaudel on 3/27/25.
//

#ifndef SEQUANT_EVAL_FWD_HPP
#define SEQUANT_EVAL_FWD_HPP

#include <SeQuant/core/eval_expr.hpp>

namespace sequant {

class CacheManager;
class Result;

///
/// \brief Managed pointer to the result of an evaluation.
///
using ResultPtr = std::shared_ptr<Result>;

namespace meta {

///
/// \brief Satisfied by a type with a method named \code annot that returns
///        a non-void type.
///
template <typename T>
concept has_annot = requires(T t) {
  t.annot();
  requires !std::is_void_v<decltype(t.annot())>;
};

///
/// \brief Satisfied by an eval_node whose dereferenced type satisfies the
///        has_annot method.
/// \example
///          * \code static_assert(!meta::can_evaluate<EvalNode<EvalExpr>>);
///          * \code static_assert(meta::can_evaluate<EvalNode<EvalExprTA>>);
///
template <typename T>
concept can_evaluate = eval_node<T> && requires(T n) {
  { *n } -> has_annot;
};

///
/// \brief Satisfied by a range type of objects satisfying can_evaluate.
///
template <typename Rng>
concept can_evaluate_range =
    std::ranges::range<Rng> && can_evaluate<std::ranges::range_value_t<Rng>>;

///
/// \brief \tparam F is a leaf node evaluator of type \tparam Node if
///        an object (a function object) of type \tparam F returns ResultPtr
///        when called with the single argument of const ref type to
///        \tparam Node and the \tparam Node satisfies can_evaluate.
///
template <typename Node, typename F>
concept leaf_node_evaluator =
    can_evaluate<Node> && requires(F f, Node const& n) {
      { f(n) } -> std::same_as<ResultPtr>;
    };
}  // namespace meta

static_assert(meta::eval_node<EvalNode<EvalExpr>>);
static_assert(meta::eval_node<EvalNode<EvalExprTA>>);
static_assert(meta::eval_node<EvalNode<EvalExprBTAS>>);

static_assert(!meta::can_evaluate<EvalNode<EvalExpr>>);
static_assert(meta::can_evaluate<EvalNode<EvalExprTA>>);
static_assert(meta::can_evaluate<EvalNode<EvalExprBTAS>>);

}  // namespace sequant
#endif  // SEQUANT_EVAL_FWD_HPP
