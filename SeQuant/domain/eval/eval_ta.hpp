#ifndef SEQUANT_EVAL_EVAL_TA_HPP
#define SEQUANT_EVAL_EVAL_TA_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/wstring.hpp>
#include <SeQuant/domain/eval/eval.hpp>

#include <TiledArray/expressions/einsum.h>
#include <TiledArray/expressions/index_list.h>
#include <tiledarray.h>
#include <range/v3/all.hpp>

namespace sequant::eval {

///
/// Given an iterable of Index objects, generate a string annotation
/// that can be used for TiledArray tensor expressions.
/// Tensor-of-tensors also supported.
template <typename Indices>
std::string braket_to_annot(Indices&& indices) {
  using ranges::find;
  using ranges::views::filter;
  using ranges::views::intersperse;
  using ranges::views::join;

  // make a comma-separated string out of an iterable of strings
  auto add_commas = [](auto const& strs) -> std::string {
    return strs | intersperse(",") | join | ranges::to<std::string>;
  };

  container::svector<std::string> idxs{}, pidxs{};
  for (auto&& idx : indices) {
    idxs.emplace_back(sequant::to_string(idx.label()));
    for (auto&& pidx : idx.proto_indices())
      pidxs.emplace_back(sequant::to_string(pidx.label()));
  }

  if (pidxs.empty()) {
    // not a tensor-of-tensor type expression
    return add_commas(idxs);
  } else {
    ranges::stable_sort(pidxs);
    ranges::actions::unique(pidxs);
    auto not_in_pidx = [&pidxs](auto&& l) {
      return find(pidxs, l) == pidxs.end();
    };
    return add_commas(pidxs) + ";" +
           add_commas(idxs | filter(not_in_pidx) | ranges::to<decltype(idxs)>);
  }
}

class EvalExprTA final : public EvalExpr {
 public:
  ///
  /// annotation for TiledArray
  ///
  [[nodiscard]] std::string const& annot() const;

  ///
  /// Whether this object represents a tensor-of-tensor kind expression
  ///
  [[nodiscard]] bool tot() const;

  explicit EvalExprTA(Tensor const&);

  EvalExprTA(EvalExprTA const&, EvalExprTA const&, EvalOp);

 private:
  std::string annot_;

  bool tot_;
};

using EvalNodeTA = FullBinaryNode<EvalExprTA>;

template <typename DistArrayTot, typename DistArray>
using tot_result_t = std::variant<DistArrayTot, DistArray>;

EvalNodeTA to_eval_node_ta(ExprPtr const& expr);

EvalNodeTA to_eval_node_ta(EvalNode const& node);

namespace detail {

template <typename RngOfOrdinals>
auto ords_to_annot(RngOfOrdinals const& ords) {
  using ranges::views::intersperse;
  using ranges::views::join;
  using ranges::views::transform;
  auto to_str = [](auto x) { return std::to_string(x); };
  return ords | transform(to_str) | intersperse(std::string{","}) | join |
         ranges::to<std::string>;
}

template <typename Tensor_t>
Tensor_t eval_inode(EvalNodeTA const& node, Tensor_t const& leval,
                    Tensor_t const& reval) {
  assert((node->op() == EvalOp::Sum || node->op() == EvalOp::Prod) &&
         "unsupported intermediate operation");

  auto assert_imaginary_zero = [](sequant::Constant const& c) {
    assert(c.value().imag() == 0 && "complex scalar unsupported");
  };

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  auto const lscal = node.left()->scalar().value().real();
  auto const rscal = node.right()->scalar().value().real();
  auto const& this_annot = node->annot();
  auto const& lannot = node.left()->annot();
  auto const& rannot = node.right()->annot();

  Tensor_t result;
  {
    if (node->op() == EvalOp::Prod) {
      // prod
      result(this_annot) = (lscal * rscal) * leval(lannot) * reval(rannot);
    } else {
      // sum
      result(this_annot) = lscal * leval(lannot) + rscal * reval(rannot);
    }
  }
  Tensor_t::wait_for_lazy_cleanup(result.world());

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "eval_inode: evaluated "
            << to_string(to_latex_align(to_expr(node)))
            << " worldobj.id=" << result.id();
#ifdef TA_TENSOR_MEM_PROFILE
  std::cout << " TA::Tensor allocated {"
            << "hw="
            << TA::hostEnv::instance()->host_allocator_getActualHighWatermark()
            << ","
            << "cur="
            << TA::hostEnv::instance()->host_allocator().getCurrentSize() << ","
            << "act="
            << TA::hostEnv::instance()->host_allocator().getActualSize() << "}"
            << " bytes";
#endif  // TA_TENSOR_MEM_PROFILE
  std::cout << std::endl;
#endif  // SEQUANT_EVAL_TRACE

  return result;
}

template <typename Tensor_t, typename Yielder>
Tensor_t eval_single_node(EvalNodeTA const& node, Yielder&& leaf_evaluator,
                          CacheManager<Tensor_t const>& cache_manager) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

  auto const key = node->hash_value();

  if (auto&& exists = cache_manager.access(key); exists && exists.value())
    return *exists.value();

  Tensor_t result;
  {
    Tensor_t value;
    if (node.leaf())
      value = leaf_evaluator(node->tensor());
    else {
      Tensor_t leval = eval_single_node(
          node.left(), std::forward<Yielder>(leaf_evaluator), cache_manager);

      Tensor_t reval = eval_single_node(
          node.right(), std::forward<Yielder>(leaf_evaluator), cache_manager);

      value = eval_inode(node, std::move(leval), std::move(reval));
    }
    result = *cache_manager.store(key, std::move(value));
  }
  Tensor_t::wait_for_lazy_cleanup(result.world());

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "eval_single_node: evaluated "
            << to_string(to_latex_align(to_expr(node)))
            << " worldobj.id=" << result.id();
#ifdef TA_TENSOR_MEM_PROFILE
  std::cout << " TA::Tensor allocated {"
            << "hw="
            << TA::hostEnv::instance()->host_allocator_getActualHighWatermark()
            << ","
            << "cur="
            << TA::hostEnv::instance()->host_allocator().getCurrentSize() << ","
            << "act="
            << TA::hostEnv::instance()->host_allocator().getActualSize() << "}"
            << " bytes";
#endif  // TA_TENSOR_MEM_PROFILE
  std::cout << std::endl;
#endif  // SEQUANT_EVAL_TRACE

  return result;
}

}  // namespace detail

///
///
/// Evaluate expressions from left to right as they appear in @c Expr.
/// The @c Expr corresponding a @c EvalNodeTA can be generated using @c to_expr
/// function.
///
/// A node a full-binary tree @see EvalNodeTA. The evaluation
/// occurs in a left-to-right order. This can also be thought of in terms of
/// post-order traversal (a binary tree traversal technique, see Wikipedia).
/// If the current EvalNodeTA object has not yet been evaluated earlier (tested
/// by checking the cache), then the left node is first evaluated, followed
/// by the right node. The result of the left and the right evaluations are
/// then used to evaluate this node.
///
/// \tparam Tensor_t The data tensor type. eg. TA::TArrayD from TiledArray.
/// \param node sequant::binary_node<sequant::EvalExpr> object
///            that is the evaluation tree in essence.
/// \param target_indx_labels The string labels iterable of the target indices
///                           that appear in the node->tensor().const_braket().
///                               set(node->tensor().const_braket()) =
///                                         set(target_indx_labels)
///                           This paramter is to allow to set the final
///                           physical layout of the evaluated tensor in
///                           the desired permutation.
/// \param yielder That returns Tensor_t for leaf SeQuant Tensor(g, f, t, ...).
/// \param man The cache manager.
/// \return Tensor_t
/// @see @c to_expr
///
template <typename Tensor_t, typename Iterable, typename Yielder>
auto eval(EvalNodeTA const& node, Iterable const& target_indx_labels,
          Yielder&& yielder, CacheManager<Tensor_t const>& man) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

#ifndef NDEBUG
  auto ti_sorted_input =
      target_indx_labels | ranges::to<container::svector<std::string>>;
  ranges::sort(ti_sorted_input);
  auto ti_sorted_node = node->tensor().const_braket() |
                        ranges::views::transform(
                            [](auto const& idx) { return idx.to_string(); }) |
                        ranges::to<container::svector<std::string>>;
  ranges::sort(ti_sorted_node);

  assert(ti_sorted_input == ti_sorted_node && "Invalid target indices");
#endif

  Tensor_t scaled;
  {
    auto result =
        detail::eval_single_node(node, std::forward<Yielder>(yielder), man);

    std::string lannot = ranges::front(target_indx_labels);
    for (auto const& lbl : ranges::views::tail(target_indx_labels))
      lannot += std::string{','} + lbl;

    scaled = decltype(result){};
    scaled(lannot) = node->scalar().value().real() * result(node->annot());
  }
  Tensor_t::wait_for_lazy_cleanup(scaled.world());

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "eval: evaluated " << to_string(to_latex_align(to_expr(node)))
            << " worldobj.id=" << scaled.id();
#ifdef TA_TENSOR_MEM_PROFILE
  std::cout << " TA::Tensor allocated {"
            << "hw="
            << TA::hostEnv::instance()->host_allocator_getActualHighWatermark()
            << ","
            << "cur="
            << TA::hostEnv::instance()->host_allocator().getCurrentSize() << ","
            << "act="
            << TA::hostEnv::instance()->host_allocator().getActualSize() << "}"
            << " bytes";
#endif  // TA_TENSOR_MEM_PROFILE
  std::cout << std::endl;
#endif  // SEQUANT_EVAL_TRACE

  return scaled;
}

///
/// Evaluate expressions from left to right as they appear in @c Expr and
/// particle-symmetrize the result.
///
/// \tparam Tensor_t The data tensor type. eg. TA::TArrayD from TiledArray.
/// \param node sequant::binary_node<sequant::EvalExpr> object
///            that is the evaluation tree in essence.
/// \param target_indx_labels The string labels iterable of the target indices
///                           that appear in the node->tensor().const_braket().
///                               set(node->tensor().const_braket()) =
///                                         set(target_indx_labels)
///                           This paramter is to allow to set the final
///                           physical layout of the evaluated tensor in
///                           the desired permutation.
/// \param yielder That returns Tensor_t for leaf SeQuant Tensor(g, f, t, ...).
/// \param man The cache manager.
/// \return Tensor_t
/// @see @c eval
///
template <typename Tensor_t, typename Iterable, typename Yielder>
auto eval_symm(EvalNodeTA const& node, Iterable const& target_indx_labels,
               Yielder&& yielder, CacheManager<Tensor_t const>& man) {
  Tensor_t symm_result;
  {
    auto result =
        eval(node, target_indx_labels, std::forward<Yielder>(yielder), man);

    symm_result = decltype(result){result.world(), result.trange()};
    symm_result.fill(0);

    auto const lannot = detail::ords_to_annot(
        ranges::views::iota(size_t{0}, result.trange().rank()) |
        ranges::to_vector);

    auto sym_impl = [&result, &symm_result, &lannot](auto const& perm) {
      symm_result(lannot) += result(detail::ords_to_annot(perm));
      Tensor_t::wait_for_lazy_cleanup(symm_result.world());
    };

    symmetrize_tensor(result.trange().rank(), sym_impl);
  }
  Tensor_t::wait_for_lazy_cleanup(symm_result.world());

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "eval_symm: evaluated "
            << to_string(to_latex_align(to_expr(node)))
            << " worldobj.id=" << symm_result.id();
#ifdef TA_TENSOR_MEM_PROFILE
  std::cout << " TA::Tensor allocated {"
            << "hw="
            << TA::hostEnv::instance()->host_allocator_getActualHighWatermark()
            << ","
            << "cur="
            << TA::hostEnv::instance()->host_allocator().getCurrentSize() << ","
            << "act="
            << TA::hostEnv::instance()->host_allocator().getActualSize() << "}"
            << " bytes";
#endif  // TA_TENSOR_MEM_PROFILE
  std::cout << std::endl;
#endif  // SEQUANT_EVAL_TRACE

  return symm_result;
}

///
/// Evaluate expressions from left to right as they appear in @c Expr and
/// particle-antisymmetrize the result.
///
/// \tparam Tensor_t The data tensor type. eg. TA::TArrayD from TiledArray.
/// \param node sequant::binary_node<sequant::EvalExpr> object
///            that is the evaluation tree in essence.
/// \param target_indx_labels The string labels iterable of the target indices
///                           that appear in the node->tensor().const_braket().
///                               set(node->tensor().const_braket()) =
///                                         set(target_indx_labels)
///                           This paramter is to allow to set the final
///                           physical layout of the evaluated tensor in
///                           the desired permutation.
/// \param yielder That returns Tensor_t for leaf SeQuant Tensor(g, f, t, ...).
/// \param man The cache manager.
/// \return Tensor_t
/// @see @c eval
///
template <typename Tensor_t, typename Iterable, typename Yielder>
auto eval_antisymm(EvalNodeTA const& node, Iterable const& target_indx_labels,
                   Yielder&& yielder, CacheManager<Tensor_t const>& man) {
  Tensor_t antisymm_result;
  {
    auto result =
        eval(node, target_indx_labels, std::forward<Yielder>(yielder), man);

    antisymm_result = decltype(result){result.world(), result.trange()};
    antisymm_result.fill(0);

    auto const lannot = detail::ords_to_annot(
        ranges::views::iota(size_t{0}, result.trange().rank()) |
        ranges::to_vector);

    auto asym_impl = [&result, &antisymm_result,
                      &lannot](auto const& pwp) {  // pwp = perm with phase
      antisymm_result(lannot) +=
          pwp.phase * result(detail::ords_to_annot(pwp.perm));
      Tensor_t::wait_for_lazy_cleanup(antisymm_result.world());
    };

    antisymmetrize_tensor(result.trange().rank(), asym_impl);
  }
  Tensor_t::wait_for_lazy_cleanup(antisymm_result.world());

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "eval_antisymm: evaluated "
            << to_string(to_latex_align(to_expr(node)))
            << " worldobj.id=" << antisymm_result.id();
#ifdef TA_TENSOR_MEM_PROFILE
  std::cout << " TA::Tensor allocated {"
            << "hw="
            << TA::hostEnv::instance()->host_allocator_getActualHighWatermark()
            << ","
            << "cur="
            << TA::hostEnv::instance()->host_allocator().getCurrentSize() << ","
            << "act="
            << TA::hostEnv::instance()->host_allocator().getActualSize() << "}"
            << " bytes";
#endif  // TA_TENSOR_MEM_PROFILE
  std::cout << std::endl;
#endif  // SEQUANT_EVAL_TRACE

  return antisymm_result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EVAL_TA_HPP