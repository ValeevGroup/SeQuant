#ifndef SEQUANT_EVAL_EVAL_TA_HPP
#define SEQUANT_EVAL_EVAL_TA_HPP

#include <SeQuant/domain/eval/eval.hpp>

#include <tiledarray.h>
#include <TiledArray/expressions/index_list.h>
#include <TiledArray/expressions/einsum.h>
#include <range/v3/all.hpp>

namespace sequant::eval {

namespace detail {

auto const ords_to_annot = [](auto const& ords) {
  using ranges::accumulate;
  using ranges::views::intersperse;
  using ranges::views::transform;
  auto to_str = [](auto x) { return std::to_string(x); };
  return ranges::accumulate(
      ords | transform(to_str) | intersperse(std::string{","}), std::string{},
      std::plus{});
};  // ords_to_annot

template <typename Tensor_t>
Tensor_t eval_inode(EvalNode const& node, Tensor_t const& leval,
                    Tensor_t const& reval) {
  assert((node->op() == EvalOp::Sum || node->op() == EvalOp::Prod) &&
         "unsupported intermediate operation");

  auto assert_imaginary_zero = [](sequant::Constant const& c) {
    assert(c.value().imag() == 0 &&
           "complex scalar unsupported");
  };

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  auto const lscal = node.left()->scalar().value().real();
  auto const rscal = node.right()->scalar().value().real();
  auto const& this_annot = node->annot();
  auto const& lannot = node.left()->annot();
  auto const& rannot = node.right()->annot();

  auto result = Tensor_t{};
  if (node->op() == EvalOp::Prod) {
    // prod
    result(this_annot) = (lscal * rscal) * leval(lannot) * reval(rannot);
  } else {
    // sum
    result(this_annot) = lscal * leval(lannot) + rscal * reval(rannot);
  }

  TA::get_default_world().gop.fence();
  return result;
}

template <typename Tensor_t>
Tensor_t eval_inode_tot(EvalNode const& node, Tensor_t const& leval,
                        Tensor_t const& reval) {
  assert((node->op() == EvalOp::Sum || node->op() == EvalOp::Prod) &&
         "unsupported intermediate operation");

  auto assert_imaginary_zero = [](sequant::Constant const& c) {
    assert(c.value().imag() == 0 &&
           "complex scalar unsupported");
  };

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  auto const lscal = node.left()->scalar().value().real();
  auto const rscal = node.right()->scalar().value().real();
  auto const& this_annot = node->annot();
  auto const& lannot = node.left()->annot();
  auto const& rannot = node.right()->annot();

  auto result = Tensor_t{};
  if (node->op() == EvalOp::Prod) {
    // prod
    // result(this_annot) = (lscal * rscal) * leval(lannot) * reval(rannot);
    decltype(result) unscaled = TA::expressions::einsum(leval(lannot),
                                                        reval(rannot), this_annot);
    result(this_annot) = (lscal * rscal) * unscaled(this_annot);
  } else {
    // sum
    result(this_annot) = lscal * leval(lannot) + rscal * reval(rannot);
  }
  TA::get_default_world().gop.fence();
  return result;
}

template <typename Tensor_t, typename Yielder>
Tensor_t eval_single_node(EvalNode const& node, Yielder&& leaf_evaluator,
                          CacheManager<Tensor_t const>& cache_manager) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

  auto const key = node->hash();

  if (auto&& exists = cache_manager.access(key); exists && exists.value())
    return *exists.value();

  return node.leaf()
             ? *cache_manager.store(key, leaf_evaluator(node->tensor()))
             : *cache_manager.store(
                   key,
                   eval_inode(
                       node,
                       eval_single_node(node.left(),
                                        std::forward<Yielder>(leaf_evaluator),
                                        cache_manager),
                       eval_single_node(node.right(),
                                        std::forward<Yielder>(leaf_evaluator),
                                        cache_manager)));
}

template <typename Tensor_t, typename Yielder>
Tensor_t eval_single_node_tot(EvalNode const& node, Yielder&& leaf_evaluator,
                          CacheManager<Tensor_t const>& cache_manager) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

  auto const key = node->hash();

  if (auto&& exists = cache_manager.access(key); exists && exists.value())
    return *exists.value();

  return node.leaf()
             ? *cache_manager.store(key, leaf_evaluator(node->tensor()))
             : *cache_manager.store(
                   key,
                   eval_inode_tot(
                       node,
                       eval_single_node_tot(node.left(),
                                        std::forward<Yielder>(leaf_evaluator),
                                        cache_manager),
                       eval_single_node_tot(node.right(),
                                        std::forward<Yielder>(leaf_evaluator),
                                        cache_manager)));
}

}  // namespace detail

///
/// Evaluate a node.
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
template <typename Tensor_t, typename Iterable, typename Yielder>
auto eval(EvalNode const& node, Iterable const& target_indx_labels,
          Yielder&& yielder, CacheManager<Tensor_t const>& man) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

  auto ti_sorted_input =
      target_indx_labels | ranges::to<container::svector<std::string>>;
  ranges::sort(ti_sorted_input);
  auto ti_sorted_node = node->tensor().const_braket() |
                        ranges::views::transform([](auto const& idx) {
                          return idx.ascii_label();
                        }) |
                        ranges::to<container::svector<std::string>>;
  ranges::sort(ti_sorted_node);

  assert(ti_sorted_input == ti_sorted_node && "Invalid target indices");

  auto result =
      detail::eval_single_node(node, std::forward<Yielder>(yielder), man);

  std::string lannot = ranges::front(target_indx_labels);
  for (auto const& lbl : ranges::views::tail(target_indx_labels))
    lannot += std::string{','} + lbl;

  auto scaled = decltype(result){};
  scaled(lannot) = node->scalar().value().real() * result(node->annot());
  TA::get_default_world().gop.fence();
  return scaled;
}

template <typename Tensor_t,
          typename Iterable1,
          typename Iterable2,
          typename Yielder,
          std::enable_if_t<std::is_invocable_r_v<Tensor_t,
                                                 Yielder,
                                                 sequant::Tensor const&>,
                           bool> = true>
Tensor_t eval_tot(EvalNode const& node,
          Iterable1 const& outer_indx_labels,
          Iterable2 const& inner_indx_labels,
          Yielder&& yielder, CacheManager<Tensor_t const>& man) {
    auto bpindx_rcvd = TA::expressions::BipartiteIndexList{outer_indx_labels, inner_indx_labels};
    auto bpindx_exst = TA::expressions::BipartiteIndexList{node->annot()};
    assert(bpindx_exst.first().is_permutation(bpindx_exst.first())
            && bpindx_exst.second().is_permutation(bpindx_exst.second())
            && "Invalid target index labels");

    auto result = detail::eval_single_node_tot(node, std::forward<Yielder>(yielder), man);

    auto const lannot = bpindx_rcvd.first().string()
                        + ";"
                        + bpindx_rcvd.second().string();
    auto scaled = decltype(result){};
    scaled(lannot) = node->scalar().value().real() * result(node->annot());
    TA::get_default_world().gop.fence();
    return scaled;
}

///
/// Evaluate a node and symmetrize the result.
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
template <typename Tensor_t, typename Iterable, typename Yielder>
auto eval_symm(EvalNode const& node, Iterable const& target_indx_labels,
               Yielder&& yielder, CacheManager<Tensor_t const>& man) {
  auto result =
      eval(node, target_indx_labels, std::forward<Yielder>(yielder), man);

  auto symm_result = decltype(result){result.world(), result.trange()};
  symm_result.fill(0);

  auto const lannot = detail::ords_to_annot(
      ranges::views::iota(size_t{0}, result.trange().rank()) |
      ranges::to_vector);

  auto sym_impl = [&result, &symm_result, &lannot](auto const& perm) {
    symm_result(lannot) += result(detail::ords_to_annot(perm));
  };

  symmetrize_tensor(result.trange().rank(), sym_impl);
  TA::get_default_world().gop.fence();
  return symm_result;
}

///
/// Evaluate a node and anit-symmetrize the result.
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
template <typename Tensor_t, typename Iterable, typename Yielder>
auto eval_antisymm(EvalNode const& node, Iterable const& target_indx_labels,
                   Yielder&& yielder, CacheManager<Tensor_t const>& man) {
  auto result =
      eval(node, target_indx_labels, std::forward<Yielder>(yielder), man);

  auto antisymm_result = decltype(result){result.world(), result.trange()};
  antisymm_result.fill(0);

  auto const lannot = detail::ords_to_annot(
      ranges::views::iota(size_t{0}, result.trange().rank()) |
      ranges::to_vector);

  auto asym_impl = [&result, &antisymm_result,
                    &lannot](auto const& pwp) {  // pwp = perm with phase
    antisymm_result(lannot) += pwp.phase * result(detail::ords_to_annot(pwp.perm));
  };

  antisymmetrize_tensor(result.trange().rank(), asym_impl);
  TA::get_default_world().gop.fence();
  return antisymm_result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EVAL_TA_HPP
