#ifndef SEQUANT_EVAL_EVAL_TA_HPP
#define SEQUANT_EVAL_EVAL_TA_HPP

#include <SeQuant/domain/eval/eval.hpp>

#include <tiledarray.h>
#include <range/v3/all.hpp>

namespace sequant::eval::ta {

namespace detail {

///
/// Given an iterable of Index objects, generate a string annotation
/// that can be used for TiledArray tensor expressions. Tensor-of-tensors also
/// supported.
template <typename Indices_t>
std::string braket_to_annot(Indices_t const& indices) {
  // make a comma-separated and concatenated string out of an iterable of strings
  auto add_commas = [](auto const& strs) -> std::string {
    std::string result{ranges::front(strs)};
    for (auto&& s: ranges::views::tail(strs))
      result += "," + s;
    return result;
  };

  container::vector<std::string> outer_labels{}, inner_labels{};
  for (auto&& idx: indices) {
    inner_labels.emplace_back(idx.ascii_label());
    for (auto&& pidx: idx.proto_indices())
      outer_labels.emplace_back(pidx.ascii_label());
  }

  if (outer_labels.empty())
    return add_commas(inner_labels);

  // support CSV methods
  ranges::sort(outer_labels);
  ranges::sort(inner_labels);
  auto outer_labels_updated = ranges::views::set_difference(outer_labels, inner_labels);
  return add_commas(outer_labels_updated) + ";" + add_commas(inner_labels);
}

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
           "complex scalar unsupported for real tensor");
  };

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  auto const this_annot = braket_to_annot(node->tensor().const_braket());
  auto const lannot = braket_to_annot(node.left()->tensor().const_braket());
  auto const rannot = braket_to_annot(node.right()->tensor().const_braket());

  auto const lscal = node.left()->scalar().value().real();
  auto const rscal = node.right()->scalar().value().real();

  auto result = Tensor_t{};
  if (node->op() == EvalOp::Prod) {
    // prod
    result(this_annot) = (lscal * rscal) * leval(lannot) * reval(rannot);
  } else {
    // sum
    assert(node->op() == EvalOp::Sum && "unsupported operation for eval");
    result(this_annot) = lscal * leval(lannot) + rscal * reval(rannot);
  }

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

  auto const rannot = detail::braket_to_annot(node->tensor().const_braket());

  std::string lannot = ranges::front(target_indx_labels);
  for (auto const& lbl : ranges::views::tail(target_indx_labels))
    lannot += std::string{','} + lbl;

  auto scaled = decltype(result){};
  scaled(lannot) = node->scalar().value().real() * result(rannot);
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
  return antisymm_result;
}

}  // namespace sequant::eval::ta

#endif  // SEQUANT_EVAL_EVAL_TA_HPP
