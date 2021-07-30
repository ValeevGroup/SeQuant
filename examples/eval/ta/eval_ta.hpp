#ifndef SEQUANT_EVAL_TA_HPP
#define SEQUANT_EVAL_TA_HPP

#include "eval.hpp"

#include <tiledarray.h>
#include <range/v3/all.hpp>

namespace sequant::eval {

namespace detail {
auto const braket_to_annot = [](auto const& bk) {
  std::string annot;
  for (auto& idx : bk) {
    annot += idx.string_label() + ",";
  }
  annot.pop_back();
  return annot;
};  // braket_to_annot

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
Tensor_t evaluate_inode_ta(EvalNode const& node, Tensor_t const& leval,
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

} // namespace

template <typename Tensor_t, typename Yielder>
auto eval(EvalNode const& node, Yielder&& yielder, CacheManager<Tensor_t>& man) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

//  Tensor_t (*evaluator)(EvalNode const&, Tensor_t const&, Tensor_t const&)
//       = detail::evaluate_inode_ta;

  auto result = detail::evaluate_single_node(node, std::forward<Yielder>(yielder),
      detail::evaluate_inode_ta<Tensor_t>, man);
  // NOTE:
  // At this point the physical layout of `result`
  // maybe off from what is expected in the residual tensors
  // pre-symmetrization or anti-symmetrization
  //
  // eg.
  //       i_2, i_3, i_1
  // Result
  //       a_1, a_2, a_3
  //
  // we now permute it to the layout:
  //       i_1, i_2, i_3
  // Result
  //       a_1, a_2, a_3
  //
  auto sorted_bra = node->tensor().bra() | ranges::to_vector;
  ranges::sort(sorted_bra, Index::LabelCompare{});
  auto sorted_ket = node->tensor().ket() | ranges::to_vector;
  ranges::sort(sorted_ket, Index::LabelCompare{});

  auto const rannot = detail::braket_to_annot(node->tensor().const_braket());
  auto const lannot = detail::braket_to_annot(ranges::views::concat(sorted_bra,sorted_ket));
  auto scaled = decltype(result){};
  scaled(lannot) = node->scalar().value().real() * result(rannot);
  return scaled;
}

template <typename Tensor_t, typename Yielder>
auto eval_symm(EvalNode const& node, Yielder&& yielder, CacheManager<Tensor_t>& man) {
  auto result = eval(node, std::forward<Yielder>(yielder), man);

  auto symm_result = decltype(result){result.world(), result.trange()};
  symm_result.fill(0);

  auto const lannot =
      ords_to_annot(ranges::views::iota(size_t{0}, result.trange().rank()) |
                    ranges::to_vector);

  auto sym_impl = [&result, &symm_result, &lannot](auto const& perm) {
    symm_result(lannot) += result(ords_to_annot(perm));
  };

  symmetrize_tensor(result.trange().rank(), sym_impl);
  return symm_result;
}

template <typename Tensor_t, typename Yielder>
auto eval_antisymm(EvalNode const& node, Yielder&& yielder, CacheManager<Tensor_t>& man) {
  auto result = eval(node, std::forward<Yielder>(yielder), man);

  auto asymm_result = decltype(result){result.world(), result.trange()};
  asymm_result.fill(0);

  auto const lannot =
      detail::ords_to_annot(ranges::views::iota(size_t{0}, result.trange().rank())
                    | ranges::to_vector);

  auto asym_impl = [&result, &asymm_result,
      &lannot](auto const& pwp) {  // pwp = perm with phase
    asymm_result(lannot) += pwp.phase * result(detail::ords_to_annot(pwp.perm));
  };

  eval::detail::antisymmetrize_tensor(result.trange().rank(), asym_impl);
  return asymm_result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_TA_HPP
