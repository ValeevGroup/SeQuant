#ifndef SEQUANT_EVAL_BTAS_HPP
#define SEQUANT_EVAL_BTAS_HPP

#include <btas/btas.h>
#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <range/v3/all.hpp>

namespace sequant::eval::btas {

namespace detail {
///
/// \param bk iterable of Index objects.
/// \return vector of long-type hash values
///         of the labels of indices in \c bk
///
auto index_hash = [](auto const& bk) {
  return ranges::views::transform(bk,
                                  [](auto const& idx) {
                                    size_t seed = 0;
                                    sequant::hash::combine(seed, idx.label());
                                    return static_cast<long>(seed);  //
                                  }) |
         ranges::to_vector;
};  // index_hash

template <typename Tensor_t>
Tensor_t eval_inode(EvalNode const& node, Tensor_t const& leval,
                    Tensor_t const& reval) {
  assert((node->op() == EvalOp::Sum || node->op() == EvalOp::Prod) &&
         "unsupported intermediate operation");

  auto const assert_imaginary_zero = [](sequant::Constant const& c) {
    assert(c.value().imag() == 0 &&
           "complex scalar unsupported for real tensor");
  };

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  if (node->op() == EvalOp::Prod) {
    auto scalar = node.left()->scalar().value().real() *
                  node.right()->scalar().value().real();

    auto lannot = index_hash(node.left()->tensor().const_braket());
    auto rannot = index_hash(node.right()->tensor().const_braket());
    auto this_annot = index_hash(node->tensor().const_braket());

    Tensor_t prod;

    ::btas::contract(   //
        scalar,         //
        leval, lannot,  //
        reval, rannot,  //
        0.0,            //
        prod, this_annot);
    return prod;
  }

  else {  // sum

    auto const post_annot = index_hash(node->tensor().const_braket());
    auto permute_and_scale = [&post_annot](auto const& btensor,
                                           auto const& child_seqt, auto scal) {
      auto pre_annot = index_hash(child_seqt.const_braket());
      Tensor_t result;
      ::btas::permute(btensor, pre_annot, result, post_annot);
      ::btas::scal(scal, result);
      return result;
    };

    auto lscal = node.left()->scalar().value().real();
    auto rscal = node.right()->scalar().value().real();

    auto sum = permute_and_scale(leval, node.left()->tensor(), lscal);
    sum += permute_and_scale(reval, node.right()->tensor(), rscal);

    return sum;
  }
}

template <typename Tensor_t, typename Yielder>
auto eval_single_node(EvalNode const& node, Yielder&& leaf_evaluator,
                      CacheManager<Tensor_t const>& cache_manager) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

  auto const key = node->hash_value();
  if (auto&& exists = cache_manager.access(key); exists && exists.value())
    return exists.value();
  return node.leaf()
             ? cache_manager.store(key, leaf_evaluator(node->tensor()))
             : cache_manager.store(
                   key,
                   eval_inode(
                       node,
                       *eval_single_node(node.left(),
                                         std::forward<Yielder>(leaf_evaluator),
                                         cache_manager),
                       *eval_single_node(node.right(),
                                         std::forward<Yielder>(leaf_evaluator),
                                         cache_manager)));
}

}  // namespace detail

template <typename Tensor_t, typename Yielder>
auto eval(EvalNode const& node, Yielder&& yielder,
          CacheManager<Tensor_t>& man) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

  auto result =
      *detail::eval_single_node(node, std::forward<Yielder>(yielder), man);
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

  auto const pre_annot = detail::index_hash(node->tensor().const_braket());
  auto const post_annot =
      detail::index_hash(ranges::views::concat(sorted_bra, sorted_ket));

  auto scaled = decltype(result){};
  ::btas::permute(result, pre_annot, scaled, post_annot);
  ::btas::scal(node->scalar().value().real(), scaled);
  return scaled;
}

template <typename Tensor_t, typename Yielder>
auto eval_symm(EvalNode const& node, Yielder&& yielder,
               CacheManager<Tensor_t>& man) {
  auto pre_symm = eval(node, std::forward<Yielder>(yielder), man);
  auto result = decltype(pre_symm){pre_symm.range()};
  result.fill(0.);

  auto const lannot = ranges::views::iota(size_t{0}, pre_symm.rank()) |
                      ranges::to<eval::perm_type>;

  auto symm_impl = [&result, &pre_symm, &lannot](auto const& annot) {
    decltype(result) temp;
    ::btas::permute(pre_symm, lannot, temp, annot);
    result += temp;
  };

  symmetrize_tensor(pre_symm.rank(), symm_impl);
  return result;
}

template <typename Tensor_t, typename Yielder>
auto eval_antisymm(EvalNode const& node, Yielder&& yielder,
                   CacheManager<Tensor_t>& man) {
  auto pre_antisymm = eval(node, std::forward<Yielder>(yielder), man);
  auto result = decltype(pre_antisymm){pre_antisymm.range()};
  result.fill(0.);

  auto const lannot = ranges::views::iota(size_t{0}, pre_antisymm.rank()) |
                      ranges::to<sequant::eval::perm_type>;

  auto asymm_impl = [&result, &pre_antisymm,
                     &lannot](auto const& pwp) {  // pwp = phase with perm
    decltype(result) temp;
    ::btas::permute(pre_antisymm, lannot, temp, pwp.perm);
    ::btas::scal(pwp.phase, temp);
    result += temp;
  };

  antisymmetrize_tensor(pre_antisymm.rank(), asymm_impl);
  return result;
}

}  // namespace sequant::eval::btas

#endif  // SEQUANT_EVAL_BTAS_HPP
