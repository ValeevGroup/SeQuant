#ifndef SEQUANT_EVAL_BTAS_HPP
#define SEQUANT_EVAL_BTAS_HPP

#include "eval.hpp"

#include <btas/btas.h>
#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/tensor.hpp>
#include <range/v3/all.hpp>

namespace sequant::eval {

template <typename Tensor_t>
Tensor_t inode_evaluate_btas(EvalNode const& node, Tensor_t const& leval,
                             Tensor_t const& reval) {
  assert((node->op() == EvalOp::Sum ||
          node->op() == EvalOp::Prod) &&
         "unsupported intermediate operation");

  auto const assert_imaginary_zero = [](sequant::Constant const& c) {
    assert(c.value().imag() == 0 &&
           "complex scalar unsupported for real tensor");
  };

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  auto index_hash = [](auto const& bk) {
    return ranges::views::transform(bk,
                                    [](auto const& idx) {
                                      size_t seed = 0;
                                      sequant::hash::combine(seed, idx.label());
                                      return static_cast<long>(seed);  //
                                    }) |
           ranges::to_vector;
  };  // index_hash

  if (node->op() == EvalOp::Prod) {
    auto scalar = node.left()->scalar().value().real() *
                  node.right()->scalar().value().real();

    auto lannot = index_hash(node.left()->tensor().const_braket());
    auto rannot = index_hash(node.right()->tensor().const_braket());
    auto this_annot = index_hash(node->tensor().const_braket());

    Tensor_t prod;

    btas::contract(     //
        scalar,         //
        leval, lannot,  //
        reval, rannot,  //
        0.0,            //
        prod, this_annot);
    return prod;
  }

  else {  // sum

    auto const post_annot = index_hash(node->tensor().const_braket());
    auto permute_and_scale = [&post_annot, &index_hash](auto const& btensor,
                                                        auto const& child_seqt,
                                                        auto scal) {
      auto pre_annot = index_hash(child_seqt.const_braket());
      Tensor_t result;
      btas::permute(btensor, pre_annot, result, post_annot);
      btas::scal(scal, result);
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
Tensor_t evaluate_btas(EvalNode const& node, Yielder& yielder,
                       sequant::utils::cache_manager<Tensor_t>& cman) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);
  auto const key = node->hash();

  auto exists = cman.access(key);
  if (exists && exists.value()) return *exists.value();

  return node.leaf()
             ? *cman.store(key, yielder(node->tensor()))
             : *cman.store(key,
                           inode_evaluate_btas(
                               node, evaluate_btas(node.left(), yielder, cman),
                               evaluate_btas(node.right(), yielder, cman)));
}

struct eval_instance_btas {
  EvalNode const& node;

  template <typename Tensor_t, typename Fetcher>
  auto evaluate(Fetcher& yielder,
                sequant::utils::cache_manager<Tensor_t>& cman) {
    static_assert(
        std::is_invocable_r_v<Tensor_t, Fetcher, sequant::Tensor const&>);
    auto result = evaluate_btas(node, yielder, cman);

    btas::scal(node->scalar().value().real(), result);

    return result;
  }

  template <typename Tensor_t, typename Fetcher>
  auto evaluate_symm(Fetcher& yielder,
                     sequant::utils::cache_manager<Tensor_t>& cman) {
    auto pre_symm = evaluate(yielder, cman);
    auto result = decltype(pre_symm){pre_symm.range()};
    result.fill(0.);

    auto const lannot = ranges::views::iota(size_t{0}, pre_symm.rank()) |
                        ranges::to<sequant::eval::perm_type>;

    auto symm_impl = [&result, &pre_symm, &lannot](auto const& annot) {
      decltype(result) temp;
      btas::permute(pre_symm, lannot, temp, annot);
      result += temp;
    };

    sequant::eval::symmetrize_tensor(pre_symm.rank(), symm_impl);
    return result;
  }

  template <typename Tensor_t, typename Fetcher>
  auto evaluate_asymm(Fetcher& yielder,
                      sequant::utils::cache_manager<Tensor_t>& cman) {
    auto pre_asymm = evaluate(yielder, cman);
    auto result = decltype(pre_asymm){pre_asymm.range()};
    result.fill(0.);

    auto const lannot = ranges::views::iota(size_t{0}, pre_asymm.rank()) |
                        ranges::to<sequant::eval::perm_type>;

    auto asymm_impl = [&result, &pre_asymm,
                       &lannot](auto const& pwp) {  // pwp = phase with perm
      decltype(result) temp;
      btas::permute(pre_asymm, lannot, temp, pwp.perm);
      btas::scal(pwp.phase, temp);
      result += temp;
    };

    sequant::eval::antisymmetrize_tensor(pre_asymm.rank(), asymm_impl);
    return result;
  }

};  // eval_instance_btas

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_BTAS_HPP
