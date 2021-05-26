#ifndef SEQUANT_EVAL_TA_HPP
#define SEQUANT_EVAL_TA_HPP

#include "eval.hpp"

#include <tiledarray.h>
#include <range/v3/all.hpp>

namespace sequant::eval {

auto const braket_to_annot = [](auto const& bk) {
  using ranges::views::join;
  using ranges::views::transform;
  using ranges::views::intersperse;
  return join(bk | transform([](auto const& idx) { return idx.label(); }) |
              intersperse(L",")) |
         ranges::to<std::string>;
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
Tensor_t inode_evaluate_ta(EvalNode const& node, Tensor_t const& leval,
                           Tensor_t const& reval) {
  assert((node->op() == EvalExpr::EvalOp::Sum ||
          node->op() == EvalExpr::EvalOp::Prod) &&
         "unsupported intermediate operation");

  auto assert_imaginary_zero = [](sequant::Constant const& c) {
    assert(c.value().imag() == 0 &&
           "complex scalar unsupported for real tensor");
  };

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  auto this_annot = braket_to_annot(node->tensor().const_braket());
  auto lannot = braket_to_annot(node.left()->tensor().const_braket());
  auto rannot = braket_to_annot(node.right()->tensor().const_braket());

  auto lscal = node.left()->scalar().value().real();
  auto rscal = node.right()->scalar().value().real();

  auto result = Tensor_t{};
  if (node->op() == EvalExpr::EvalOp::Prod) {
    // prod
    result(this_annot) = (lscal * rscal) * leval(lannot) * reval(rannot);
  } else {
    // sum
    result(this_annot) = lscal * leval(lannot) + rscal * reval(rannot);
  }

  return result;
}

template <typename Tensor_t, typename Yielder>
Tensor_t evaluate_ta(EvalNode const& node, Yielder& yielder,
                     sequant::utils::cache_manager<Tensor_t>& cman) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

  auto const key = node->hash();

  if (auto&& exists = cman.access(key); exists && exists.value())
    return *exists.value();

  return node.leaf()
             ? *cman.store(key, yielder(node->tensor()))
             : *cman.store(key,
                           inode_evaluate_ta(
                               node, evaluate_ta(node.left(), yielder, cman),
                               evaluate_ta(node.right(), yielder, cman)));
}

struct eval_instance_ta {
  EvalNode const& node;

  template <typename Tensor_t, typename Fetcher>
  auto evaluate(Fetcher& f, sequant::utils::cache_manager<Tensor_t>& man) {
    static_assert(
        std::is_invocable_r_v<Tensor_t, Fetcher, sequant::Tensor const&>);

    auto result = evaluate_ta(node, f, man);
    auto const annot = braket_to_annot(node->tensor().const_braket());
    auto scaled = decltype(result){};
    scaled(annot) = node->scalar().value().real() * result(annot);
    return scaled;
  }

  template <typename Tensor_t, typename Fetcher>
  auto evaluate_asymm(Fetcher& f,
                      sequant::utils::cache_manager<Tensor_t>& man) {
    auto result = evaluate(f, man);

    auto asymm_result = decltype(result){result.world(), result.trange()};
    asymm_result.fill(0);

    auto const lannot =
        ords_to_annot(ranges::views::iota(size_t{0}, result.trange().rank()) |
                      ranges::to_vector);

    auto asym_impl = [&result, &asymm_result,
                      &lannot](auto const& pwp) {  // pwp = perm with phase
      asymm_result(lannot) += pwp.phase * result(ords_to_annot(pwp.perm));
    };

    sequant::eval::antisymmetrize_tensor(result.trange().rank(), asym_impl);
    return asymm_result;
  }

  template <typename Tensor_t, typename Fetcher>
  auto evaluate_symm(Fetcher& f, sequant::utils::cache_manager<Tensor_t>& man) {
    auto result = evaluate(f, man);

    auto symm_result = decltype(result){result.world(), result.trange()};
    symm_result.fill(0);

    auto const lannot =
        ords_to_annot(ranges::views::iota(size_t{0}, result.trange().rank()) |
                      ranges::to_vector);

    auto sym_impl = [&result, &symm_result, &lannot](auto const& perm) {
      symm_result(lannot) += result(ords_to_annot(perm));
    };

    sequant::eval::symmetrize_tensor(result.trange().rank(), sym_impl);
    return symm_result;
  }
};  // eval_instance_ta

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_TA_HPP
