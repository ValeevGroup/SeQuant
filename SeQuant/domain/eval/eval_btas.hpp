#ifndef SEQUANT_EVAL_BTAS_HPP
#define SEQUANT_EVAL_BTAS_HPP

#include "eval.hpp"

#include <btas/btas.h>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <range/v3/all.hpp>

namespace sequant::eval {

template <typename Tensor_t>
class yield_leaf {
 private:
  Tensor_t const &fock, &eri, &t_vo, &t_vvoo;

  size_t const nocc, nvirt;

  auto range1_limits(sequant::Tensor const& tensor) {
    return tensor.const_braket() |
           ranges::views::transform([this](auto const& idx) {
             auto ao = sequant::IndexSpace::active_occupied;
             auto au = sequant::IndexSpace::active_unoccupied;
             auto sp = idx.space();
             assert(sp == ao || sp == au);

             return sp == ao ? nocc : nvirt;
           });
  }

 public:
  yield_leaf(size_t no, size_t nv, Tensor_t const& F, Tensor_t const& G,
             Tensor_t const& ampl_vo, Tensor_t const& ampl_vvoo)
      : nocc{no},
        nvirt{nv},
        fock{F},
        eri{G},
        t_vo{ampl_vo},
        t_vvoo{ampl_vvoo} {}

  Tensor_t operator()(sequant::Tensor const& texpr) {
    auto rank = texpr.bra_rank() + texpr.ket_rank();

    if (texpr.label() == L"t") {
      assert(rank == 2 || rank == 4 && "only t_vo and t_vvoo supported");
      return rank == 2 ? t_vo : t_vvoo;
    }

    assert((texpr.label() == L"g" || texpr.label() == L"f") &&
           "unsupported tensor label encountered");

    auto&& big_tensor = texpr.label() == L"g" ? eri : fock;

    auto r1_limits = range1_limits(texpr);
    auto iter_limits = r1_limits | ranges::views::transform([this](auto x) {
                         return x == nocc ? std::pair{size_t{0}, nocc}
                                          : std::pair{nocc, nocc + nvirt};
                       });

    auto slice = Tensor_t{btas::Range{r1_limits | ranges::to_vector}};

    if (iter_limits.size() == 2) {
      auto loop1 = iter_limits[0];
      auto loop2 = iter_limits[1];
      for (auto i = loop1.first; i < loop1.second; ++i)
        for (auto j = loop2.first; j < loop2.second; ++j)
          slice(i - loop1.first, j - loop2.first) = big_tensor(i, j);

    } else {  // iter_limits.size() == 4 true
      auto loop1 = iter_limits[0];
      auto loop2 = iter_limits[1];
      auto loop3 = iter_limits[2];
      auto loop4 = iter_limits[3];
      for (auto i = loop1.first; i < loop1.second; ++i)
        for (auto j = loop2.first; j < loop2.second; ++j)
          for (auto k = loop3.first; k < loop3.second; ++k)
            for (auto l = loop4.first; l < loop4.second; ++l)
              slice(i - loop1.first,    //
                    j - loop2.first,    //
                    k - loop3.first,    //
                    l - loop4.first) =  //
                  big_tensor(i, j, k, l);
    }

    return slice;
  }
};  // tensor yield

template <typename Tensor_t>
Tensor_t inode_evaluate_btas(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    Tensor_t const& leval, Tensor_t const& reval) {
  assert((node->op() == sequant::utils::eval_expr::eval_op::Sum ||
          node->op() == sequant::utils::eval_expr::eval_op::Prod) &&
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

  if (node->op() == sequant::utils::eval_expr::eval_op::Prod) {
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
Tensor_t evaluate_btas(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    Yielder& yielder, sequant::utils::cache_manager<Tensor_t>& cman) {
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

struct eval_instance {
  sequant::utils::binary_node<sequant::utils::eval_expr> const& node;

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

};  // eval_instance

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_BTAS_HPP
