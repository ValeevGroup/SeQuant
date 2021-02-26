#ifndef SEQUANT_EVAL_TILEDARRAY_HPP
#define SEQUANT_EVAL_TILEDARRAY_HPP

#include <SeQuant/domain/optimize/optimize.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>

namespace sequant::eval {

namespace detail {
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
  auto to_str = [](auto x) { return std::to_string(x); };
  return ranges::accumulate(
      ords | ranges::views::transform(to_str) | intersperse(std::string{","}),
      std::string{}, std::plus{});
};  // ords_to_annot

template <typename DataTensor>
struct tensor_cache {
  size_t nrepeat;
  std::unique_ptr<DataTensor> ptr;
};

template <typename TensorTA, typename TensorYielder>
TensorTA _evaluate(utils::binary_node<utils::eval_expr> const& node,
                   TensorYielder&& yielder) {
  //
  auto hash_to_cache = container::map<size_t, tensor_cache<TensorTA>>{};
  node.visit_internal([&hash_to_cache](auto const& n) {
    if (auto&& found = hash_to_cache.find(n.hash());
        found != hash_to_cache.end()) {
      ++found->second.nrepeat;
    } else {
      hash_to_cache.emplace(n.hash(), tensor_cache<TensorTA>{1, nullptr});
    }
  });

  ranges::remove_if(hash_to_cache,
                    [](auto&& pair) { return pair.second.nrepeat <= 1; });

  struct {
    container::map<size_t, tensor_cache<TensorTA>> cache;
    TensorYielder&& yield;

    TensorTA operator()(utils::binary_node<utils::eval_expr> const& node) {
      return yield(node->tensor());
    }

    TensorTA operator()(utils::binary_node<utils::eval_expr> const& node,
                        TensorTA const& leval, TensorTA const& reval) {
      using ranges::views::filter;
      using ranges::views::intersperse;
      using ranges::views::transform;

      auto&& found_cache = cache.find(node->hash());
      if (found_cache != cache.end() && found_cache->second.ptr) {
        --found_cache->second.nrepeat;
        if (found_cache->second.nrepeat == 0) {
          auto to_return = std::move(*found_cache->second.ptr);
          cache.erase(found_cache);
          return to_return;
        }
        return *found_cache->second.ptr;
      }

      auto this_annot = braket_to_annot(node->tensor().braket());
      auto lannot = braket_to_annot(node.left()->tensor().braket());
      auto rannot = braket_to_annot(node.right()->tensor().braket());
      auto lscal =
          node.left()->scalar().value().real();  // assert imaginary zero: TODO
      auto rscal = node.right()->scalar().value().real();

      auto result = TensorTA{};
      if (node->op() == utils::eval_expr::eval_op::Sum) {
        result(this_annot) = lscal * leval(lannot) + rscal * reval(rannot);
      } else {
        result(this_annot) = lscal * leval(lannot) * rscal * reval(rannot);
      }

      if (found_cache != cache.end()) {
        found_cache->second.ptr = std::make_unique<TensorTA>(std::move(result));
        --found_cache->second.nrepeat;
        return *found_cache->second.ptr;
      } else {
        return result;
      }
    }
  } evaluator{std::move(hash_to_cache),
              std::forward<TensorYielder>(yielder)};  // evaluator struct
  return node.evaluate(evaluator);
}

struct perm_with_phase {
  int phase;
  container::svector<size_t> const& perm;
};

template <typename F>
void permute_ords(container::svector<size_t>& ords, F&& callback,
                  size_t beg = 0, size_t swaps = 0) {
  static_assert(std::is_invocable_v<F, perm_with_phase const&>,
                "F(perm_with_phase) not possible");

  if (beg + 1 == ords.size())
    callback(perm_with_phase{swaps % 2 == 0 ? 1 : -1, ords});

  for (auto ii = beg; ii < ords.size(); ++ii) {
    std::swap(ords[beg], ords[ii]);
    permute_ords(ords, std::forward<F>(callback), beg + 1,
                 ii == beg ? swaps : swaps + 1);
    std::swap(ords[ii], ords[beg]);
  }
}

}  // namespace detail

template <typename TensorTA>
TensorTA symmetrize(TensorTA const& tensor) {
  using ranges::views::concat;
  using ranges::views::iota;
  using ranges::views::transform;

  auto result = TensorTA{tensor.world(), tensor.trange()};
  result.fill(0);

  // returns a lambda that transforms a given range by adding
  // n to the elements
  auto n_adder = [](size_t n) { return [n](auto x) { return x + n; }; };

  auto const rank = tensor.range().rank();
  auto const lannot = detail::ords_to_annot(iota(size_t{0}, size_t{rank}));
  auto perm_vec =
      iota(size_t{0}, rank / 2) | ranges::to<container::svector<size_t>>;
  do {
    auto const rannot = detail::ords_to_annot(
        concat(perm_vec, perm_vec | transform(n_adder(rank / 2))));

    result(lannot) += tensor(rannot);

  } while (std::next_permutation(perm_vec.begin(), perm_vec.end()));

  return result;
}

template <typename TensorTA>
TensorTA antisymmetrize(TensorTA const& tensor) {
  auto const rank = tensor.range().rank();
  auto const ords = ranges::views::iota(size_t{0}, size_t{rank});

  auto result = TensorTA{tensor.world(), tensor.trange()};
  result.fill(0);
  auto asym_impl = [lannot = detail::ords_to_annot(ords), &tensor, &result](
                       int phase, std::string_view rannot) -> void {
    result(lannot) += phase * tensor(rannot.data());
    // std::cout << "lannot = " << lannot << "  rannot = " << rannot << "\n";
  };  // asym_impl

  auto pred = [&asym_impl, rank](auto const& bra_pwp) {
    auto pred_inner = [&asym_impl, &bra_pwp](auto const& ket_pwp) {
      auto const phase = bra_pwp.phase * ket_pwp.phase;
      auto const annot = detail::ords_to_annot(
          ranges::views::concat(bra_pwp.perm, ket_pwp.perm));
      asym_impl(phase, annot);
    };

    auto ket_ords = ranges::views::iota(rank / 2, rank) |
                    ranges::to<container::svector<size_t>>;

    detail::permute_ords(ket_ords, pred_inner);
  };

  auto bra_ords = ords | ranges::views::take(rank / 2) |
                  ranges::to<container::svector<size_t>>;
  detail::permute_ords(bra_ords, pred);

  return result;
}

template <typename TensorTA, typename TensorYielder>
TensorTA evaluate(ExprPtr const& expr, bool opt, TensorYielder&& yielder) {
  static_assert(
      std::is_invocable_r_v<TensorTA, TensorYielder, sequant::Tensor const&>,
      "TensorYielder(sequant::Tensor const&) should return TensorTA");

  using ranges::views::filter;
  using ranges::views::transform;

  if (!(expr->is<Sum>() or expr->is<Product>()))
    throw std::runtime_error(
        "evaluate called on non-Sum and non-Product type expression");

  expr->visit(
      [](ExprPtr const& xpr) {
        if (auto const lbl = xpr->as<Tensor>().label();
            lbl == L"A" || lbl == L"S")
          throw std::runtime_error(
              "A and S tensors should be factored out before evaluation");
      },
      true);  // visit only atoms

  const size_t nocc = 2, nvirt = 20;

  if (expr->is<Product>()) {
    auto eval_node =
        opt ? std::move(*optimize::single_term_opt(expr->as<Product>(), nocc,
                                                    nvirt, {}, true)
                             .optimal_seqs.begin())
            : utils::binarize_expr(expr);

    auto res = TensorTA{};
    auto annot = detail::braket_to_annot(eval_node->tensor().const_braket());
    res(annot) = eval_node->scalar().value().real() *
                 detail::_evaluate<TensorTA>(
                     eval_node, std::forward<TensorYielder>(yielder))(annot);
    return res;
  } else {
    auto summands = *expr | transform([opt](ExprPtr p) {
      return (p->is<Tensor>() or !opt)
                 ? utils::binarize_expr(p)
                 : std::move(*optimize::single_term_opt(p->as<Product>(), nocc,
                                                         nvirt, {}, true)
                                  .optimal_seqs.begin());
    }) | transform([](auto const& node) {
      return utils::debinarize_eval_expr(node);
    });

    auto eval_node = utils::binarize_expr(
        ex<Sum>(ranges::begin(summands), ranges::end(summands)));

    return detail::_evaluate<TensorTA>(eval_node,
                                       std::forward<TensorYielder>(yielder));
  }
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_TILEDARRAY_HPP
