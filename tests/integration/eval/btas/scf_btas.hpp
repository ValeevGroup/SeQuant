//
// Created by Bimal Gaudel on 7/31/21.
//

#ifndef SEQUANT_EVAL_SCF_BTAS_HPP
#define SEQUANT_EVAL_SCF_BTAS_HPP

#include <btas/data_world_btas.hpp>
#include <calc_info.hpp>
#include <data_info.hpp>
#include <scf.hpp>

#include <SeQuant/domain/eval/cache_manager.hpp>
#include <SeQuant/domain/eval/eval.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/parse.hpp>

#include <btas/btas.h>
#include <memory>

namespace sequant::eval::btas {

template <typename Tensor_t>
class SequantEvalScfBTAS final : public SequantEvalScf {
 public:
  using ExprT = EvalExprBTAS;
  using EvalNodeBTAS = EvalNode<ExprT>;

 private:
  container::vector<container::vector<EvalNodeBTAS>> nodes_;
  CacheManager cman_;
  DataWorldBTAS<Tensor_t> data_world_;

  Tensor_t const& f_vo() const {
    static Tensor_t tnsr =
        data_world_(parse_expr(L"f{a1;i1}", Symmetry::nonsymm)->as<Tensor>());
    return tnsr;
  }

  Tensor_t const& g_vvoo() const {
    static Tensor_t tnsr = data_world_(
        parse_expr(L"g{a1,a2;i1,i2}", Symmetry::nonsymm)->as<Tensor>());
    return tnsr;
  }

  double energy_spin_orbital() {
    static const std::wstring_view energy_expr =
        L"f{i1;a1} * t{a1;i1} + g{i1,i2;a1,a2} * "
        L"(1/4 * t{a1,a2;i1,i2} + 1/2 t{a1;i1} * t{a2;i2})";
    static auto const node =
        binarize<EvalExprBTAS>(parse_expr(energy_expr, Symmetry::antisymm));

    return evaluate(node, data_world_)->template get<double>();
  }

  double energy_spin_free_orbital() const {
    auto const& T1 = data_world_.amplitude(1);
    auto const& T2 = data_world_.amplitude(2);
    auto const& G_vvoo = g_vvoo();
    auto const& F_vo = f_vo();

    Tensor_t temp;
    ::btas::contract(1.0, T1, {'a', 'i'}, T1, {'b', 'j'}, 0.0, temp,
                     {'a', 'b', 'i', 'j'});

    return 2.0 * (::btas::dot(F_vo, T1) + ::btas::dot(G_vvoo, T2) +
                  ::btas::dot(G_vvoo, temp)) -
           ::btas::dot(G_vvoo, temp);
  }

  void reset_cache() override { cman_.reset(); }

  double norm() const override {
    // todo use all Ts instead of only T2
    auto const& T2 = data_world_.amplitude(2);
    return std::sqrt(::btas::dot(T2, T2));
  }

  double solve() override {
    using ranges::views::concat;
    using ranges::views::repeat_n;
    using ranges::views::transform;
    using ranges::views::zip;

    using Evaluator = ResultPtr (*)(
        std::vector<EvalNodeBTAS> const&, EvalExprBTAS::annot_t const&,
        DataWorldBTAS<Tensor_t> const&, CacheManager&);

    constexpr auto funcs =
        std::array{std::array<Evaluator, 2>{evaluate_antisymm<Trace::Off>,
                                            evaluate_antisymm<Trace::On>},
                   std::array<Evaluator, 2>{evaluate_symm<Trace::Off>,
                                            evaluate_symm<Trace::On>}};

    auto sorted_annot = [](Tensor const& tnsr) {
      auto b = tnsr.bra() | ranges::to_vector;
      auto k = tnsr.ket() | ranges::to_vector;
      ranges::sort(b, Index::LabelCompare{});
      ranges::sort(k, Index::LabelCompare{});
      return EvalExprBTAS::index_hash(concat(b, k)) |
             ranges::to<EvalExprBTAS::annot_t>;
    };

    auto rs = repeat_n(Tensor_t{}, info_.eqn_opts.excit) | ranges::to_vector;
    auto st = info_.eqn_opts.spintrace;
    auto log = Logger::instance().eval.level > 0;

    for (auto&& [r, n] : zip(rs, nodes_)) {
      auto const& target_indices = ranges::front(n)->annot();
      r = funcs[st][log](n, target_indices, data_world_, cman_)
              ->template get<Tensor_t>();
    }

    data_world_.update_amplitudes(rs);
    return info_.eqn_opts.spintrace ? energy_spin_free_orbital()
                                    : energy_spin_orbital();
  }

 public:
  SequantEvalScfBTAS(CalcInfo const& calc_info)
      : SequantEvalScf{calc_info},
        cman_{{}},
        data_world_{calc_info.fock_eri, calc_info.eqn_opts.excit} {
    assert(info_.eqn_opts.excit >= 2 &&
           "At least double excitation (CCSD) is required!");
    // todo time it
    auto const exprs = info_.exprs();

    nodes_ = info_.nodes<ExprT>(exprs);

    cman_ = info_.cache_manager_scf(nodes_);
  }
};

}  // namespace sequant::eval::btas

#endif  // SEQUANT_EVAL_SCF_BTAS_HPP
