//
// Created by Bimal Gaudel on 7/18/21.
//

#ifndef SEQUANT_EVAL_SCF_TA_HPP
#define SEQUANT_EVAL_SCF_TA_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>
#include <SeQuant/domain/eval/eval.hpp>

#include <calc_info.hpp>
#include <data_info.hpp>
#include <scf.hpp>
#include <ta/data_world_ta.hpp>

#include <chrono>

namespace sequant::eval {

template <typename Tensor_t>
class SequantEvalScfTA final : public SequantEvalScf {
 public:
  using ExprT = EvalExprTA;
  using EvalNodeTA = EvalNode<ExprT>;

 private:
  container::vector<container::vector<EvalNodeTA>> nodes_;
  CacheManager cman_;
  DataWorldTA<Tensor_t> data_world_;

  Tensor_t const& f_vo() const {
    static Tensor_t tnsr =
        data_world_(parse_expr(L"f{a1;i1}", Symmetry::Nonsymm)->as<Tensor>());
    return tnsr;
  }

  Tensor_t const& g_vvoo() const {
    static Tensor_t tnsr = data_world_(
        parse_expr(L"g{a1,a2;i1,i2}", Symmetry::Nonsymm)->as<Tensor>());
    return tnsr;
  }

  double energy_spin_orbital() {
    static const std::wstring_view energy_expr =
        L"f{i1;a1} * t{a1;i1} + g{i1,i2;a1,a2} * "
        L"(1/4 * t{a1,a2;i1,i2} + 1/2 t{a1;i1} * t{a2;i2})";
    static auto const node =
        binarize<EvalExprTA>(parse_expr(energy_expr, Symmetry::Antisymm));

    return evaluate(node, data_world_)->template get<double>();
  }

  double energy_spin_free_orbital() const {
    auto const& T1 = data_world_.amplitude(1);
    auto const& T2 = data_world_.amplitude(2);
    auto const& G_vvoo = g_vvoo();
    auto const& F_vo = f_vo();

    Tensor_t tau{};
    tau("a,b,i,j") = T2("a,b,i,j") + T1("a,i") * T1("b,j");

    return 2.0 * TA::dot(F_vo("a,i"), T1("a,i")) +
           TA::dot(2 * G_vvoo("a,b,i,j") - G_vvoo("a,b,j,i"), tau("a,b,i,j"));
  }

  void reset_cache() override { cman_.reset(); }

  double norm() const override {
    // todo use all Ts instead of only T2
    return data_world_.amplitude(2)("0,1,2,3").norm().get();
  }

  double solve() override {
    using ranges::views::concat;
    using ranges::views::intersperse;
    using ranges::views::join;
    using ranges::views::repeat_n;
    using ranges::views::transform;
    using ranges::views::zip;

    using Evaluator =
        ResultPtr (*)(std::vector<EvalNodeTA> const&, std::string const&,
                      DataWorldTA<Tensor_t> const&, CacheManager&);

    constexpr auto funcs =
        std::array{std::array<Evaluator, 2>{evaluate_antisymm<Trace::Off>,
                                            evaluate_antisymm<Trace::On>},
                   std::array<Evaluator, 2>{evaluate_symm<Trace::Off>,
                                            evaluate_symm<Trace::On>}};

    auto rs = repeat_n(Tensor_t{}, info_.eqn_opts.excit) | ranges::to_vector;
    auto st = info_.eqn_opts.spintrace;
    auto log = Logger::instance().eval.level > 0;

    for (auto&& [r, n] : zip(rs, nodes_)) {
      SEQUANT_ASSERT(!n.empty());
      auto trank = ranges::front(n)->as_tensor().bra_rank();
      SEQUANT_ASSERT(trank == 1 || trank == 2);
      // update_amplitudes assumes that the residuals are in specific layout
      std::string target_indices = trank == 1 ? "a_1,i_1" : "a_1,a_2,i_1,i_2";
      r = funcs[st][log](n, target_indices, data_world_, cman_)
              ->template get<Tensor_t>();
    }

    auto E = info_.eqn_opts.spintrace ? energy_spin_free_orbital()
                                      : energy_spin_orbital();
    data_world_.update_amplitudes(rs);
    return E;
  }

 public:
  SequantEvalScfTA(TA::World& ta_world, CalcInfo const& calc_info)
      : SequantEvalScf{calc_info},
        cman_{{}},
        data_world_{ta_world, calc_info.fock_eri, calc_info.eqn_opts.excit} {
    assert(info_.eqn_opts.excit >= 2 &&
           "At least double excitation (CCSD) is required!");

    auto const exprs = info_.exprs();

    nodes_ = info_.nodes<ExprT>(exprs);

    cman_ = info_.cache_manager_scf(nodes_);
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_SCF_TA_HPP
