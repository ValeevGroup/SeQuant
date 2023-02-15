//
// Created by Bimal Gaudel on 7/18/21.
//

#ifndef SEQUANT_EVAL_SCF_TA_HPP
#define SEQUANT_EVAL_SCF_TA_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>
#include <SeQuant/domain/eval/eval_ta.hpp>

#include "examples/eval/calc_info.hpp"
#include "examples/eval/data_info.hpp"
#include "examples/eval/scf.hpp"
#include "examples/eval/ta/data_world_ta.hpp"

#include <chrono>

namespace sequant::eval {

template <typename Tensor_t>
class SequantEvalScfTA final : public SequantEvalScf {
 private:
  container::vector<EvalNodeTA> nodes_;
  CacheManager<Tensor_t const> cman_;
  DataWorldTA<Tensor_t> data_world_;

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

  double energy_spin_orbital() const {
    auto const& T1 = data_world_.amplitude(1);
    auto const& T2 = data_world_.amplitude(2);
    auto const& G_vvoo = g_vvoo();
    auto const& F_vo = f_vo();

    Tensor_t tau_scaled{};
    tau_scaled("a,b,i,j") = 0.25 * T2("a,b,i,j") + 0.5 * T1("a,i") * T1("b,j");

    return TA::dot(F_vo("a,i"), T1("a,i")) +
           TA::dot(G_vvoo("a,b,i,j"), tau_scaled("a,b,i,j"));
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

  void reset_cache_decaying() override { cman_.reset_decaying(); }

  void reset_cache_all() override { cman_.reset_all(); }

  double norm() const override {
    // todo use all Ts instead of only T2
    return data_world_.amplitude(2)("0,1,2,3").norm().get();
  }

  double solve() override {
    auto bk_to_labels_sorted =
        [](auto const& bk) -> container::svector<std::string> {
      auto vec = bk | ranges::views::transform([](auto const& idx) {
                   return idx.to_string();
                 }) |
                 ranges::to<container::svector<std::string>>;
      ranges::sort(vec);
      return vec;
    };

    auto tnsr_to_bk_labels_sorted =
        [&bk_to_labels_sorted](
            Tensor const& tnsr) -> container::svector<std::string> {
      auto const bra_sorted = bk_to_labels_sorted(tnsr.bra());
      auto const ket_sorted = bk_to_labels_sorted(tnsr.ket());
      return ranges::views::concat(bra_sorted, ket_sorted) |
             ranges::to<container::svector<std::string>>;
    };

    auto rs = ranges::views::repeat_n(Tensor_t{}, info_.eqn_opts.excit) |
              ranges::to_vector;
    for (auto&& [r, n] : ranges::views::zip(rs, nodes_)) {
      auto const target_indices = tnsr_to_bk_labels_sorted(n->tensor());
      r = info_.eqn_opts.spintrace
              ? eval::eval_symm(n, target_indices, data_world_, cman_)
              : eval::eval_antisymm(n, target_indices, data_world_, cman_);
    }
    data_world_.update_amplitudes(rs);
    return info_.eqn_opts.spintrace ? energy_spin_free_orbital()
                                    : energy_spin_orbital();
  }

 public:
  SequantEvalScfTA(TA::World& ta_world, CalcInfo const& calc_info)
      : SequantEvalScf{calc_info},
        cman_{{}, {}},
        data_world_{ta_world, calc_info.fock_eri, calc_info.eqn_opts.excit} {
    assert(info_.eqn_opts.excit >= 2 &&
           "At least double excitation (CCSD) is required!");

    using HRC = std::chrono::high_resolution_clock;

    auto const exprs = info_.exprs();

    auto ns = info_.nodes(exprs);

    nodes_ =
        ns |
        ranges::views::transform([](auto&& n) { return to_eval_node_ta(n); }) |
        ranges::to<decltype(nodes_)>;

    cman_ = info_.cache_manager_scf<Tensor_t const>(ns);
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_SCF_TA_HPP
