//
// Created by Bimal Gaudel on 7/18/21.
//

#ifndef SEQUANT_EVAL_SCF_TA_HPP
#define SEQUANT_EVAL_SCF_TA_HPP

#include "examples/eval/calc_info.hpp"
#include "examples/eval/data_info.hpp"
#include "examples/eval/scf.hpp"
#include "examples/eval/ta/eval_ta.hpp"
#include "examples/eval/ta/data_world_ta.hpp"

#include <SeQuant/domain/eval/cache_manager.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/parse_expr.hpp>

#include <chrono>

namespace sequant::eval::ta {

template <typename Tensor_t>
class SequantEvalScfTA final: public SequantEvalScf {
 private:
  container::vector<EvalNode> nodes_;
  CacheManager<Tensor_t const> cman_;
  DataWorldTA<Tensor_t> data_world_;

 private:

  Tensor_t const& f_vo() const {
    static Tensor_t tnsr = data_world_(parse_expr(L"f{a1;i1}", Symmetry::nonsymm)->as<Tensor>());
    return tnsr;
  }

  Tensor_t const& g_vvoo() const {
   static Tensor_t tnsr = data_world_(parse_expr(L"g{a1,a2;i1,i2}", Symmetry::nonsymm)->as<Tensor>());
   return tnsr;
  }

  double energy_spin_orbital() const {
    auto const& T1 = data_world_.amplitude(1);
    auto const& T2 = data_world_.amplitude(2);
    auto const& G_vvoo = g_vvoo();
    auto const& F_vo   = f_vo();

    Tensor_t tau_scaled{};
    tau_scaled("a,b,i,j") =
          0.25  * T2("a,b,i,j") + 0.5 * T1("a,i") * T1("b,j");

    return TA::dot(F_vo("a,i"), T1("a,i"))
           + TA::dot(G_vvoo("a,b,i,j"), tau_scaled("a,b,i,j"));
  }

  double energy_spin_free_orbital() const {
    auto const& T1 = data_world_.amplitude(1);
    auto const& T2 = data_world_.amplitude(2);
    auto const& G_vvoo = g_vvoo();
    auto const& F_vo   = f_vo();

    Tensor_t tau{};
    tau("a,b,i,j") = T2("a,b,i,j") + T1("a,i") * T1("b,j");

    return 2.0 * TA::dot(F_vo("a,i"), T1("a,i"))
             + TA::dot(2 * G_vvoo("a,b,i,j") - G_vvoo("a,b,j,i"), tau("a,b,i,j"));
  }

  void reset_cache_decaying() override {
    cman_.reset_decaying();
  }

  void reset_cache_all() override {
    cman_.reset_all();
  }

  double norm() const override {
    // todo use all Ts instead of only T2
    return data_world_.amplitude(2)("0,1,2,3").norm().get();
  }

  double solve() override {
   auto rs = ranges::views::repeat_n(Tensor_t{}, info_.eqn_opts.excit)
             | ranges::to_vector;
   for (auto&& [r, n] : ranges::views::zip(rs, nodes_))
     r = info_.eqn_opts.spintrace ? eval::ta::eval_symm(n, data_world_, cman_)
                                : eval::ta::eval_antisymm(n, data_world_, cman_);
   data_world_.update_amplitudes(rs);
   return info_.eqn_opts.spintrace ? energy_spin_free_orbital()
                                 : energy_spin_orbital();
  }

 public:
  SequantEvalScfTA(CalcInfo const& calc_info, TA::World& ta_world)
      : SequantEvalScf{calc_info}, cman_{{}, {}},
        data_world_{ta_world, calc_info.eqn_opts.excit, calc_info.fock_eri} {

    assert(info_.eqn_opts.excit >= 2
           && "At least double excitation (CCSD) is required!");

    using HRC = std::chrono::high_resolution_clock;

    auto const exprs = info_.exprs();

    nodes_ = info_.nodes(exprs);

    cman_  = info_.cache_manager<Tensor_t const>(nodes_);
  }

};

} // namespace

#endif  // SEQUANT_EVAL_SCF_TA_HPP
