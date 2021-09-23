//
// Created by Conner Masteran on 8/16/21.
//

#ifndef SEQUANT_SINGLE_REF_UCCF12_H
#define SEQUANT_SINGLE_REF_UCCF12_H
#include "../transcorrelated/three_body_decomp.hpp"
#include "../transcorrelated/simplifications.h"
#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/sr/sr.hpp>

#include <clocale>

using namespace sequant;
using namespace sequant::mbpt::sr::so;

class uccf12{

  public:
   bool sr;
   bool fock;
   unsigned int op_rank;
  //TODO implement logic for non-default variables. should also include logic for spin-orbital expressions.
  uccf12(bool single_reference = true, bool fock_approx = true, unsigned int max_op_rank = 2){ sr = single_reference; fock = fock_approx; op_rank = max_op_rank;
    sequant::set_default_context(SeQuant(Vacuum::Physical, IndexSpaceMetric::Unit, BraKetSymmetry::conjugate,
                                         SPBasis::spinfree));
//    mbpt::set_default_convention();
    sequant::detail::OpIdRegistrar op_id_registrar;
    TensorCanonicalizer::register_instance(std::make_shared<DefaultTensorCanonicalizer>());
  }
  //[[e1,e2],e3]_12
  ExprPtr compute_double_com(ExprPtr e1, ExprPtr e2, ExprPtr e3){
    auto first_com = do_wick((e1 * e2) - (e2 * e1));
    auto first_com_clone = first_com->clone();
    auto second_com_1 = do_wick((first_com_clone * e3));
    auto second_com_2 = do_wick(e3 * first_com);
    auto second_com = second_com_1 - second_com_2;
    simplify(second_com);
    second_com = keep_up_to_3_body_terms(second_com);
    second_com = decompositions::three_body_substitution(second_com,2);
    simplify(second_com);
    return second_com;
  }

  ExprPtr keep_up_to_3_body_terms(const ExprPtr& input) {
    if (input->is<Sum>()) {auto filtered_summands = input->as<Sum>().summands() |
          ranges::views::remove_if([](const ExprPtr& ptr) {assert(ptr->is<Product>());
            bool keep = false;
            bool found_operator = false;
            for (auto&& factor : ptr->as<Product>().factors()) {
              if (factor->is<FNOperator>()) {
                assert(!found_operator);
                found_operator = true;
                const auto rank = factor->as<FNOperator>().rank();
                keep = (rank <= 3);
              }
            }
            return !keep;
          });
      auto result = ex<Sum>(ranges::begin(filtered_summands),
                            ranges::end(filtered_summands));
      return result;
    }
    else if (input->is<Product>()){
      for(auto&& factor : input->as<Product>().factors()){
        if(factor->is<FNOperator>()){
          if(factor->as<FNOperator>().rank() <= 3){
            return input;
          }
          else{
            return ex<Constant>(0);
          }
        }
      }
    }
  }

  ExprPtr do_wick(ExprPtr expr) {
    using sequant::FWickTheorem;
    FWickTheorem wick{expr};
    wick.spinfree(false).full_contractions(false);
    auto result = wick.compute();
    simplify(result);
    return result;
  }

  std::pair<ExprPtr,ExprPtr> compute(bool print = false) {
    auto gg_space = IndexSpace::occupied;  // Geminal-generating space: active occupieds is the normal choice, all orbitals is the reference-independent (albeit expensive) choice
                                      // start transformation

    auto h = H(false);
    auto r = R12(gg_space);
    auto r_1 = R12(gg_space);

    auto A = r - adjoint(r_1);
    auto H_A = do_wick(ex<Constant>(1.) *((h * A) - (A * h)));
    auto H_A_3 = keep_up_to_3_body_terms(H_A);
    auto H_A_2 = decompositions::three_body_substitution(H_A_3,2);
    simplify(H_A_2);
    auto com_1 = simplification::hamiltonian_based(H_A_2);

    auto fFF = compute_double_com(F(),r,r_1);
    auto fFFt = compute_double_com(F(),r,ex<Constant>(-1.) * adjoint(r_1));
    auto fFtFt = compute_double_com(F(),ex<Constant>(-1.) * adjoint(r),ex<Constant>(-1.) * adjoint(r_1));
    auto fFtF = compute_double_com(F(),ex<Constant>(-1.) * adjoint(r),r_1);

    auto fFF_sim = simplification::fock_based(fFF);
    //std::wcout << "FF: " << to_latex_align(fFF_sim.second,20,2) << std::endl;
    auto fFFt_sim = simplification::fock_based(fFFt);
    //std::wcout << "FFt: " << to_latex_align(fFFt_sim.second,20,2) << std::endl;
    auto fFtFt_sim = simplification::fock_based(fFtFt);
    //std::wcout << "FtFt: " << to_latex_align(fFtFt_sim.second,20,2) << std::endl;
    auto fFtF_sim = simplification::fock_based(fFtF);
    //std::wcout << "FtF: " << to_latex_align(fFtF_sim.second,20,2) << std::endl;


    auto one_body = com_1.first + ex<Constant>(.5) * (fFF_sim.first  +fFFt_sim.first + fFtFt_sim.first + fFtF_sim.first);
    auto two_body = com_1.second + ex<Constant>(.5) * (fFF_sim.second + fFFt_sim.second + fFtFt_sim.second + fFtF_sim.second);

    non_canon_simplify(one_body);
    non_canon_simplify(two_body);

    if (print){
      std::wcout << "one body terms: " << to_latex_align(one_body,20,2) << std::endl;
      std::wcout << "two body terms: " << to_latex_align(two_body,20,2) << std::endl;
    }
    return std::pair<ExprPtr, ExprPtr>{one_body, two_body};
  }
};

#endif  // SEQUANT_SINGLE_REF_UCCF12_H
