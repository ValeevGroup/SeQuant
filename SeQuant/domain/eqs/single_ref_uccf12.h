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
    auto second_com_1 = do_wick((first_com * e3));
    auto second_com_2 = do_wick(e3 * first_com);
    auto second_com = second_com_1 - second_com_2;
    simplify(second_com);
    second_com = keep_up_to_3_body_terms(second_com);
    //std::wcout << to_latex_align(second_com,20,2) << std::endl;
    second_com = second_com + ex<Constant>(0.);//make a sum to avoid heavy code duplication for product and sum variants.
    second_com = simplification::overlap_with_obs(second_com);
    //std::wcout << to_latex_align(second_com,20,2) << std::endl;
    second_com = second_com + ex<Constant>(0.);
    second_com = simplification::screen_F12_and_density(second_com);
    //std::wcout << to_latex_align(second_com,20,2) << std::endl;
    second_com = simplification::tens_to_FNOps(second_com);
    second_com = decompositions::three_body_substitution(second_com,2);
    second_com = ex<Constant>(1./1) * second_com;
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

  std::pair<ExprPtr,ExprPtr> compute(std::string gg_label, bool print = false) {
    //auto gg_space = IndexSpace::active_occupied;  // Geminal-generating space: active occupieds is the normal choice, all orbitals is the reference-independent (albeit expensive) choice

    auto gg_space = IndexSpace::frozen_occupied;
    if(gg_label == "act_occ"){
      gg_space = IndexSpace::active_occupied;
    }
    else if(gg_label == "occ"){
      gg_space = IndexSpace::occupied;
    }
    else if(gg_label == "all"){
      gg_space = IndexSpace::all;
    }
    else if(gg_label == "fz"){
      gg_space = IndexSpace::frozen_occupied;
    }
    else if(gg_label == "uocc"){
      gg_space = IndexSpace::unoccupied;
    }
    else {
      throw " USUPPORTED SPACE LABEL! CHECK ABOVE FOR VALID ENTRIES";
    }

    auto h = H(false);
    auto r = R12(gg_space);
    auto r_1 = R12(gg_space);

    auto A = r - adjoint(r);
    auto H_A = do_wick(ex<Constant>(1.) * ((h * A) - (A * h)));
    auto H_A_3 = keep_up_to_3_body_terms(H_A);
    //std::wcout << "pre decomp: " << to_latex_align(single_Comm,20,2) << std::endl;
    H_A_3 = simplification::overlap_with_obs(H_A_3);
    H_A_3 = H_A_3 + ex<Constant>(0.);
    H_A_3 = simplification::screen_F12_and_density(H_A_3);
    // std::wcout << to_latex_align(H_A_3,20,2) << std::endl;
    H_A_3 = simplification::tens_to_FNOps(H_A_3);
    auto H_A_2 = decompositions::three_body_substitution(H_A_3,2);
    simplify(H_A_2);
    auto com_1 = simplification::hamiltonian_based(H_A_2);

    auto fFF = ex<Constant>(1./2) * compute_double_com(F(),r,r_1);
    non_canon_simplify(fFF);
    auto fFFt = ex<Constant>(1./2) * compute_double_com(F(),r,ex<Constant>(-1.) * adjoint(r_1));
    non_canon_simplify(fFFt);
    auto fFtFt = ex<Constant>(1./2) * compute_double_com(F(),ex<Constant>(-1.) * adjoint(r),ex<Constant>(-1.) * adjoint(r_1));
    non_canon_simplify(fFtFt);
    auto fFtF = ex<Constant>(1./2) * compute_double_com(F(),ex<Constant>(-1.) * adjoint(r),r_1);
    non_canon_simplify(fFtF);

    auto fFF_sim = simplification::fock_based(fFF);
    //std::wcout << "FF: " << to_latex_align(fFF_sim.second,20,2) << std::endl;
    auto fFFt_sim = simplification::fock_based(fFFt);
    //std::wcout << "FFt: " << to_latex_align(fFFt_sim.second,20,2) << std::endl;
    auto fFtFt_sim = simplification::fock_based(fFtFt);
    //std::wcout << "FtFt: " << to_latex_align(fFtFt_sim.second,20,2) << std::endl;
    auto fFtF_sim = simplification::fock_based(fFtF);
    //std::wcout << "FtF: " << to_latex_align(fFtF_sim.second,20,2) << std::endl;


    auto one_body = com_1.first +  (fFF_sim.first  +fFFt_sim.first + fFtFt_sim.first + fFtF_sim.first);
    auto two_body = com_1.second +  (fFF_sim.second + fFFt_sim.second + fFtFt_sim.second + fFtF_sim.second);

    //cannot use non_canon_simplify here because of B term.
    non_canon_simplify(one_body);
    non_canon_simplify(two_body);
    int term_count = 0;
    for (auto i =0; i < one_body->as<Sum>().summands().size(); i++){
      term_count +=1;
    }
    for (auto i =0; i < two_body->as<Sum>().summands().size(); i++){
      term_count +=1;
    }
    std::cout << "number of terms: " << term_count << std::endl;

    if (print){
      std::wcout << "one body terms: " << to_latex_align(one_body,20,2) << std::endl;
      std::wcout << "two body terms: " << to_latex_align(two_body,20,2) << std::endl;
    }
    return std::pair<ExprPtr, ExprPtr>{one_body, two_body};
  }
};

#endif  // SEQUANT_SINGLE_REF_UCCF12_H
