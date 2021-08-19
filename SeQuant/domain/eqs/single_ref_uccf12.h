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
 private:
  ExprPtr single_ref_approx(const ExprPtr& ex_){//densities become deltas
    auto result = ex<Constant>(0);
    for (auto&& product : ex_->as<Sum>().summands()){
      auto new_product = ex<Constant>(1);
      for (auto&& factor : product->as<Product>().factors()){
        if (factor->as<Tensor>().label() == L"\\gamma" | factor->as<Tensor>().label() == L"\\Gamma"){
          for (size_t i = 0; i < factor->as<Tensor>().bra().size(); i++){
            new_product = new_product * make_overlap(factor->as<Tensor>().bra()[i], factor->as<Tensor>().ket()[i]);
          }
        }
        else{new_product = factor * new_product;}
      }
      result = result + product;
    }
    return result;
  }
  public:
   bool sr;
   bool fock;
   unsigned int op_rank;

  uccf12(bool single_reference = true, bool fock_approx = true, unsigned int max_op_rank = 2){ sr = single_reference; fock = fock_approx; op_rank = max_op_rank;
    sequant::set_default_context(SeQuant(Vacuum::Physical, IndexSpaceMetric::Unit, BraKetSymmetry::conjugate,
                                         SPBasis::spinfree));
    mbpt::set_default_convention();
    sequant::detail::OpIdRegistrar op_id_registrar;
    TensorCanonicalizer::register_instance(std::make_shared<DefaultTensorCanonicalizer>());
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
  }

  ExprPtr do_wick(ExprPtr expr) {
    using sequant::FWickTheorem;
    FWickTheorem wick{expr};
    wick.spinfree(false).full_contractions(false);
    auto result = wick.compute();
    simplify(result);
    return result;
  }

  ExprPtr compute(bool print = false) {
    auto gg_space = IndexSpace::active_occupied;  // Geminal-generating space: active occupieds is the normal choice, all orbitals is the reference-independent (albeit expensive) choice
                                      // start transformation

    auto h = H(false);
    auto r = R12(gg_space);

    // start single comutator compute.
    auto hr_comm = do_wick(ex<Constant>(1) * (h * r - r * h));

    auto hr_comm_12 = keep_up_to_3_body_terms(hr_comm);
    auto cumu_decomp = decompositions::three_body_substitution(hr_comm_12, 2);  // apply the cumulant decompostition to the expression
    simplify(cumu_decomp);

    auto compare_h = simplification::hamiltonian_based(cumu_decomp);// simplify the expression to remove zero terms not caught by sequant core.
    auto first_comutator = compare_h + adjoint(compare_h);
    // end single comutator compute


    // end of of first commutator. now comes the second comutator one at a time, i.e. (1/2)[[f,A],A] = (1/2)[[f,R]_12, R]_12 + complex conjugates and their versions. if we get this term we can interpret the rest from properties of comutators.
    auto com_f_a =do_wick(ex<Constant>(-0.5) * (F() * adjoint(r) - adjoint(r) * F()));
    auto com_f_a_3 = keep_up_to_3_body_terms(com_f_a);
    auto com_f_a_2 = decompositions::three_body_substitution(com_f_a_3, 2);
    simplify(com_f_a_2);
    // start double comutator
    auto r_ = R12(gg_space);  // a new unique expression pointer is needed here.
    // odd, but wick needs these terms evaluated separatly
    auto comr_ = do_wick(com_f_a_2 * r_);
    auto r_com = do_wick(r_ * com_f_a_2);
    auto double_com = comr_ - r_com;
    simplify(double_com);
    auto double_com_3 = keep_up_to_3_body_terms(double_com);
    auto double_com_2 =decompositions::three_body_substitution(double_com_3, 2);
    simplify(double_com_2);
    auto fr_dag_r = simplification::fock_based(double_com_2);  // fock matrix based simplifcation which has additional screening of fock matrix from Brillouin's Theory

    auto com_f_r =do_wick(ex<Constant>(-0.5) * (F() * r - r * F()));
    auto com_f_r_3 = keep_up_to_3_body_terms(com_f_r);
    auto com_f_r_2 = decompositions::three_body_substitution(com_f_r_3, 2);
    simplify(com_f_r_2);

    // start double comutator
    // odd, but wick needs these terms evaluated separatly
    auto com_adr_ = do_wick(com_f_r_2 * adjoint(r_));
    auto r_ad_com = do_wick(adjoint(r_) * com_f_r_2);
    auto double_com_sec = com_adr_ - r_ad_com;
    simplify(double_com_sec);
    auto double_com_sec_3 = keep_up_to_3_body_terms(double_com_sec);
    auto double_com_sec_2 =decompositions::three_body_substitution(double_com_sec_3, 2);
    simplify(double_com_sec_2);
    auto frr_dag = simplification::fock_based(double_com_sec_2);

    auto total_double_com = frr_dag + fr_dag_r;
    // end double comutator.

    //add the two comutators together
    auto total = total_double_com + first_comutator;
    non_canon_simplify(total);
    if (sr) {
      total = single_ref_approx(total);
      non_canon_simplify(total);
    }
    if (print){
      std::wcout << "single reference Hbar: " << to_latex_align(total,20,3) << std::endl;
    }
    return total;
  }
};

#endif  // SEQUANT_SINGLE_REF_UCCF12_H
