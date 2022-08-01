//
// Created by Conner Masteran on 8/16/21.
//

#ifndef SEQUANT_DOMAIN_SINGLE_REF_UCCF12_H
#define SEQUANT_DOMAIN_SINGLE_REF_UCCF12_H

#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/sr/sr.hpp>
#include "../transcorrelated/simplifications.h"
#include "../transcorrelated/three_body_decomp.hpp"

#include <clocale>

using namespace sequant;
using namespace sequant::mbpt::sr::so;

class uccf12 {
 public:
  bool sr;
  bool fock;
  unsigned int op_rank;
  //no default constructor
  IndexSpace::TypeAttr gg_space = IndexSpace::all;
  IndexSpace::TypeAttr singles_gg_space = IndexSpace::all;
  int ansatz_;
  int comutator_order_;
  bool print_;
  bool singles_;
  bool doubles_;

  //default constructor
  uccf12(){}

  // TODO implement logic for non-default variables. should also include logic
  // for spin-orbital expressions.
  uccf12(std::string gg_label, int ansatz = 2,
         bool print = false, bool singles = false, std::string singles_gg = "all",
         bool doubles = true, int comutator_order = 2) {
    ansatz_ = ansatz;
    print_ = print;
    singles_ = singles;
    doubles_ = doubles;
    comutator_order_ = comutator_order;



    sequant::set_default_context(
        SeQuant(Vacuum::Physical, IndexSpaceMetric::Unit,
                BraKetSymmetry::conjugate, SPBasis::spinfree));
    // mbpt::set_default_convention();
    std::setlocale(LC_ALL, "en_US.UTF-8");
    std::wcout.precision(std::numeric_limits<double>::max_digits10);
    std::wcerr.precision(std::numeric_limits<double>::max_digits10);
    std::wcout.sync_with_stdio(false);
    std::wcerr.sync_with_stdio(false);
    std::wcout.imbue(std::locale("en_US.UTF-8"));
    std::wcerr.imbue(std::locale("en_US.UTF-8"));
    std::wcout.sync_with_stdio(true);
    std::wcerr.sync_with_stdio(true);
    sequant::detail::OpIdRegistrar op_id_registrar;
    TensorCanonicalizer::register_instance(
        std::make_shared<DefaultTensorCanonicalizer>());

    // auto gg_space = IndexSpace::active_occupied;  // Geminal-generating
    // space: active occupieds is the normal choice, all orbitals is the
    // reference-independent (albeit expensive) choice
    assert(singles_ == true || doubles_ == true);
    //doubles space options
    if (gg_label == "act_occ") {
      gg_space = IndexSpace::active_occupied;
    } else if (gg_label == "occ") {
      gg_space = IndexSpace::occupied;
    } else if (gg_label == "all") {
      gg_space = IndexSpace::all;
    } else if (gg_label == "fz") {
      gg_space = IndexSpace::frozen_occupied;
    } else if (gg_label == "uocc") {
      gg_space = IndexSpace::unoccupied;
    }
    else if (gg_label == "act_obs") {
      gg_space = IndexSpace::all_active;
    } else {
      throw std::runtime_error(
          "uccf12::compute(gg_label) unsupported space label");
    }
    //singles space options
    if (singles_gg == "all"){
      singles_gg_space = IndexSpace::all;
    }
    else if (singles_gg == "occ") {
      singles_gg_space = IndexSpace::occupied;
    }
    else if (singles_gg == "all_active"){
      singles_gg_space = IndexSpace::all_active;
    }
    else if(singles_gg == "occ_active"){
      singles_gg_space = IndexSpace::active_occupied;
    }
    else{throw "singles index space not supported!";}
  }
  //[[e1,e2],e3]_12
  ExprPtr compute_double_com(ExprPtr e1, ExprPtr e2, ExprPtr e3,
                             int ansatz = 2) {
    auto first_com = do_wick((e1 * e2) - (e2 * e1));
    non_canon_simplify(first_com);
    auto second_com_1 = first_com * e3;
    non_canon_simplify(second_com_1);
    second_com_1 = do_wick(second_com_1);
    auto second_com_2 = e3 * first_com;
    non_canon_simplify(second_com_2);
    second_com_2 = do_wick(second_com_2);
    auto second_com = second_com_1 - second_com_2;
    non_canon_simplify(second_com);
    if (ansatz == 2) {
      second_com = keep_up_to_3_body_terms(second_com);
      simplify(second_com);
      second_com = simplification::detail::tens_to_FNOps(second_com);
      second_com = decompositions::three_body_substitution(second_com, 2);
      simplify(second_com);
      second_com =
          second_com +
          ex<Constant>(0.);  // make a sum to avoid heavy code duplication for
                             // product and sum variants.
      second_com = simplification::overlap_with_obs(second_com);
      second_com = second_com + ex<Constant>(0.);
      second_com = simplification::screen_F12_proj(second_com, 2);
      //std::wcout << to_latex_align(second_com,20,3) << std::endl;
      return second_com;
    }
    if (ansatz == 1) {
      second_com = keep_up_to_2_body_terms(second_com);
      // std::wcout << to_latex_align(second_com,20,2) << std::endl;
      second_com =
          second_com +
          ex<Constant>(0.);  // make a sum to avoid heavy code duplication for
                             // product and sum variants.
      second_com = simplification::overlap_with_obs(second_com);
      // std::wcout << to_latex_align(second_com,20,2) << std::endl;
      second_com = second_com + ex<Constant>(0.);
      second_com = simplification::screen_F12_proj(second_com, 1);
      // std::wcout << to_latex_align(second_com,20,2) << std::endl;
      second_com = simplification::detail::tens_to_FNOps(second_com);
      simplify(second_com);
      return second_com;
    }
  }

  ExprPtr keep_up_to_3_body_terms(const ExprPtr& input) {
    if (input->is<Sum>()) {
      auto filtered_summands =
          input->as<Sum>().summands() |
          ranges::views::remove_if([](const ExprPtr& ptr) {
            assert(ptr->is<Product>());
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
    } else if (input->is<Product>()) {
      for (auto&& factor : input->as<Product>().factors()) {
        if (factor->is<FNOperator>()) {
          if (factor->as<FNOperator>().rank() <= 3) {
            return input;
          } else {
            return ex<Constant>(0);
          }
        }
      }
    } else {
      return input;
    }
  }
  ExprPtr keep_up_to_2_body_terms(const ExprPtr& input) {
    if (input->is<Sum>()) {
      auto filtered_summands =
          input->as<Sum>().summands() |
          ranges::views::remove_if([](const ExprPtr& ptr) {
            assert(ptr->is<Product>());
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
    } else if (input->is<Product>()) {
      for (auto&& factor : input->as<Product>().factors()) {
        if (factor->is<FNOperator>()) {
          if (factor->as<FNOperator>().rank() <= 2) {
            return input;
          } else {
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

  ExprPtr do_com (ExprPtr t1, ExprPtr t2, int ansatz = 2){
    // wick's theorm
    auto expression = (t1 * t2) - (t2 * t1);
    non_canon_simplify(expression);
    auto first_com = do_wick(expression);
    // simplification
    non_canon_simplify(first_com);
    return first_com;
  }

  ExprPtr post_wick_sim(ExprPtr first_com, int ansatz = 2){
    if (ansatz == 2) {
      first_com = keep_up_to_3_body_terms(first_com);
      first_com = first_com + ex<Constant>(0.);// make a sum to avoid heavy code duplication for product and sum variants.
      first_com = simplification::overlap_with_obs(first_com);
      first_com = first_com + ex<Constant>(0.);
      simplify(first_com);
      first_com = simplification::screen_F12_proj(first_com, 2);
      non_canon_simplify(first_com);
      first_com = simplification::detail::tens_to_FNOps(first_com);
      non_canon_simplify(first_com);
      first_com = decompositions::three_body_substitution(first_com, 2);
      non_canon_simplify(first_com);
      auto temp = ex<Constant>(0.0);
      for (auto&& product : first_com->as<Sum>().summands()){
        auto temp_temp = product + ex<Constant>(0.0);
        simplify(temp_temp);
        temp = temp_temp + temp;
      }
      simplify(temp);
      return temp;
    }
    if (ansatz == 1) {
      first_com = keep_up_to_2_body_terms(first_com);
      first_com =
          first_com +
          ex<Constant>(0.);
      first_com = simplification::overlap_with_obs(first_com);
      first_com = first_com + ex<Constant>(0.);
      first_com = simplification::screen_F12_proj(first_com, 1);
      first_com = simplification::detail::tens_to_FNOps(first_com);
      simplify(first_com);
      return first_com;
    }
  }

  // produces a uniquely indexed version of the given expression.
  // assumes same number of upper and lower indices for operators and tensors
  // do not simplify(expr) after use! this will cannonicalize labeling, undoing
  // this work.
  ExprPtr relable(ExprPtr expr) {
    if (expr->is<Sum>()) {
      auto new_sum = ex<Constant>(0.0);
      for (auto&& product : expr->as<Sum>().summands()) {
        auto new_product = relable(product);
        new_sum = new_sum + new_product;
      }
      return new_sum;
    }

    // product does not benefit from recursion
    //  must reproduce same connectivity to produce identical expressions.
    else if (expr->is<Product>()) {
      std::vector<Index> changed_indices;   // list of original indices
      std::vector<Index> original_indices;  // list of new indices
      auto new_product = ex<Constant>(expr->as<Product>().scalar());
      for (auto&& factor : expr->as<Product>().factors()) {
        std::pair<std::vector<Index>, std::vector<Index>> new_up_low;
        if (factor->is<Tensor>()) {
          for (int i = 0; i < factor->as<Tensor>().bra().size(); i++) {
            auto in_where_bra = simplification::detail::in_list(
                factor->as<Tensor>().bra()[i], original_indices);
            if (in_where_bra.first) {
              new_up_low.first.push_back(changed_indices[in_where_bra.second]);
            } else {
              original_indices.push_back(factor->as<Tensor>().bra()[i]);
              changed_indices.push_back(
                  Index::make_tmp_index(IndexSpace::instance(
                      factor->as<Tensor>().bra()[i].space().attr())));
              new_up_low.first.push_back(
                  changed_indices[changed_indices.size() - 1]);
            }
            auto in_where_ket = simplification::detail::in_list(
                factor->as<Tensor>().ket()[i], original_indices);
            if (in_where_ket.first) {
              new_up_low.second.push_back(changed_indices[in_where_ket.second]);
            } else {
              original_indices.push_back(factor->as<Tensor>().ket()[i]);
              changed_indices.push_back(
                  Index::make_tmp_index(IndexSpace::instance(
                      factor->as<Tensor>().ket()[i].space().attr())));
              new_up_low.second.push_back(
                  changed_indices[changed_indices.size() - 1]);
            }
          }
          auto new_factor = ex<Tensor>(factor->as<Tensor>().label(),
                                       new_up_low.first, new_up_low.second);
          new_product = new_product * new_factor;
        } else if (factor->is<FNOperator>()) {
          for (int i = 0; i < factor->as<FNOperator>().nannihilators(); i++) {
            auto in_where_ann = simplification::detail::in_list(
                factor->as<FNOperator>().annihilators()[i].index(),
                original_indices);
            if (in_where_ann.first) {
              new_up_low.first.push_back(changed_indices[in_where_ann.second]);
            } else {
              original_indices.push_back(
                  factor->as<FNOperator>().annihilators()[i].index());
              changed_indices.push_back(Index::make_tmp_index(
                  IndexSpace::instance(factor->as<FNOperator>()
                                           .annihilators()[i]
                                           .index()
                                           .space()
                                           .attr())));
              new_up_low.first.push_back(
                  changed_indices[changed_indices.size() - 1]);
            }
            auto in_where_cre = simplification::detail::in_list(
                factor->as<FNOperator>().creators()[i].index(),
                original_indices);
            if (in_where_cre.first) {
              new_up_low.second.push_back(changed_indices[in_where_cre.second]);
            } else {
              original_indices.push_back(
                  factor->as<FNOperator>().creators()[i].index());
              changed_indices.push_back(Index::make_tmp_index(
                  IndexSpace::instance(factor->as<FNOperator>()
                                           .creators()[i]
                                           .index()
                                           .space()
                                           .attr())));
              new_up_low.second.push_back(
                  changed_indices[changed_indices.size() - 1]);
            }
          }
          auto new_factor = ex<FNOperator>(new_up_low.second, new_up_low.first);
          new_product = new_product * new_factor;
        } else {
          throw "unsupported factor type";
        }
      }
      return new_product;
    } else if (expr->is<Tensor>()) {
      std::pair<std::vector<Index>, std::vector<Index>> new_bra_ket;
      for (int i = 0; i < expr->as<Tensor>().bra().size(); i++) {
        new_bra_ket.first.push_back(Index::make_tmp_index(
            IndexSpace::instance(expr->as<Tensor>().bra()[i].space().attr())));
        new_bra_ket.second.push_back(Index::make_tmp_index(
            IndexSpace::instance(expr->as<Tensor>().ket()[i].space().attr())));
      }
      return ex<Tensor>(expr->as<Tensor>().label(), new_bra_ket.first,
                        new_bra_ket.second);
    } else if (expr->is<FNOperator>()) {
      std::pair<std::vector<Index>, std::vector<Index>> new_ann_cre;
      for (int i = 0; i < expr->as<FNOperator>().nannihilators(); i++) {
        new_ann_cre.first.push_back(Index::make_tmp_index(IndexSpace::instance(
            expr->as<FNOperator>().annihilators()[i].index().space().attr())));
        new_ann_cre.second.push_back(Index::make_tmp_index(IndexSpace::instance(
            expr->as<FNOperator>().creators()[i].index().space().attr())));
      }
      return ex<FNOperator>(new_ann_cre.first, new_ann_cre.second);
    } else if (expr->is<Constant>()) {
      return expr;
    }
  }

  std::pair<ExprPtr, ExprPtr> compute() {


    auto single = ex<Constant>(0.0);
    if (singles_) {
      auto C = ex<Tensor>(
          L"C",
          std::initializer_list<Index>{
              Index::make_tmp_index(IndexSpace::instance(singles_gg_space))},
          std::initializer_list<Index>{Index::make_tmp_index(
              IndexSpace::instance(IndexSpace::other_unoccupied))});
      auto E_pa = ex<FNOperator>(
          std::initializer_list<Index>{C->as<Tensor>().bra()[0]},
          std::initializer_list<Index>{C->as<Tensor>().ket()[0]});
      auto C_Epa = C * E_pa;
      auto anti_herm_C =  C_Epa - adjoint(C_Epa);
      single = single + anti_herm_C;
      }

    if (ansatz_ == 2) {
      auto h = H(false);
      auto r = R12(gg_space);
      auto r_1 = R12(gg_space);

      ExprPtr A = ex<Constant>(0.0);
      if (doubles_) {
        A = A + (r - adjoint(r_1)) + single;
        simplify(A);
      } else {
        A = A + single;
        non_canon_simplify(A);
      }
        auto full_single_com = do_com(h,A);
        simplify(full_single_com);
        full_single_com = post_wick_sim(full_single_com);
        auto com_1 = simplification::hamiltonian_based_projector_2(full_single_com);
        auto one_body = com_1.first;
        auto two_body = com_1.second;
      // double commutator in eq. 9. Chem. Phys. 136, 084107 (2012).
      if (comutator_order_ >= 2) {
        auto full_double_com = ex<Constant>(1. / 2) * do_com(do_com(F(),relable(A)),relable(A));
        simplify(full_double_com);
        full_double_com = post_wick_sim(full_double_com);
        auto double_com = simplification::fock_based_projector_2(full_double_com);
        one_body = one_body + double_com.first;
        two_body = two_body + double_com.second;
        if(comutator_order_ >= 3){
          auto full_triple_com =  ex<Constant>(1. / 6) * do_com(do_com(do_com(F(),A),relable(A)),relable(A));
          simplify(full_triple_com);
          full_triple_com = post_wick_sim(full_triple_com);
          auto triple_com = simplification::fock_based_projector_2(full_triple_com);
          one_body = one_body + triple_com.first;
          two_body = two_body + triple_com.second;
          if (comutator_order_ >= 4){
            auto full_quadruple_com = ex<Constant>(1. / 24) * do_com(do_com(do_com(do_com(F(),A),relable(A)),relable(A)),relable(A));
            simplify(full_quadruple_com);
            full_quadruple_com = post_wick_sim(full_quadruple_com);
            auto quadruple_com = simplification::fock_based_projector_2(full_quadruple_com);
            one_body = one_body + quadruple_com.first;
            two_body = two_body + quadruple_com.second;
            if (comutator_order_ >= 5){
              auto full_quintuple_com = ex<Constant>(1. / 120) * do_com(do_com(do_com(do_com(do_com(F(),A),relable(A)),relable(A)),relable(A)),relable(A));
              simplify(full_quintuple_com);
              full_quintuple_com = post_wick_sim(full_quintuple_com);
              auto quintuple_com = simplification::fock_based_projector_2(full_quintuple_com);
              one_body = one_body + quintuple_com.first;
              two_body = two_body + quintuple_com.second;
              if (comutator_order_ >= 6){
                auto full_six_com = ex<Constant>(1. / 720) * do_com(do_com(do_com(do_com(do_com(do_com(F(),A),relable(A)),relable(A)),relable(A)),relable(A)),relable(A));
                simplify(full_six_com);
                full_six_com = post_wick_sim(full_six_com);
                auto six_com = simplification::fock_based_projector_2(full_six_com);
                one_body = one_body + six_com.first;
                two_body = two_body + six_com.second;
              }
            }
          }
        }
      }
      non_canon_simplify(one_body);
      non_canon_simplify(two_body);
      int term_count = 0;
      if(!one_body->is<Sum>()){
        term_count += 1;
      }
      else {
        for (auto i = 0; i < one_body->as<Sum>().summands().size(); i++) {
          term_count += 1;
        }
      }
      if(!two_body->is<Sum>()){
        term_count += 1;
      }
      else {
        for (auto i = 0; i < two_body->as<Sum>().summands().size(); i++) {
          term_count += 1;
        }
      }
      std::cout << "number of terms: " << term_count << std::endl;

      if (print_) {
        std::wcout << "one body terms: " << to_latex_align(one_body, 20, 2)
                   << std::endl;
        std::wcout << "two body terms: " << to_latex_align(two_body, 20, 2)
                   << std::endl;
      }
      return std::pair<ExprPtr, ExprPtr>{one_body, two_body};
    }
    // If we use the 2 body approximation, all terms with Density fall out since
    // they will happen to contain off diagonal G elements. we would get the
    // same result if we kept the decomposition and simplified, but this should
    // save time.
    if (ansatz_ == 1) {
      auto h = H(false);
      auto r = R12(gg_space);
      auto r_1 = R12(gg_space);

      auto A = r - adjoint(r);
      auto H_A = do_wick(ex<Constant>(1.) * ((h * A) - (A * h)));
      auto H_A_3 = keep_up_to_2_body_terms(H_A);

      H_A_3 = simplification::overlap_with_obs(H_A_3);
      H_A_3 = H_A_3 + ex<Constant>(0.);
      H_A_3 = simplification::screen_F12_proj(H_A_3, 1);
      H_A_3 = simplification::detail::tens_to_FNOps(H_A_3);
      simplify(H_A_3);
      auto com_1 = simplification::hamiltonian_based_projector_1(H_A_3);

      auto fFF = ex<Constant>(1. / 2) * compute_double_com(F(), r, r_1, 1);
      simplify(fFF);
      auto fFFt =
          ex<Constant>(1. / 2) *
          compute_double_com(F(), r, ex<Constant>(-1.) * adjoint(r_1), 1);
      simplify(fFFt);
      auto fFtFt = ex<Constant>(1. / 2) *
                   compute_double_com(F(), ex<Constant>(-1.) * adjoint(r),
                                      ex<Constant>(-1.) * adjoint(r_1), 1);
      simplify(fFtFt);
      auto fFtF =
          ex<Constant>(1. / 2) *
          compute_double_com(F(), ex<Constant>(-1.) * adjoint(r), r_1, 1);
      simplify(fFtF);

      auto fFF_sim = simplification::fock_based_projector_1(fFF);
      // std::wcout << "FF: " << to_latex_align(fFF_sim.second,20,2) <<
      // std::endl;
      auto fFFt_sim = simplification::fock_based_projector_1(fFFt);
      // std::wcout << "FFt: " << to_latex_align(fFFt_sim.second,20,2) <<
      // std::endl;
      auto fFtFt_sim = simplification::fock_based_projector_1(fFtFt);
      // std::wcout << "FtFt: " << to_latex_align(fFtFt_sim.second,20,2) <<
      // std::endl;
      auto fFtF_sim = simplification::fock_based_projector_1(fFtF);
      // std::wcout << "FtF: " << to_latex_align(fFtF_sim.second,20,2) <<
      // std::endl;

      auto one_body = com_1.first + (fFF_sim.first + fFFt_sim.first +
                                     fFtFt_sim.first + fFtF_sim.first);
      auto two_body = com_1.second + (fFF_sim.second + fFFt_sim.second +
                                      fFtFt_sim.second + fFtF_sim.second);
      simplify(one_body);
      simplify(two_body);
      int term_count = 0;
      for (auto i = 0; i < one_body->as<Sum>().summands().size(); i++) {
        term_count += 1;
      }
      for (auto i = 0; i < two_body->as<Sum>().summands().size(); i++) {
        term_count += 1;
      }
      std::cout << "number of terms: " << term_count << std::endl;

      if (print_) {
        std::wcout << "one body terms: " << to_latex_align(one_body, 20, 2)
                   << std::endl;
        std::wcout << "two body terms: " << to_latex_align(two_body, 20, 2)
                   << std::endl;
      }
      return std::pair<ExprPtr, ExprPtr>{one_body, two_body};
    }
  }
};

#endif  // SEQUANT_DOMAIN_SINGLE_REF_UCCF12_H
