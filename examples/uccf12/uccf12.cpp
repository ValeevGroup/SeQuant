#include <SeQuant/core/wick.hpp>
#include "../sequant_setup.hpp"
#include "../antisymmetrizer_test/three_body_decomp.hpp"
#include <locale>
#include <codecvt>
#include "simplifications.h"


using namespace sequant;

void try_main();

int main(int argc, char* argv[]) {
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

  sequant::set_default_context(SeQuant(Vacuum::Physical, IndexSpaceMetric::Unit, BraKetSymmetry::conjugate,
                                       SPBasis::spinfree));
  mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  // WARNING some code is not thread safe ...
  //set_num_threads(1);

  try {
    try_main();
  }
  catch(std::exception& ex) {
    std::cerr << "caught a std::exception: " << ex.what();
  }
  catch(...) {
    std::cerr << "caught an unknown exception, ouch";
  }

  return 0;
}

void
try_main() {

  auto do_wick = [](ExprPtr expr) {
    using sequant::FWickTheorem;
    FWickTheorem wick{expr};
    wick.spinfree(false)
        .full_contractions(false);
    auto result = wick.compute();
   simplify(result);
    return result;
  };

  auto keep_1_body_terms = [](const ExprPtr& input){
    assert(input->is<Sum>());
    auto filtered_summands =
        input->as<Sum>().summands() | ranges::views::remove_if([](const ExprPtr &ptr) {
          assert(ptr->is<Product>());
          bool keep = false;
          bool found_operator = false;
          for(auto&& factor : ptr->as<Product>().factors()) {
            if (factor->is<FNOperator>()) {
              assert(!found_operator);
              found_operator = true;
              const auto rank = factor->as<FNOperator>().rank();
              keep = (rank <= 1);
            }
          }
          return !keep;
        });
    auto result = ex<Sum>(ranges::begin(filtered_summands),
                          ranges::end(filtered_summands));
    return result;
  };
  // this keeps 1- and 2-body terms w.r.t. physical vacuum
  auto keep_1_and_2_body_terms = [](const ExprPtr& input) {
    assert(input->is<Sum>());
    auto filtered_summands =
        input->as<Sum>().summands() | ranges::views::remove_if([](const ExprPtr &ptr) {
          assert(ptr->is<Product>());
          bool keep = false;
          bool found_operator = false;
          for(auto&& factor : ptr->as<Product>().factors()) {
            if (factor->is<FNOperator>()) {
              assert(!found_operator);
              found_operator = true;
              const auto rank = factor->as<FNOperator>().rank();
              keep = (rank >= 1 && rank <=2);
            }
          }
          return !keep;
        });
    auto result = ex<Sum>(ranges::begin(filtered_summands),
                          ranges::end(filtered_summands));
    return result;
  };

  auto keep_up_to_3_body_terms = [](const ExprPtr& input) {
    if(input->is<Sum>()) {
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
    }
    else return input;
  };

  auto gg_space = IndexSpace::active_occupied;  // Geminal-generating space: active occupieds is the normal choice, all orbitals is the reference-independent (albeit expensive) choice

  auto remove_const = [](const ExprPtr ex_){
    auto new_expression = ex<Constant>(0);
    if (ex_->is<Sum>()){
      for (auto&& product : ex_->as<Sum>().summands()){
        bool has_fnop = false;
        for (auto&& factor : product->as<Product>().factors()){
          if (factor->is<FNOperator>()){
            has_fnop = true;
          }
        }
        if (has_fnop){ new_expression = new_expression + product;}
      }
    }
    return new_expression;
  };
    // based on the expressions in the paper, the tensor in the A operator never exists with lower indicies unoccupied and upper indicies occoupied
    auto revert_tens_adj = [](ExprPtr ex_, std::wstring label){
    if (ex_->is<Sum>()){
      for (auto&& product : ex_->as<Sum>().summands()){
        for (auto&& factor : product->as<Product>().factors()){
          if(factor->is<Tensor>() && factor->as<Tensor>().label() == label){
            factor = adjoint(factor);
          }
        }
      }
    }
    else if (ex_->is<Product>()){
      for(auto&& factor : ex_->as<Product>()){
        if(factor->is<Tensor>() && factor->as<Tensor>().label() == label){
          factor = adjoint(factor);
        }
      }
    }
    else if(ex_->is<Product>() && ex_->as<Tensor>().label() == label){
      ex_ = adjoint(ex_);
    }
    return ex_;
  };
    auto new_idx = [](int pos, Index idx){
      auto new_number = std::to_wstring(pos + 6);
      auto new_string = IndexSpace::base_key(idx.space()) + L'_' + new_number;
      auto new_label = Index{new_string};
      return new_label;
    };
    auto replace_idx = [](ExprPtr ex_, Index og, Index newer){
        assert(ex_->is<Product>());
        auto constant = ex_->as<Product>().scalar();
        auto new_product = ex<Constant>(1);
        for(auto&& factor : ex_->as<Product>().factors()){
          if(factor->is<Tensor>()){
            std::vector<Index> new_bras;
            for(auto&& bra : factor->as<Tensor>().bra()){
              if(bra.label() == og.label()){
                new_bras.push_back(newer);
              }
              else{new_bras.push_back(bra);}
            }
            std::vector<Index> new_kets;
            for(auto&& ket : factor->as<Tensor>().ket()){
              if(ket.label() == og.label()){
                new_kets.push_back(newer);
              }
              else{new_kets.push_back(ket);}
            }
            auto new_tensor = ex<Tensor>(factor->as<Tensor>().label(), new_bras, new_kets);
            new_product = new_tensor * new_product;
          }
          if(factor->is<FNOperator>()){
            std::vector<Index> new_cres;
            for (auto&& cre : factor->as<FNOperator>().creators()){
              if(cre.index().label() == og.label()){
                new_cres.push_back(newer);
              }
              else{new_cres.push_back(cre.index());}
            }
            std::vector<Index> new_anns;
            for (auto&& ann : factor->as<FNOperator>().annihilators()){
              if(ann.index().label() == og.label()){
                new_anns.push_back(newer);
              }
              else{new_anns.push_back(ann.index());}
            }
            auto new_op = ex<FNOperator>(new_cres,new_anns);
            new_product = new_product * new_op;
          }
        }
        auto result = (ex<Constant>(constant) * new_product);
        return result->as<Product>();
    };
    auto op_to_tens = [](ExprPtr ex_){
      assert(ex_->is<FNOperator>());
      //vector upper
      //adjoint lower
      //creators lower
      // annihilators upper.
      //bra_lower
      //ket_upper
      std::vector<Index> bra_indices;
      std::vector<Index> ket_indices;
      for(auto&& cre : ex_->as<FNOperator>().creators()){
        bra_indices.push_back(cre.index());
      }
      for(auto&& ann : ex_->as<FNOperator>().annihilators()){
        ket_indices.push_back(ann.index());
      }
      auto label = get_default_context().spbasis() == SPBasis::spinfree ? L"E" : L"a";
      auto result = ex<Tensor>(label, bra_indices, ket_indices);
      return result;
    };
    auto tens_to_op = [](ExprPtr ex_) {
      assert(ex_->is<Tensor>());
      auto result = ex<FNOperator>(ex_->as<Tensor>().ket(),ex_->as<Tensor>().bra());
      return result;
    };
    auto screen_F_tensors = [](ExprPtr ex_) {// F tensors must contain contain indices in the bra with space > all. this includes complete, completeunoccupied, and inactiveunoccupied.
      assert(ex_->is<Tensor>());
      assert(ex_->as<Tensor>().label() == L"F");
      bool good = false;
      for (auto&& bra : ex_->as<Tensor>().bra()){
        if(bra.space().type() == IndexSpace::complete || bra.space().type() == IndexSpace::complete_unoccupied || bra.space().type() == IndexSpace::inactive_unoccupied){
          good = true;
        }
      }

      for (auto&& ket : ex_->as<Tensor>().ket()){
        if(ket.space().type() == IndexSpace::complete || ket.space().type() == IndexSpace::complete_unoccupied || ket.space().type() == IndexSpace::inactive_unoccupied){
          good = true;
        }
      }
      if(good){
        return ex_;
      }
      else{
        return ex<Constant>(0);
      }
    };
    auto screen_densities = [](ExprPtr ex_){// densities probably should be non-zero if each index has a chance to be occupied, in other words, screen out densities containing unoccupied labels.
      assert(ex_->is<Tensor>());
      assert(ex_->as<Tensor>().label() == L"\\Gamma");
      bool good = true;
      for (auto&& bra : ex_->as<Tensor>().bra()){
        if (bra.space().type() == IndexSpace::unoccupied || bra.space().type() == IndexSpace::complete_unoccupied){
          good = false;
        }
      }
      for (auto&& ket : ex_->as<Tensor>().ket()){
        if (ket.space().type() == IndexSpace::unoccupied || ket.space().type() == IndexSpace::complete_unoccupied){
          good = false;
        }
      }
      if(good){
        return ex_;
      }
      else{return ex<Constant>(0);}
    };
  auto treat_fock = [&] (ExprPtr ex_){
    auto new_ex_ = ex<Constant>(0);
    for (auto&& product : ex_->as<Sum>().summands()){
      auto new_product = ex<Constant>(1);
      auto scalar = product->as<Product>().scalar();
      bool contains_density = false;
      for (auto&& factor : product->as<Product>().factors()){
         if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\Gamma"){
           contains_density = true;
         }
         if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"f" && contains_density && (factor->as<Tensor>().bra()[0].space().type() != IndexSpace::complete_unoccupied || factor->as<Tensor>().ket()[0].space().type() != IndexSpace::complete_unoccupied)){
          auto overlap = make_overlap(factor->as<Tensor>().bra()[0],factor->as<Tensor>().ket()[0]);
          new_product = overlap * factor * new_product;
        }
        /* else if(factor->is<Tensor>() && (factor->as<Tensor>().label() == L"E" || factor->as<Tensor>().label() == L"a")){
          new_product = new_product * tens_to_op(factor);
        }*/
        else new_product = new_product * factor;
      }
      new_product = ex<Constant>(scalar) * new_product;
      std::wcout << to_latex_align(new_product) << std::endl;
      non_canon_simplify(new_product);
      new_ex_ = new_ex_ + new_product;
    }
    //std::wcout << to_latex_align(new_ex_,20, 3) << std::endl;
    non_canon_simplify(new_ex_);
    FWickTheorem wick{new_ex_};
    wick.reduce(new_ex_);
    std::wcout << to_latex_align(new_ex_,20, 3) << std::endl;
    non_canon_simplify(new_ex_);
    return new_ex_;
  };
  // single commutator, needs symmetrization
  {
    auto h = H(false);
    std::wcout << "H = " << to_latex_align(h, 20)<< std::endl;
    auto r = R12(gg_space);
    std::wcout << "r = " << to_latex_align(r, 20)<< std::endl;
    auto hr_comm = do_wick( ex<Constant>(1) * (h*r - r*h) );  // this assumes symmetrization includes 1/2
    //auto ht_comm = do_wick(ex<Constant>(1) * (h*r - r*h));
    std::wcout << "[H,R] = " << to_latex_align(hr_comm, 20)
               << std::endl;

    auto hr_comm_12 = keep_up_to_3_body_terms(hr_comm);
    auto cumu_decomp = decompositions::three_body_substitution(hr_comm_12,2); // apply the cumulant decompostition to the expression
    simplify(cumu_decomp);

    auto compare_h = simplification::hamiltonian_based(cumu_decomp);
    std::wcout << "[H,R]12" << to_latex_align(compare_h, 20, 3) << std::endl;

    //end of of first commutator. now comes the comutator one at a time, i.e. [[F,A],A] = [[F,R]_12, R]_12 + complex conjugates and their versions. if we get this term we can interpret the rest from properties of comutators.
    auto com_f_a = do_wick(ex<Constant>(0.5) * (F() * adjoint(r) - adjoint(r) * F()));
    auto com_f_a_3 = keep_up_to_3_body_terms(com_f_a);
    auto com_f_a_2 = decompositions::three_body_substitution(com_f_a_3,2);
    simplify(com_f_a_2);
    auto r_ = R12(gg_space);
    auto comr_ = do_wick(com_f_a_2 * r_);
    auto r_com = do_wick(r_ * com_f_a_2);
    auto double_com = comr_ - r_com;
    simplify(double_com);
    auto double_com_3 = keep_up_to_3_body_terms(double_com);
    auto double_com_2 = decompositions::three_body_substitution(double_com_3,2);
    simplify(double_com_2);
    auto compare_f = simplification::fock_based(double_com_2);
    std::wcout << "[[H,R]12,R]12" << to_latex_align(compare_f,20,3) << std::endl;

    auto two_body = do_wick(ex<Constant>(0.5) * (F() * adjoint(r) - adjoint(r) * F()));
    two_body = keep_1_and_2_body_terms(two_body);
    auto asdfa = do_wick(two_body * r_);
    auto lkuahsdf = do_wick(r_ * two_body);
    auto duble_2_body = asdfa - lkuahsdf;
    simplify(duble_2_body);
    duble_2_body = keep_1_and_2_body_terms(duble_2_body);
    //duble_2_body = simplification::fock_based(duble_2_body);
    std::wcout << "terms without densities" << to_latex_align(duble_2_body,20,3) << std::endl;

    /*
    //now I need to assert deltas between each E op. such that every E^{ }_{ } -> E^{p}_{q} and E^{}{}_{}{} -> E^{p}{r}_{q}{s}
    // maximal contractions from E^{q}{s}_{p}{r} E^{k1}{k2}_{k3}{k4} E^{q}{s}_{p}{r} -> d^p_k1 d^r_k2 d^q_k3 d^s_k4 E^{q}{s}_{p}{r}. adjoint leads to correct result
    auto E2_ = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_4"), Index(L"p_3")}), std::initializer_list<Index>({Index(L"p_1"), Index(L"p_2")}));
    auto E1_ = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_2")}), std::initializer_list<Index>({Index(L"p_1")}));
    auto transform_E2 = keep_1_and_2_body_terms(do_wick(H2(false) * E2_));
    simplify(transform_E2);
    std::wcout << to_latex_align(E2_) << std::endl;
    simplify(E2_);
    auto H2_E2 = H2(false) * E2_;
    simplify(H2_E2);
    std::wcout << to_latex_align(H2_E2) << std::endl;

    auto temper = E2_ * transform_E2;
    simplify(temper);
    auto fin = adjoint(keep_1_and_2_body_terms(do_wick(temper)));
    simplify(fin);
    std::wcout << "wick ErE = " << to_latex_align(fin, 20,3) << std::endl;
    auto transform_E1 = keep_1_body_terms(do_wick(H1() * E1_));
    simplify(transform_E1);
    auto temp_E1 = E1_ * transform_E1;
    simplify(temp_E1);
    auto fin_E1 = adjoint(keep_1_body_terms(do_wick(temp_E1)));
    std::wcout << " wick EhE = " << to_latex_align(fin_E1,20,3) << std::endl;
     */
  }

  // double commutator, also needs symmetrization
  {/*
    mbpt::sr::so::make_op f_(1,1,sequant::mbpt::OpType::f,false);
    auto ao_space = IndexSpace::all;
    auto unocc = IndexSpace::all;
    auto f = f_(unocc, ao_space,false);
    //std::wcout << "F: " << to_latex_align(f) << std::endl;
    auto r = revert_tens_adj(R12(gg_space),L"F");
    //std::wcout << "r: " << to_latex_align(r) << std::endl;
    auto a = ex<Constant>(0.5) * r - revert_tens_adj(adjoint(r), L"F");
    auto fa_comm = do_wick( ex<Constant>(1) * (f*a - a*f) );

    {
      auto r = revert_tens_adj(R12(gg_space),L"F");  // second instance of R
      auto a = ex<Constant>(0.5) * r - revert_tens_adj(adjoint(r), L"F");
      auto comm2 = do_wick(ex<Constant>(0.5) * (fa_comm * a - a * fa_comm));


      auto comm2_12 = keep_up_to_3_body_terms(comm2);
      //comm2_12 = decompositions::three_body_substitution(comm2_12,2);
      comm2_12->canonicalize();
      simplify(comm2_12);
      std::wcout << "[[[F,A],A]_{1,2} = " << to_latex_align(comm2_12, 30) << std::endl;
    }*/
    // compute terms individually
    /*auto r_ = R12(gg_space);// need to make second instance.

    auto fr = f * r;
    auto frr = do_wick(ex<Constant>(1) * (fr * r_));
    auto frr_12 = keep_up_to_3_body_terms(frr);
    frr_12 = decompositions::three_body_substitution(frr_12,2);
    std::wcout <<  "frr term =  " << to_latex_align(frr_12,20) << std::endl;

    auto rf = do_wick(ex<Constant>(1)* r * f);
    auto rfr = do_wick(ex<Constant>(1) * (rf * r_));
    auto rfr_12 = keep_up_to_3_body_terms(rfr);
    rfr_12 = decompositions::three_body_substitution(rfr,2);
    std::wcout <<  "rfr term =  " << to_latex_align(rfr_12,20) << std::endl;

    auto rr = do_wick(ex<Constant>(1)* r * r_);
    auto rrf = do_wick(ex<Constant>(1) * (rr * f));
    auto rrf_12 = keep_up_to_3_body_terms(rrf);
    rrf_12 = decompositions::three_body_substitution(rrf,2);
    std::wcout <<  "rrf term =  " << to_latex_align(rrf_12,20) << std::endl;

    auto rr_dag = do_wick(ex<Constant>(1) * (r * adjoint(r_)));
    auto frr_dag = do_wick(ex<Constant>(1) * f * (rr_dag));
    auto frr_dag_12 = keep_up_to_3_body_terms(frr_dag);
    frr_dag_12 = decompositions::three_body_substitution(frr_dag_12,2);
    std::wcout <<  "$frr_{adj}$ term =  " << to_latex_align(frr_dag_12,20) << std::endl;

    auto r_dagf = do_wick(ex<Constant>(1) * adjoint(r) * f);
    auto r_dagfr = do_wick(ex<Constant>(1) * r_dagf * r_);
    auto r_dagfr_12 = keep_up_to_3_body_terms(r_dagfr);
    r_dagfr_12 = decompositions::three_body_substitution(r_dagfr_12,2);
    std::wcout <<  "$r_{adj}fr $term =  " << to_latex_align(r_dagfr_12,20) << std::endl;

    auto aut_rf = do_wick(ex<Constant>(1) * r_ * f);
    auto rfr_dag = do_wick(ex<Constant>(1) * aut_rf * adjoint(r));
    auto rfr_dag_12 = keep_up_to_3_body_terms(rfr_dag);
    rfr_dag_12 = decompositions::three_body_substitution(rfr_dag_12,2);
    std::wcout <<  "$rfr_{adj}$ term =  " << to_latex_align(rfr_dag_12,20) << std::endl;

    auto rr_dagf = do_wick(ex<Constant>(1) * rr_dag * f);
    auto rr_dagf_12 = keep_up_to_3_body_terms(rr_dagf);
    rr_dagf_12 = decompositions::three_body_substitution(rr_dagf_12, 2);
    std::wcout <<  "$rr_{adj}f $term =  " << to_latex_align(rr_dagf_12,20) << std::endl;

    auto r_dagr_dag = do_wick(ex<Constant>(1) * adjoint(r) * adjoint(r_));
    auto r_dagr_dagf = do_wick(ex<Constant>(1) * r_dagr_dag * f);
    auto r_dagr_dagf_12 = keep_up_to_3_body_terms(r_dagr_dagf);
    r_dagr_dagf_12 = decompositions::three_body_substitution(r_dagr_dagf_12,2);
    std::wcout <<  "$r_{adj}r_{adj}f$ term =  " << to_latex_align(r_dagr_dagf_12,20) << std::endl;


    auto r_dagr = do_wick(ex<Constant>(1) * adjoint(r_) * r);
    auto r_dagrf = do_wick(ex<Constant>(1) * r_dagr * f);
    auto r_dagrf_12 = keep_up_to_3_body_terms(r_dagrf);
    r_dagrf_12 = decompositions::three_body_substitution(r_dagrf_12,2);
    std::wcout <<  "$r_{adj}rf$ =  " << to_latex_align(r_dagrf_12,20) << std::endl;

    auto fr_dag = do_wick(ex<Constant>(1) * f * adjoint(r_));
    auto fr_dagr = do_wick(ex<Constant>(1) * fr_dag * r);
    auto fr_dagr_12 = keep_up_to_3_body_terms(fr_dagr);
    fr_dagr_12 = decompositions::three_body_substitution(fr_dagr_12,2);
    auto fr_dag_12 = keep_up_to_3_body_terms(fr_dag);
    fr_dag_12 = decompositions::three_body_substitution(fr_dag_12,2);
    std::wcout <<  "$fr_{adj}r$ =  " << to_latex_align(fr_dagr_12,20) << std::endl;

    std::wcout << "FINAL EXPRESSION: " << std::endl;
    auto total_expression = fr_dag_12 + fr_dagr_12 + r_dagrf_12;*/

  }

}