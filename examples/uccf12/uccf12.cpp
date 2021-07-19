#include <SeQuant/core/wick.hpp>
#include "../sequant_setup.hpp"
#include "../antisymmetrizer_test/three_body_decomp.hpp"

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
              keep = (rank <=3);
            }
          }
          return !keep;
        });
    auto result = ex<Sum>(ranges::begin(filtered_summands),
                          ranges::end(filtered_summands));
    return result;
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
    auto enforce_obs = [&](ExprPtr ex_){
      assert(!ex_->is<Product>());
      auto E2 = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_1"), Index(L"p_2")}), std::initializer_list<Index>({Index(L"p_3"), Index(L"p_4")}));
      auto E2_ = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_5"), Index(L"p_6")}), std::initializer_list<Index>({Index(L"p_3"), Index(L"p_4")}));
      auto E1 = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_1")}), std::initializer_list<Index>({Index(L"p_3")}));
      auto E1_ = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_5")}), std::initializer_list<Index>({Index(L"p_3")}));
      int fn_rank;
      for (auto&& factor : ex_->as<Product>().factors()){
        if (factor->is<FNOperator>()){
          fn_rank = ex_->as<FNOperator>().rank();
        }
      }
      if(fn_rank == 2){
        auto transform_E2 = keep_1_and_2_body_terms(do_wick(E2 * ex_));
        auto temper = transform_E2 * E2_;
        simplify(temper);
        auto fin = keep_1_and_2_body_terms(do_wick(temper));// two terms are left and we only want one of them.
        for (auto&& product : fin->as<Sum>().summands()){
          for (auto&& factor : product->as<Product>().factors()){
            if(factor->is<FNOperator>() && factor->as<FNOperator>().creators()[0].index().label() == L"p1"){
              ex_ = product;
              break;
            }
          }
        }
      }
      if(fn_rank == 1){
        auto transform_E1 = keep_1_body_terms(do_wick(E1 * ex_));
        auto temper = transform_E1 * E1_;
        simplify(temper);
        ex_ = keep_1_body_terms(do_wick(temper));
      }
      return ex_;
    };
  // single commutator, needs symmetrization
  {
    auto h = H(false);
    std::wcout << "H = " << to_latex_align(h, 20)<< std::endl;
    auto r = R12(gg_space,2);
    std::wcout << "r = " << to_latex_align(r, 20)<< std::endl;
    auto hr_comm = do_wick( ex<Constant>(2) * (h*r - r*h) );  // this assumes symmetrization includes 1/2
    //auto ht_comm = do_wick(ex<Constant>(1) * (h*r - r*h));
    std::wcout << "[H,R] = " << to_latex_align(hr_comm, 20)
               << std::endl;

    auto hr_comm_12 = keep_up_to_3_body_terms(hr_comm);
    //auto cumu_decomp = decompositions::three_body_substitution(hr_comm_12,2);
    //simplify(cumu_decomp);
    //cumu_decomp = remove_const(cumu_decomp);
    /*for (auto&& product : cumu_decomp->as<Sum>().summands()){
      product = enforce_obs(product);
    }*/
    //now I need to assert deltas between each E op. such that every E^{ }_{ } -> E^{p}_{q} and E^{}{}_{}{} -> E^{p}{r}_{q}{s}
    // maximal contractions from E^{q}{s}_{p}{r} E^{k1}{k2}_{k3}{k4} E^{q}{s}_{p}{r} -> d^p_k1 d^r_k2 d^q_k3 d^s_k4 E^{q}{s}_{p}{r}. adjoint leads to correct result
    auto E2 = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_3"), Index(L"p_4")}), std::initializer_list<Index>({Index(L"p_1"), Index(L"p_2")}));
    auto E2_ = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_3"), Index(L"p_4")}), std::initializer_list<Index>({Index(L"p_1"), Index(L"p_2")}));
    auto E1 = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_1")}), std::initializer_list<Index>({Index(L"p_3")}));
    auto E1_ = ex<FNOperator>(std::initializer_list<Index>({Index(L"p_5")}), std::initializer_list<Index>({Index(L"p_3")}));
    auto transform_E2 = keep_1_and_2_body_terms(do_wick(r * E2_));
    std::wcout << "E2: " << to_latex_align(E2) << std::endl;
    auto sym = E2_->as<FNOperator>()._symmetry();
    for (auto&& product : transform_E2->as<Sum>().summands()){
      for(auto&& factor : product->as<Product>().factors()){
        if (factor->is<FNOperator>()){
          assert(factor->as<FNOperator>()._symmetry() == Symmetry::nonsymm);
        }
      }
    }
    simplify(transform_E2);
    //std::wcout << "ErE = " << to_latex_align(H2(false) * E2_, 20,3) << std::endl;
    for (auto&& product : transform_E2->as<Sum>().summands()){
      for(auto&& factor : product->as<Product>().factors()){
        if (factor->is<FNOperator>()){
          assert(factor->as<FNOperator>()._symmetry() == Symmetry::nonsymm);
        }
      }
    }
    std::wcout << " wick rE = " << to_latex_align(transform_E2, 20,3) << std::endl;
    auto temper = E2 * transform_E2;
    simplify(temper);
    std::wcout << "ErE = " << to_latex_align(temper, 20,3) << std::endl;
    for (auto&& product : temper->as<Sum>().summands()){
      for(auto&& factor : product->as<Product>().factors()){
        if (factor->is<FNOperator>()){
          assert(factor->as<FNOperator>()._symmetry() == Symmetry::nonsymm);
        }
      }
    }
    auto fin = adjoint(keep_1_and_2_body_terms(do_wick(temper)));
    for (auto&& product : fin->as<Sum>().summands()){
      for(auto&& factor : product->as<Product>().factors()){
        if (factor->is<FNOperator>()){
          assert(factor->as<FNOperator>()._symmetry() == Symmetry::nonsymm);
        }
      }
    }
    std::wcout << "wick ErE = " << to_latex_align(fin, 20,3) << std::endl;
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