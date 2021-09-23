#include "../../SeQuant/domain/transcorrelated/three_body_decomp.hpp"
#include "../../SeQuant/domain/transcorrelated/simplifications.h"
#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/sr/sr.hpp>

#include <clocale>

using namespace sequant;
using namespace sequant::mbpt::sr::so;


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
  auto compute_double_com = [&](ExprPtr e1, ExprPtr e2, ExprPtr e3){
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
  };

  auto gg_space = IndexSpace::occupied;  // Geminal-generating space: active occupieds is the normal choice, all orbitals is the reference-independent (albeit expensive) choice
  //start transformation
  {
    auto h = H(false);
    std::wcout << "H = " << to_latex_align(h, 20)<< std::endl;
    auto r = R12(gg_space);
    auto r_1 = R12(gg_space);
    std::wcout << "r = " << to_latex_align(r, 20)<< std::endl;
    //way 1
    auto A = r - adjoint(r);
    auto H_A = do_wick((h * A) - (A * h));
    auto H_A_adj = do_wick((h * adjoint(r)) - (adjoint(r_1) * h));
     auto H_A_3 = keep_up_to_3_body_terms(H_A);
    auto H_A_2 = decompositions::three_body_substitution(H_A_3,2);
    simplify(H_A_2);
    auto com_1 = simplification::hamiltonian_based(H_A_2);

    std::wcout << "h A one body: " << to_latex_align(com_1.first,20,2) << std::endl;
    std::wcout << "h A two body: " << to_latex_align(com_1.second,20,2) << std::endl;

    auto fFF = compute_double_com(F(),r,r_1);
    auto fFFt = compute_double_com(F(),r,ex<Constant>(-1.) * adjoint(r_1));
    auto fFtFt = compute_double_com(F(),ex<Constant>(-1.) * adjoint(r),ex<Constant>(-1.) * adjoint(r_1));
    auto fFtF = compute_double_com(F(),ex<Constant>(-1.) * adjoint(r),r_1);

    auto fFF_sim = simplification::fock_based(fFF);
   // std::wcout << "FF: " << to_latex_align(fFF_sim.second,20,2) << std::endl;
    auto fFFt_sim = simplification::fock_based(fFFt);
    std::wcout << "FFt one body: " << to_latex_align(fFFt_sim.first,20,2) << std::endl;
    std::wcout << "FFt two body: " << to_latex_align(fFFt_sim.second,20,2) << std::endl;
    auto fFtFt_sim = simplification::fock_based(fFtFt);
    //std::wcout << "FtFt: " << to_latex_align(fFtFt_sim.second,20,2) << std::endl;
    auto fFtF_sim = simplification::fock_based(fFtF);
    std::wcout << "FtF one body: " << to_latex_align(fFtF_sim.first,20,2) << std::endl;
    std::wcout << "FtF two body: " << to_latex_align(fFtF_sim.second,20,2) << std::endl;


    auto one_body = com_1.first + ex<Constant>(0.5) * (fFF_sim.first + fFFt_sim.first + fFtFt_sim.first + fFtF_sim.first);
    auto two_body = com_1.second + ex<Constant>(0.5) * (fFF_sim.second + fFFt_sim.second + fFtFt_sim.second + fFtF_sim.second);
    non_canon_simplify(one_body);
    non_canon_simplify(two_body);
    std::wcout << "one body terms: " << to_latex_align(one_body,20,2) << std::endl;
    std::wcout << "two body terms: " << to_latex_align(two_body,20,2) << std::endl;

  }//end transformation

}