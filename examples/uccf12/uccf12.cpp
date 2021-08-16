#include <SeQuant/core/wick.hpp>
#include "../sequant_setup.hpp"
#include "../../SeQuant/domain/transcorrelated/three_body_decomp.hpp"
#include <locale>
#include <codecvt>
#include "../../SeQuant/domain/transcorrelated/simplifications.h"


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
  //start transformation
  {
    auto h = H(false);
    std::wcout << "H = " << to_latex_align(h, 20)<< std::endl;
    auto r = R12(gg_space);
    std::wcout << "r = " << to_latex_align(r, 20)<< std::endl;

    //start single comutator compute.
    auto hr_comm = do_wick( ex<Constant>(1) * (h*r - r*h) );
    std::wcout << "[H,R] = " << to_latex_align(hr_comm, 20)
               << std::endl;

    auto hr_comm_12 = keep_up_to_3_body_terms(hr_comm);
    auto cumu_decomp = decompositions::three_body_substitution(hr_comm_12,2); // apply the cumulant decompostition to the expression
    simplify(cumu_decomp);

    auto compare_h = simplification::hamiltonian_based(cumu_decomp);// simplify the expression to remove zero terms not caught by sequant core.
    std::wcout << "$[H,R]_12$" << to_latex_align(compare_h, 20, 3) << std::endl;
    //end single comutator compute

    //end of of first commutator. now comes the second comutator one at a time, i.e. (1/2)[[f,A],A] = (1/2)[[f,R]_12, R]_12 + complex conjugates and their versions. if we get this term we can interpret the rest from properties of comutators.
    auto com_f_a = do_wick(ex<Constant>(0.5) * (F() * adjoint(r) - adjoint(r) * F()));
    auto com_f_a_3 = keep_up_to_3_body_terms(com_f_a);
    auto com_f_a_2 = decompositions::three_body_substitution(com_f_a_3,2);
    simplify(com_f_a_2);
    //start double comutator
    auto r_ = R12(gg_space); //a new unique expression pointer is needed here.
    //odd, but wick needs these terms evaluated separatly
    auto comr_ = do_wick(com_f_a_2 * r_);
    auto r_com = do_wick(r_ * com_f_a_2);
    auto double_com = comr_ - r_com;
    simplify(double_com);
    auto double_com_3 = keep_up_to_3_body_terms(double_com);
    auto double_com_2 = decompositions::three_body_substitution(double_com_3,2);
    simplify(double_com_2);
    auto compare_f = simplification::fock_based(double_com_2); // fock matrix based simplifcation which has additional screening of fock matrix from Brillouin's Theory
    //end double comutator.
    std::wcout << "$[[H,adj(R)]_12,R]_12$" << to_latex_align(compare_f,20,3) << std::endl;

  }//end transformation

}