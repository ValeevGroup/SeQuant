#include "antisymmetrizer.hpp"
#include <SeQuant/core/wick.hpp>
#include <iostream>
#include "../sequant_setup.hpp"

using namespace sequant;

void try_main() {

  std::wcout << "START ANTISYMM_TEST: " << std::endl;
  const auto cumulant =ex<Tensor>(L"\\lambda", WstrList{L"i_1"}, WstrList{L"a_1"});
 // const auto a =ex<Tensor>(L"a",WstrList{L"i_2", L"i_3"},WstrList{L"a_2", L"a_3"});
  const auto a = ex<FNOperator>(std::initializer_list<Index>({Index(L"i_2"), Index(L"i_3")}), std::initializer_list<Index>({Index(L"a_2"), Index(L"a_3")}));
  auto a_cumulant = cumulant * a;
  std::wcout << "a_cumulant " << to_latex_align(a_cumulant) << std::endl;
  antisymmetrize _a_cumulant(a_cumulant);
  std::wcout << to_latex_align(_a_cumulant.result) << std::endl;

  auto cumulant2 = ex<Tensor>(L"\\lambda", WstrList{L"i_2"}, WstrList{L"a_2"});
  auto cumulant3 = ex<Tensor>(L"\\lambda", WstrList{L"i_3"}, WstrList{L"a_3"});
  auto cumulant_3x = cumulant * cumulant2 * cumulant3;
  std::wcout << "cumulant_3x " << to_latex_align(cumulant_3x) << std::endl;
  antisymmetrize _cumulant_3x(cumulant_3x);
  std::wcout << to_latex_align(_cumulant_3x.result) << std::endl;

  auto a1 = ex<Tensor>(L"a", WstrList{L"i_1"}, WstrList{L"a_1"});
  auto a1_cumu1_cumu2 = a1 * cumulant2 * cumulant3;
  std::wcout << "a1 y1 y2 " << to_latex_align(a1_cumu1_cumu2) << std::endl;
  antisymmetrize _a1_cumu1_cumu2(a1_cumu1_cumu2);
  std::wcout << to_latex_align(_a1_cumu1_cumu2.result) << std::endl;

  auto two_body_cumu = ex<Tensor>(L"\\lambda", WstrList{L"i_2", L"i_3"}, WstrList{L"a_2", L"a_3"});
  auto a1_cumu2 = a1 * two_body_cumu;
  std::wcout << " a1 y2 " << to_latex_align(a1_cumu2) << std::endl;
  antisymmetrize _a1_cumu2(a1_cumu2);
  std::wcout << to_latex_align(_a1_cumu2.result) << std::endl;

  auto cumu1_cumu2 = cumulant * two_body_cumu;
  std::wcout << " y1 y2 " << to_latex_align(cumu1_cumu2) << std::endl;
  antisymmetrize _cumu1_cumu2(cumu1_cumu2);
  std::wcout << to_latex_align(_cumu1_cumu2.result) << std::endl;

  auto cumu3 = ex<Tensor>(L"\\lambda", WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"});
  std::wcout << " y3 " << to_latex_align(cumu3) << std::endl;
  antisymmetrize _cumu3(cumu3);
  std::wcout << to_latex_align(_cumu3.result) << std::endl;

  std::wcout << "END ANTISYMM TEST: " << std::endl << std::endl << std::endl;

  auto cumu_to_dens = [](ExprPtr ex_){
    assert(ex_->is<Tensor>());
    assert(ex_->as<Tensor>().rank() == 1);
    assert(ex_->as<Tensor>().label() == L"\\lambda");
    auto down_0 = ex_->as<Tensor>().ket()[0];
    auto up_0 = ex_->as<Tensor>().bra()[0];

    auto density = ex<Tensor>(L"\\gamma",std::initializer_list<Index>{up_0}, std::initializer_list<Index>{down_0});
    return density;
  };

  auto cumu2_to_density = [&](ExprPtr ex_){
    assert(ex_->is<Tensor>());
    assert(ex_->as<Tensor>().rank() == 2);
    assert(ex_->as<Tensor>().label() == L"\\lambda");

    auto down_0 = ex_->as<Tensor>().ket()[0];
    auto up_0 = ex_->as<Tensor>().bra()[0];
    auto down_1 = ex_->as<Tensor>().ket()[1];
    auto up_1 = ex_->as<Tensor>().bra()[1];

    auto density2 = ex<Tensor>(L"\\gamma",std::initializer_list<Index>{up_0, up_1}, std::initializer_list<Index>{down_0, down_1});
    auto density_1 = ex<Tensor>(L"\\gamma",std::initializer_list<Index>{up_0}, std::initializer_list<Index>{down_0});
    auto density_2 = ex<Tensor>(L"\\gamma",std::initializer_list<Index>{up_1}, std::initializer_list<Index>{down_1});

    auto d1_d2 = antisymmetrize(density_1 * density_2);
    return density2 * ex<Constant>(-1) * d1_d2.result;
  };

  auto cumu3_to_density = [&](ExprPtr ex_){
    assert(ex_->is<Tensor>());
    assert(ex_->as<Tensor>().rank() == 3);
    assert(ex_->as<Tensor>().label() == L"\\lambda");

    auto down_0 = ex_->as<Tensor>().ket()[0];
    auto up_0 = ex_->as<Tensor>().bra()[0];
    auto down_1 = ex_->as<Tensor>().ket()[1];
    auto up_1 = ex_->as<Tensor>().bra()[1];
    auto down_2 = ex_->as<Tensor>().ket()[2];
    auto up_2 = ex_->as<Tensor>().bra()[2];

    auto cumulant2 = ex<Tensor>(L"\\lambda",std::initializer_list<Index>{up_1, up_2}, std::initializer_list<Index>{down_1, down_2});
    auto density_1 = ex<Tensor>(L"\\gamma",std::initializer_list<Index>{up_0}, std::initializer_list<Index>{down_0});
    auto density_2 = ex<Tensor>(L"\\gamma",std::initializer_list<Index>{up_1}, std::initializer_list<Index>{down_1});
    auto density_3 = ex<Tensor>(L"\\gamma",std::initializer_list<Index>{up_2}, std::initializer_list<Index>{down_2});
    auto density3 = ex<Tensor>(L"\\gamma",std::initializer_list<Index>{up_0, up_1,up_2}, std::initializer_list<Index>{down_0, down_1, down_2});

    auto d1_d2 = antisymmetrize(density_1 * density_2 * density_3 + density_1 * cumulant2);
    auto temp_result = density3 * ex<Constant>(-1) * d1_d2.result;

    for (auto&& product : temp_result->as<Sum>().summands()){
      for (auto&& factor : product->as<Product>().factors()){
        if (factor->is<Tensor>() && (factor->as<Tensor>().label() == L"\\lambda") && (factor->as<Tensor>().rank() == 2)){
          factor = cumu2_to_density(factor);
        }
      }
    }
    for (auto&& product : temp_result->as<Sum>().summands()){
      for (auto&& factor : product->as<Product>().factors()){
        if (factor->is<Tensor>() && factor->as<Tensor>().label() == L"\\lambda" && factor->as<Tensor>().rank() == 1){
          factor = cumu_to_dens(factor);
        }
      }
    }
    return temp_result;
  };

  auto one_body_sub = [](ExprPtr ex_){//J. Chem. Phys. 132, 234107 (2010); https://doi.org/10.1063/1.3439395 eqn 15 for
    assert(ex_->is<FNOperator>());
    assert(ex_->as<FNOperator>().rank() == 1);
    auto down_0 = ex_->as<FNOperator>().annihilators()[0].index();
    auto up_0 = ex_->as<FNOperator>().creators()[0].index();

    const auto a = ::ex<FNOperator>(std::initializer_list<Index>{up_0}, std::initializer_list<Index>{down_0});
    const auto cumu1 = ::ex<Tensor>(L"\\lambda", WstrList{up_0.label()}, WstrList{down_0.label()});

    return (a + ex<Constant>(-1) * cumu1);
  };

  auto two_body_decomp = [](ExprPtr ex_){//J. Chem. Phys. 132, 234107 (2010); https://doi.org/10.1063/1.3439395 eqn 16 for \tilde{a^{pr}_{qs}}
    assert(ex_->is<FNOperator>());
    assert(ex_->as<FNOperator>().rank() == 2);

    auto down_0 = ex_->as<FNOperator>().annihilators()[0].index();
    auto down_1 = ex_->as<FNOperator>().annihilators()[1].index();

    auto up_0 = ex_->as<FNOperator>().creators()[0].index();
    auto up_1 = ex_->as<FNOperator>().creators()[1].index();

    const auto cumu1 = ::ex<Tensor>(L"\\lambda", WstrList{up_0.label()}, WstrList{down_0.label()});
    const auto cumu2 = ::ex<Tensor>(L"\\lambda", WstrList{up_1.label()}, WstrList{down_1.label()});
    const auto a = ::ex<FNOperator>(std::initializer_list<Index>{up_1}, std::initializer_list<Index>{down_1});
    const auto double_cumu = ::ex<Tensor>(L"\\lambda", WstrList{up_0.label(), up_1.label()}, WstrList{down_0.label(), down_1.label()});

    auto term1 = cumu1 * a;
    auto term2 = cumu1 * cumu2;
    auto term3 = double_cumu;

    auto sum_of_terms = antisymmetrize( term1 + term2 + term3);
    return ex<Constant>(-1) * sum_of_terms.result;

  };

  //express 3-body term as sums of 1 and 2-body term. as described in J. Chem. Phys. 132, 234107 (2010); https://doi.org/10.1063/1.3439395 eqn 17.
  auto three_body_decomp = [&](ExprPtr ex_){
    assert(ex_->is<FNOperator>());
    assert(ex_->as<FNOperator>().rank() == 3);

    auto down_0 = ex_->as<FNOperator>().annihilators()[0].index();
    auto down_1 = ex_->as<FNOperator>().annihilators()[1].index();
    auto down_2 = ex_->as<FNOperator>().annihilators()[2].index();

    auto up_0 = ex_->as<FNOperator>().creators()[0].index();
    auto up_1 = ex_->as<FNOperator>().creators()[1].index();
    auto up_2 = ex_->as<FNOperator>().creators()[2].index();

    const auto cumulant = ::ex<Tensor>(L"\\lambda", WstrList{up_0.label()}, WstrList{down_0.label()});
    const auto a = ::ex<FNOperator>(std::initializer_list<Index>{up_1, up_2}, std::initializer_list<Index>{down_1, down_2});
    auto a_cumulant = cumulant * a;


    auto cumulant2 = ::ex<Tensor>(L"\\lambda", WstrList{up_1.label()}, WstrList{down_1.label()});
    auto cumulant3 = ::ex<Tensor>(L"\\lambda", WstrList{up_2.label()}, WstrList{down_2.label()});
    auto cumulant_3x = cumulant * cumulant2 * cumulant3;

    auto a1 = ::ex<FNOperator>( std::initializer_list<Index>{up_0}, std::initializer_list<Index>{down_0});
    auto a1_cumu1_cumu2 = a1 * cumulant2 * cumulant3;

    auto two_body_cumu = ::ex<Tensor>(L"\\lambda", WstrList{up_1.label(), up_2.label()}, WstrList{down_1.label(),down_2.label()});
    auto a1_cumu2 = a1 * two_body_cumu;

    auto cumu1_cumu2 = cumulant * two_body_cumu;

    auto cumu3 = ::ex<Tensor>(L"\\lambda", WstrList{up_0.label(), up_1.label(), up_2.label()}, WstrList{down_0.label(), down_1.label(),down_2.label()});

    auto sum_of_terms = antisymmetrize(a_cumulant + cumulant_3x + a1_cumu1_cumu2 + a1_cumu2 + cumu1_cumu2 + cumu3);

    auto temp_result = sum_of_terms.result;

    for (auto&& product : temp_result->as<Sum>().summands()){//replace all the two body terms with one body terms.
      if (product->is<Product>()) {
        for (auto&& factor : product->as<Product>().factors()) {
          if (factor->is<FNOperator>() &&
              factor->as<FNOperator>().rank() == 2) {
            factor = two_body_decomp(factor);
          }
        }
      }
    }

    for (auto&& product : temp_result->as<Sum>().summands()){//replace the one body terms with the substituted expression
      if(product->is<Product>()) {
        for (auto&& factor : product->as<Product>().factors()) {
          if (factor->is<FNOperator>() &&
              factor->as<FNOperator>().rank() == 1) {
            factor = one_body_sub(factor);
          }
        }
      }
    }
    for (auto&& product : temp_result->as<Sum>().summands()) {
      if (product->is<Product>()) {
        for(auto&& factor : product->as<Product>().factors()) {
          if (factor->is<Tensor>()) {
            if (factor->as<Tensor>().label() == L"\\lambda" &&
                factor->as<Tensor>().rank() == 3) {
              factor = cumu3_to_density(factor);
            } else if (factor->as<Tensor>().label() == L"\\lambda" &&
                       factor->as<Tensor>().rank() == 2) {
              factor = cumu2_to_density(factor);
            } else if (factor->as<Tensor>().label() == L"\\lambda") {
              factor = cumu_to_dens(factor);
            }
            else{assert(factor->as<Tensor>().label() != L"\\lambda");}
          }
        }
      }
    }
    return temp_result;
  };

  auto three_body_substitution = [&](ExprPtr &input){
    assert(input->is<Sum>());
    for(auto&& product : input->as<Sum>().summands()){
      if(product->is<Product>()){
        for (auto&& factor : product->as<Product>().factors()){
          if (factor->is<FNOperator>() && (factor->as<FNOperator>().rank() == 3)){ // find the 3-body terms
            factor = three_body_decomp(factor); // decompose that term and replace the existing term.
          }
        }
      }
    }
  input->canonicalize();
  simplify(input);

    for(auto&& product : input->as<Sum>().summands()){
      if(product->is<Product>()){
        for (auto&& factor : product->as<Product>().factors()){
          assert(factor->is<FNOperator>() | factor->is<Tensor>() | factor->is<Constant>());
          if (factor->is<Tensor>() && (factor->as<Tensor>().label() == L"\\lambda")){ // find the 3-body terms
            factor = cumu_to_dens(factor); // decompose that term and replace the existing term.
          }
        }
      }
    }
    return input;
  };

  auto a3 = ex<FNOperator>(std::initializer_list<Index>({Index(L"i_1"), Index(L"i_2"), Index(L"i_3")}), std::initializer_list<Index>({Index(L"a_1"), Index(L"a_2"), Index(L"a_3")}));
  a3 = ex<Constant>(1) * a3 + ex<Constant>(1);
  auto new_a3 = three_body_substitution(a3);
  new_a3->canonicalize();
  //simplify(new_a3);
  std::wcout << to_latex_align(new_a3) << std::endl;

}

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
  sequant::set_default_context(SeQuant(Vacuum::Physical, IndexSpaceMetric::Unit, BraKetSymmetry::conjugate /*,Spinbasis::spin_orbit*/));
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
