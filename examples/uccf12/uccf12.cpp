#include <SeQuant/core/op.hpp>
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
  // set_num_threads(1);

  try {
    try_main();
  } catch (std::exception& ex) {
    std::cerr << "caught a std::exception: " << ex.what();
  } catch (...) {
    std::cerr << "caught an unknown exception, ouch";
  }

  return 0;
}

void try_main() {
  auto do_wick = [](ExprPtr expr) {
    using sequant::FWickTheorem;
    FWickTheorem wick{expr};
    wick.spinfree(false).full_contractions(false);
    auto result = wick.compute();
    simplify(result);
    return result;
  };

  // this keeps 1- and 2-body terms w.r.t. physical vacuum
  auto keep_1_and_2_body_terms = [](const ExprPtr& input) {
    assert(input->is<Sum>());
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
              keep = (rank >= 1 && rank <= 2);
            }
          }
          return !keep;
        });
    auto result = ex<Sum>(ranges::begin(filtered_summands),
                          ranges::end(filtered_summands));
    return result;
  };

  auto gg_space =
      IndexSpace::active_occupied;  // Geminal-generating space: active
                                    // occupieds is the normal choice, all
                                    // orbitals is the reference-independent
                                    // (albeit expensive) choice

  // single commutator, needs symmetrization
  {
    auto h = H(false);
    auto r = R12(gg_space);
    auto hr_comm = do_wick( ex<Constant>(2) * (h*r - r*h) );  // this assumes symmetrization includes 1/2

    std::wcout << "[H,R] = " << to_latex_align(hr_comm, 20) << std::endl;

    auto hr_comm_12 = keep_1_and_2_body_terms(hr_comm);
    std::wcout << "[H,R]_{1,2} = " << to_latex_align(hr_comm_12, 20)
               << std::endl;
  }

  // double commutator, also needs symmetrization
  {
    auto f = F();
    auto r = R12(gg_space);
    auto a = ex<Constant>(0.5) * (r - adjoint(r));
    auto fa_comm = do_wick( ex<Constant>(1) * (f*a - a*f) );

    {
      auto r = R12(gg_space);  // second instance of R
      auto a = ex<Constant>(0.5) * (r - adjoint(r));
      auto comm2 = do_wick(ex<Constant>(0.5) * (fa_comm * a - a * fa_comm));


      auto comm2_12 = keep_1_and_2_body_terms(comm2);
      std::wcout << "[[[F,A],A]_{1,2} = " << to_latex_align(comm2_12, 20)
                 << std::endl;
    }
    // compute terms individually
    auto r_ = R12(gg_space);// need to make second instance.

    auto fr = f * r;
    auto frr = do_wick(ex<Constant>(1) * (fr * r_));
    auto frr_12 = keep_1_and_2_body_terms(frr);
    std::wcout <<  "frr term =  " << to_latex_align(frr_12,20) << std::endl;

    auto rf = do_wick(ex<Constant>(1)* r * f);
    auto rfr = do_wick(ex<Constant>(1) * (rf * r_));
    auto rfr_12 = keep_1_and_2_body_terms(rfr);
    std::wcout <<  "rfr term =  " << to_latex_align(rfr_12,20) << std::endl;

    auto rr = do_wick(ex<Constant>(1)* r * r_);
    auto rrf = do_wick(ex<Constant>(1) * (rr * f));
    auto rrf_12 = keep_1_and_2_body_terms(rrf);
    std::wcout <<  "rrf term =  " << to_latex_align(rrf_12,20) << std::endl;

    auto rr_dag = do_wick(ex<Constant>(1) * (r * adjoint(r_)));
    auto frr_dag = do_wick(ex<Constant>(1) * f * (rr_dag));
    auto frr_dag_12 = keep_1_and_2_body_terms(frr_dag);
    std::wcout <<  "$frr_{adj}$ term =  " << to_latex_align(frr_dag_12,20) << std::endl;

    auto r_dagf = do_wick(ex<Constant>(1) * adjoint(r) * f);
    auto r_dagfr = do_wick(ex<Constant>(1) * r_dagf * r_);
    auto r_dagfr_12 = keep_1_and_2_body_terms(r_dagfr);
    std::wcout <<  "$r_{adj}fr $term =  " << to_latex_align(r_dagfr_12,20) << std::endl;

    auto aut_rf = do_wick(ex<Constant>(1) * r_ * f);
    auto rfr_dag = do_wick(ex<Constant>(1) * aut_rf * adjoint(r));
    auto rfr_dag_12 = keep_1_and_2_body_terms(rfr_dag);
    std::wcout <<  "$rfr_{adj}$ term =  " << to_latex_align(rfr_dag_12,20) << std::endl;

    auto rr_dagf = do_wick(ex<Constant>(1) * rr_dag * f);
    auto rr_dagf_12 = keep_1_and_2_body_terms(rr_dagf);
    std::wcout <<  "$rr_{adj}f $term =  " << to_latex_align(rr_dagf_12,20) << std::endl;

    auto r_dagr_dag = do_wick(ex<Constant>(1) * adjoint(r) * adjoint(r_));
    auto r_dagr_dagf = do_wick(ex<Constant>(1) * r_dagr_dag * f);
    auto r_dagr_dagf_12 = keep_1_and_2_body_terms(r_dagr_dagf);
    std::wcout <<  "$r_{adj}r_{adj}f$ term =  " << to_latex_align(r_dagr_dagf_12,20) << std::endl;


    auto r_dagr = do_wick(ex<Constant>(1) * adjoint(r_) * r);
    auto r_dagrf = do_wick(ex<Constant>(1) * r_dagr * f);
    auto r_dagrf_12 = keep_1_and_2_body_terms(r_dagrf);
    std::wcout <<  "$r_{adj}rf$ =  " << to_latex_align(r_dagrf_12,20) << std::endl;

    auto fr_dag = do_wick(ex<Constant>(1) * f * adjoint(r_));
    auto fr_dagr = do_wick(ex<Constant>(1) * fr_dag * r);
    auto fr_dagr_12 = keep_1_and_2_body_terms(fr_dagr);
    auto fr_dag_12 = keep_1_and_2_body_terms(fr_dag);
    std::wcout <<  "$fr_{adj}r$ =  " << to_latex_align(fr_dagr_12,20) << std::endl;

  }
}
