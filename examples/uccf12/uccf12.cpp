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
  std::locale::global(std::locale("en_US.UTF-8"));
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(false);
  std::wcerr.sync_with_stdio(false);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  sequant::detail::OpIdRegistrar op_id_registrar;

  sequant::set_default_context(SeQuant(Vacuum::Physical, IndexSpaceMetric::Unit,
                                       BraKetSymmetry::conjugate));
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
    auto h = H();
    auto r = R12(gg_space);
    auto hr_comm =
        do_wick(ex<Constant>(2) *
                (h * r - r * h));  // this assumes symmetrization includes 1/2

    std::wcout << "[H,R] = " << to_latex_align(hr_comm, 20) << std::endl;

    auto hr_comm_12 = keep_1_and_2_body_terms(hr_comm);
    std::wcout << "[H,R]_{1,2} = " << to_latex_align(hr_comm_12, 20)
               << std::endl;
  }

  // double commutator, also needs symmetrization
  {
    auto f = F();
    auto r = R12(gg_space);
    auto fr_comm = do_wick(ex<Constant>(2) * (f * r - r * f));

    std::wcout << "[F,R] = " << to_latex_align(fr_comm, 20) << std::endl;

    {
      auto r = R12(gg_space);  // second instance of R
      auto a = r - adjoint(r);
      auto comm2 = do_wick(fr_comm * a - a * fr_comm);

      std::wcout << "[[F,R],A] = " << to_latex_align(comm2, 20) << std::endl;

      auto comm2_12 = keep_1_and_2_body_terms(comm2);
      std::wcout << "[[[F,R],A]_{1,2} = " << to_latex_align(comm2_12, 20)
                 << std::endl;
    }
  }
}
