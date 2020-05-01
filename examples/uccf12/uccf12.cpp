#include "../../src/SeQuant/core/wick.hpp"
#include "../sequant_setup.hpp"

using namespace  sequant;

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

  mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  set_num_threads(1);

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

  {
    auto h = H();
    auto x = R12();
    auto hx_comm = do_wick( h*x - x*h );

    std::wcout << "[H,F] = " << to_latex_align(hx_comm, 20)
               << std::endl;

    // keep only 1 and 2-body ints
    assert(hx_comm->is<Sum>());
    auto filtered_summands =
        hx_comm->as<Sum>().summands() | ranges::views::remove_if([](const ExprPtr &ptr) {
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
    auto hx_comm_12 = ex<Sum>(ranges::begin(filtered_summands),
        ranges::end(filtered_summands));
    std::wcout << "[H,F]_{1,2} = " << to_latex_align(hx_comm_12, 20)
               << std::endl;
  }

}