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
  //set_num_threads(1);

  try {
    try_main();
  }
  catch(...) {
    std::cerr << "caught an unknown exception, ouch";
  }

  return 0;
}

void
try_main() {

  {
    auto result = vac_av( H() * T(2), {{0, 1}} );

    std::wcout << "H*T12 -> 0 = " << to_latex_align(result, 20)
               << std::endl;
  }

}