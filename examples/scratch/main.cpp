#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <clocale>
#include <iostream>

/*
 * Scratch main()
 */

int main(int argc, char* argv[]) {
  using namespace sequant;
  // global sequant setup...
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  detail::OpIdRegistrar op_id_registrar;

  mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  Logger::get_instance().wick_stats = false;
  auto tnsr1 = ex<Tensor>(
      Tensor(L"g", {L"i_1", L"i_2"}, {L"a_1", L"a_2"}, Symmetry::antisymm));
  canonicalize(tnsr1);  // tnsr1->canonicalize();
  auto tnsr2 = ex<Tensor>(
      Tensor(L"g", {L"i_2", L"i_1"}, {L"a_1", L"a_2"}, Symmetry::antisymm));
  canonicalize(tnsr2);  // ->canonicalize();

  // std::wcout << "canon1 = " << (canon1 ? canon1->to_latex() : L"null")
  //            << "\ncanon2 = " << (canon2 ? canon2->to_latex() : L"null")
  //            << std::endl;
  std::wcout << "tnsr1 = " << tnsr1->to_latex()
             << " tnsr2 = " << tnsr2->to_latex() << std::endl;
  return 0;
}
