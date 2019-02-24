#include <clocale>
#include <iostream>
#include "../src/SeQuant2/mbpt/sr/sr.hpp"

using namespace sequant2;

int main( int argc, char* argv[] )
{
  std::setlocale(LC_ALL,"en_US.UTF-8");
  std::cout.precision(std::numeric_limits<double>::max_digits10);
  sequant2::IndexSpace::register_standard_instances();
  sequant2::detail::OpIdRegistrar op_id_registrar;
  TensorCanonicalizer::set_cardinal_tensor_labels({L"A", L"f", L"g", L"t"});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  debug_canonicalize = true;

  using namespace sequant2::mbpt::sr::so;

  // H2**T2**T2 -> R2
  auto result = vac_av( A<2>() * H2() * T_<2>() * T_<2>(), {{1, 2}, {1, 3}} );
  std::wcout << "H2**T2**T2 -> R2 = " << to_latex_align(result, 20)
             << std::endl;
  assert(result->size() == 4);

  return 0;
}

