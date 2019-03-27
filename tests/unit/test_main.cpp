//
// Created by Eduard Valeyev on 3/20/18.
//

#define CATCH_CONFIG_RUNNER
#include <clocale>
#include "../../src/SeQuant/op.hpp"
#include "../../src/SeQuant/space.hpp"
#include "../../src/SeQuant/utility.hpp"
#include "../../src/domain/mbpt/op.hpp"
#include "catch.hpp"

int main( int argc, char* argv[] )
{
  using namespace std;
  using namespace sequant;

  Catch::Session session;

  // global setup...
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(false);
  std::wcerr.sync_with_stdio(false);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  IndexSpace::register_standard_instances();
  detail::OpIdRegistrar op_id_registrar;

  // load the default cardinal tensor labels for MBPT ...
  // tensors that do not appear on this list will appear after Tensors with these
  // labels
  TensorCanonicalizer::set_cardinal_tensor_labels(mbpt::cardinal_tensor_labels);

  // uncomment to enable verbose output ...
  //Logger::set_instance(1);
  // ... or can instead selectively set/unset particular logging flags
  //Logger::get_instance().wick_contract = true;

  int result = session.run(argc, argv);

  return result;
}
