//
// Created by Eduard Valeyev on 3/20/18.
//

#define CATCH_CONFIG_RUNNER
#include <clocale>
#include "../../src/SeQuant/op.hpp"
#include "../../src/SeQuant/space.hpp"
#include "../../src/SeQuant/utility.hpp"
#include "../../src/domain/mbpt/convention.hpp"
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
  detail::OpIdRegistrar op_id_registrar;

  mbpt::set_default_convention();

  // uncomment to enable verbose output ...
  //Logger::set_instance(1);
  // ... or can instead selectively set/unset particular logging flags
  //Logger::get_instance().wick_contract = true;

  int result = session.run(argc, argv);

  return result;
}
