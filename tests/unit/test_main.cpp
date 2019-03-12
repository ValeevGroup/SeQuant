//
// Created by Eduard Valeyev on 3/20/18.
//

#define CATCH_CONFIG_RUNNER
#include <clocale>
#include "../../src/SeQuant/op.hpp"
#include "../../src/SeQuant/space.hpp"
#include "../../src/SeQuant/utility.hpp"
#include "catch.hpp"

int main( int argc, char* argv[] )
{
  using namespace std;
  using namespace sequant;

  Catch::Session session;

  // global setup...
  setlocale(LC_ALL,"en_US.UTF-8");
  cout.precision(numeric_limits<double>::max_digits10);
  IndexSpace::register_standard_instances();
  detail::OpIdRegistrar op_id_registrar;

  // uncomment to enable verbose output ...
  //Logger::set_instance(1);
  // ... or can instead selectively set/unset particular logging flags
  //Logger::get_instance().wick_contract = true;

  int result = session.run(argc, argv);

  return result;
}
