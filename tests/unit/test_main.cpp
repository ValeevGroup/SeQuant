//
// Created by Eduard Valeyev on 3/20/18.
//

#define CATCH_CONFIG_RUNNER
#include <clocale>
#include "../../src/SeQuant/op.hpp"
#include "../../src/SeQuant/space.hpp"
#include "catch.hpp"

int main( int argc, char* argv[] )
{
  Catch::Session session;

  // global setup...
  std::setlocale(LC_ALL,"en_US.UTF-8");
  std::cout.precision(std::numeric_limits<double>::max_digits10);
  sequant::IndexSpace::register_standard_instances();
  sequant::detail::OpIdRegistrar op_id_registrar;

  int result = session.run(argc, argv);

  return result;
}
