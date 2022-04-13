//
// Created by Nakul Teke on 4/13/22.
//
// Spin-orbital pCCSDT example
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/eqs/cceqs.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <clocale>

using namespace sequant;

#define runtime_assert(tf)                                             \
  if (!(tf)) {                                                         \
    std::ostringstream oss;                                            \
    oss << "failed assert at line " << __LINE__ << " in SRCC example"; \
    throw std::runtime_error(oss.str().c_str());                       \
  }

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

  using sequant::eqs::compute_all;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  auto cc_r = sequant::eqs::cceqvec{3, 3}(true, true, true, true, true);

}
