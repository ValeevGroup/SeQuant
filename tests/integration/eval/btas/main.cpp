//
// Created by Bimal Gaudel on 7/31/21.
//
#include <clocale>
#include <iostream>

#include <btas/btas.h>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <btas/scf_btas.hpp>
#include <calc_info.hpp>

#define runtime_assert(tf)                                                  \
  if (!(tf)) {                                                              \
    std::ostringstream oss;                                                 \
    oss << "failed assert at line " << __LINE__ << " in eval_btas example"; \
    throw std::runtime_error(oss.str().c_str());                            \
  }

///
/// excitation (2 = ccsd (default) through 6 supported)
/// fock(/eri) tensor data file name (default "fock.dat"/"eri.dat")
/// eri(/fock) tensor data file name (default "fock.dat"/"eri.dat")
/// \note Format of data files:
/// ----------------------------
/// size_t size_t size_t         # rank, nocc, nvirt
/// double                       # data ------
/// ...                          # data       |
/// ...                          # ....       |  no. of double entries =
/// (nocc+nvirt)^rank
/// ...                          # data       |
/// double                       # data ------
/// ----------------------------
/// \note The rank of fock tensor is 2
/// \note The rank of eri tensor is 4
/// \note The nocc and nvirt in both files must match
///
/// \todo Use markdown file formats such as yml for input of calculation params.
///
int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "\nHelp:\n"
              << "<executable> <config_file> <fock_or_eri_file> "
                 "<eri_or_fock_file> [<output_file>]"
              << "\n\n"
              << "Config file format\n"
              << "----\n"
              << sequant::eval::ParseConfigFile{}.help() << "\n----\n";
    // return 1;
  }

  sequant::set_locale();

  using namespace sequant;
  detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      {.index_space_registry_shared_ptr = mbpt::make_min_sr_spaces(),
       .vacuum = Vacuum::SingleProduct,
       .canonicalization_options =
           CanonicalizeOptions::default_options().copy_and_set(
               CanonicalizationMethod::Complete)});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  TensorCanonicalizer::set_cardinal_tensor_labels(
      sequant::mbpt::cardinal_tensor_labels());
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_minimal_registry()});

  // for optimization tests, set occupied and unoccupied index extents
  {
    auto reg = get_default_context().mutable_index_space_registry();
    auto occ = reg->retrieve_ptr(L"i");
    auto uocc = reg->retrieve_ptr(L"a");
    SEQUANT_ASSERT(occ);
    SEQUANT_ASSERT(uocc);
    occ->approximate_size(10);
    uocc->approximate_size(100);
    SEQUANT_ASSERT(uocc->approximate_size() == 100);
  }

  std::string calc_config = argc > 1 ? argv[1] : "calc.inp";
  std::string fock_file = argc > 2 ? argv[2] : "fock_so.dat";
  std::string eri_file = argc > 3 ? argv[3] : "eri_so.dat";
  // not yet implemented:
  std::string out_file = argc > 4 ? argv[4] : "";
  //

  auto const calc_info =
      eval::make_calc_info(calc_config, fock_file, eri_file, out_file);
  auto scf_btas =
      eval::btas::SequantEvalScfBTAS<btas::Tensor<double>>{calc_info};
  scf_btas.scf(std::wcout);

  double const expected{-0.07068045196165902};
  double const threshold = calc_info.scf_opts.conv;
  double const ediff = std::fabs(expected - scf_btas.energy());

  runtime_assert((ediff <= threshold));
  return 0;
}
