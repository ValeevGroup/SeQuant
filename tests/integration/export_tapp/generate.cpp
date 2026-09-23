// Generates the C code exercised by driver.c:
//   export_tapp_generate <output file> real|complex

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/tapp.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/io/shorthands.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace sequant;

namespace {

ExportNode<> to_tree(std::wstring_view spec) {
  const auto res =
      deserialize<ResultExpr>(spec, {.def_perm_symm = Symmetry::Nonsymm,
                                     .def_braket_symm = BraKetSymmetry::Nonsymm,
                                     .def_col_symm = ColumnSymmetry::Nonsymm});
  return to_export_tree(res);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <output file> real|complex\n";
    return EXIT_FAILURE;
  }

  const std::string_view kind = argv[2];
  if (kind != "real" && kind != "complex") {
    std::cerr << "Unknown scalar type '" << kind << "'\n";
    return EXIT_FAILURE;
  }
  const bool complex = kind == "complex";

  auto registry = std::make_shared<IndexSpaceRegistry>();
  registry->add(L"i", 0b01, is_vacuum_occupied, is_reference_occupied, is_hole,
                3);
  registry->add(L"a", 0b10, is_particle, 4);
  auto resetter =
      set_scoped_default_context({.index_space_registry_shared_ptr = registry});

  TAPPGeneratorContext ctx;
  ctx.set_prefix(std::string(kind));
  ctx.set_scalar_type(complex ? TAPPScalarType::Complex : TAPPScalarType::Real);
  ctx.set_tag(registry->retrieve(L"i"), "o");
  ctx.set_tag(registry->retrieve(L"a"), "v");
  ctx.set_dim(registry->retrieve(L"i"), "nocc");
  ctx.set_dim(registry->retrieve(L"a"), "nvirt");

  std::vector<ExpressionGroup<>> groups;
  if (complex) {
    groups.emplace_back("energy");
    groups.back().add(to_tree(L"E = s^* w{i1,i2;a1,a2} t{a1,a2;i1,i2}"));
  } else {
    groups.emplace_back("residual");
    groups.back().add(
        to_tree(L"R{a1,a2;i1,i2} = 1/2 g{a1,a2;a3,a4} t{a3,a4;i1,i2} + "
                L"f{a1;a3} t{a3,a2;i1,i2} + s v{a1,a2;i1,i2}"));
    // Reuses the result of the previous expression
    groups.back().add(
        to_tree(L"S{a1,a2;i1,i2} = R{a1,a2;i1,i2} - R{a2,a1;i1,i2}"));

    // Both expressions share a single contraction plan
    groups.emplace_back("singles");
    groups.back().add(to_tree(L"X{a1;i1} = f{a1;a2} p{a2;i1}"));
    groups.back().add(to_tree(L"Y{a1;i1} = 2 f{a1;a2} q{a2;i1}"));

    groups.emplace_back("energy");
    groups.back().add(
        to_tree(L"E = 1/4 w{i1,i2;a1,a2} t{a1,a2;i1,i2} + "
                L"f{i1;a1} t{a1,a2;i1,i2} u{i2;a2}"));
  }

  TAPPGenerator<> generator;
  export_groups(std::move(groups), generator, ctx);

  std::ofstream out(argv[1]);
  out << generator.get_generated_code();
  if (!out) {
    std::cerr << "Failed to write '" << argv[1] << "'\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
