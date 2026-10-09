#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include "test_export.hpp"

#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/export_expr.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/export/expression_group.hpp>
#include <SeQuant/core/export/generation_optimizer.hpp>
#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/export/julia_itensor.hpp>
#include <SeQuant/core/export/julia_tensor_kit.hpp>
#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/export/python_einsum.hpp>
#include <SeQuant/core/export/reordering_context.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include "catch2_sequant.hpp"

#include <boost/algorithm/string.hpp>

#include <filesystem>
#include <fstream>
#include <optional>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

using namespace sequant;

std::vector<std::vector<std::size_t>> twoElectronIntegralSymmetries() {
  // Symmetries of spin-summed (skeleton) two-electron integrals
  return {
      // g^{pq}_{rs}
      {0, 1, 2, 3},
      // g^{ps}_{rq}
      {0, 3, 2, 1},
      // g^{rq}_{ps}
      {2, 1, 0, 3},
      // g^{rs}_{pq}
      {2, 3, 0, 1},

      // g^{qp}_{sr}
      {1, 0, 3, 2},
      // g^{qr}_{sp}
      {1, 2, 3, 0},
      // g^{sp}_{qr}
      {3, 0, 1, 2},
      // g^{sr}_{qp}
      {3, 2, 1, 0},
  };
}

std::vector<std::filesystem::path> enumerate_export_tests() {
#ifndef SEQUANT_UNIT_TESTS_SOURCE_DIR
#error \
    "Need source dir for unit tests in order to locate directory for export test files"
#endif
  std::filesystem::path base_dir = SEQUANT_UNIT_TESTS_SOURCE_DIR;
  base_dir /= "export_tests";

  if (!std::filesystem::is_directory(base_dir)) {
    throw sequant::Exception("Invalid base dir for export tests");
  }

  std::vector<std::filesystem::path> files;

  for (std::filesystem::path current :
       std::filesystem::directory_iterator(base_dir)) {
    if (current.extension() == ".export_test") {
      files.push_back(std::move(current));
    }
  }

  return files;
}

// clang-format off
using KnownGenerators = std::tuple<
    TextGenerator<TextGeneratorContext>,
    JuliaITensorGenerator<JuliaITensorGeneratorContext>,
    JuliaTensorKitGenerator<JuliaTensorKitGeneratorContext>,
    JuliaTensorOperationsGenerator<JuliaTensorOperationsGeneratorContext>,
    NumPyEinsumGenerator,
    PyTorchEinsumGenerator,
    ItfGenerator<ItfContext>
>;
// clang-format on

template <typename Generator>
std::string get_format_name() {
  Generator g;

  return g.get_format_name();
}

template <typename... Generator>
std::set<std::string> known_format_names(std::tuple<Generator...>) {
  std::set<std::string> names;

  (names.insert(get_format_name<Generator>()), ...);

  return names;
}

void configure_context_defaults(TextGeneratorContext &) {}

void configure_context_defaults(ItfContext &ctx) {
  auto registry = get_default_context().index_basis_registry();
  IndexSpace occ = registry->retrieve("i");
  IndexSpace virt = registry->retrieve("a");
  IndexSpace aux = registry->retrieve("x");
  IndexSpace act = registry->retrieve("u");

  ctx.set_tag(occ, "c");
  ctx.set_tag(virt, "e");
  ctx.set_tag(aux, "x");
  ctx.set_tag(act, "a");

  ctx.set_name(occ, "Closed");
  ctx.set_name(virt, "External");
  ctx.set_name(aux, "Auxiliary");
  ctx.set_name(act, "Active");
}

void configure_context_defaults(JuliaTensorOperationsGeneratorContext &ctx) {
  auto registry = get_default_context().index_basis_registry();
  IndexSpace occ = registry->retrieve("i");
  IndexSpace virt = registry->retrieve("a");
  IndexSpace aux = registry->retrieve("x");

  ctx.set_dim(occ, "nocc");
  ctx.set_dim(virt, "nv");
  ctx.set_dim(aux, "naux");

  ctx.set_tag(occ, "o");
  ctx.set_tag(virt, "v");
  ctx.set_tag(aux, "x");
}

void configure_context_defaults(NumPyEinsumGeneratorContext &ctx) {
  auto registry = get_default_context().index_basis_registry();
  IndexSpace occ = registry->retrieve("i");
  IndexSpace virt = registry->retrieve("a");
  IndexSpace aux = registry->retrieve("x");

  ctx.set_shape(occ, "nocc");
  ctx.set_shape(virt, "nvirt");
  ctx.set_shape(aux, "naux");

  ctx.set_tag(occ, "o");
  ctx.set_tag(virt, "v");
  ctx.set_tag(aux, "x");
}

void configure_context_defaults(PyTorchEinsumGeneratorContext &ctx) {
  auto registry = get_default_context().index_basis_registry();
  IndexSpace occ = registry->retrieve("i");
  IndexSpace virt = registry->retrieve("a");
  IndexSpace aux = registry->retrieve("x");

  ctx.set_shape(occ, "nocc");
  ctx.set_shape(virt, "nvirt");
  ctx.set_shape(aux, "naux");

  ctx.set_tag(occ, "o");
  ctx.set_tag(virt, "v");
  ctx.set_tag(aux, "x");
}

void add_to_context(TextGeneratorContext &ctx, std::string_view key,
                    std::string_view value) {
  if (key == "batch_indices") {
    std::vector<Index> indices;
    for (auto current : std::ranges::views::split(value, ',')) {
      Index idx(std::string(current.begin(), current.end()));
      indices.push_back(idx);
    }

    ctx.set_batch_indices(indices);
  } else {
    throw Exception(
        "Unsupported key in TextGenerator context specification: '" +
        std::string(key) + "'");
  }
}
void add_to_context(ItfContext &ctx, std::string_view key,
                    std::string_view value) {
  if (key == "batch_indices") {
    std::vector<Index> indices;
    for (auto current : std::ranges::views::split(value, ',')) {
      Index idx(std::string(current.begin(), current.end()));
      indices.push_back(idx);
    }

    ctx.set_batch_indices(indices);
  } else {
    throw Exception("Unsupported key in ITF context specification: '" +
                    std::string(key) + "'");
  }
}
void add_to_context(PythonEinsumGeneratorContext &, std::string_view,
                    std::string_view) {
  // PythonEinsumGeneratorContext doesn't support context specifications
}
void add_to_context(JuliaTensorOperationsGeneratorContext &ctx,
                    std::string_view key, std::string_view value) {
  auto parse_space_map = [](std::string_view spec) {
    auto pos = spec.find("->");
    if (pos == std::string_view::npos) {
      throw sequant::Exception("Malformed space map");
    }

    std::string space(spec.substr(0, pos));
    std::string map(spec.substr(pos + 2));

    boost::trim(space);
    boost::trim(map);

    return std::make_pair(
        get_default_context().index_basis_registry()->retrieve(space),
        std::string(map));
  };

  if (key == "tag") {
    auto [space, tag] = parse_space_map(value);
    ctx.set_tag(space, tag);
  } else if (key == "dim") {
    auto [space, dim] = parse_space_map(value);
    ctx.set_dim(space, dim);
  } else {
    throw sequant::Exception(
        "Unsupported key in Julia context specification '" + std::string(key) +
        "'");
  }
}

template <typename Context>
void add_to_context(Context &ctx, const std::string &line) {
  auto pos = line.find(":");
  if (pos == std::string::npos) {
    throw sequant::Exception(
        "Malformed context specification: missing ':' in '" + line + "'");
  }

  std::string key = line.substr(0, pos);
  boost::trim(key);
  std::string value = line.substr(pos + 1);
  boost::trim(value);

  if (key.empty()) {
    throw sequant::Exception("Malformed context specification: Empty key");
  }
  if (value.empty()) {
    throw sequant::Exception("Malformed context specification: Empty value");
  }

  add_to_context(ctx, key, value);
}

std::vector<ExpressionGroup<>> parse_expression_spec(const std::string &spec,
                                                     bool retain_braket) {
  std::vector<ExpressionGroup<>> groups;

  std::istringstream in(spec);

  for (std::string line; std::getline(in, line);) {
    boost::trim(line);
    if (line.empty()) {
      continue;
    }

    if (line.starts_with("section") && line.ends_with(":")) {
      std::string name = line.substr(7, line.size() - 7 - 1);
      boost::trim(name);
      SEQUANT_ASSERT(!name.empty());
      groups.emplace_back(std::move(name));
      continue;
    }

    if (groups.empty()) {
      groups.emplace_back();
    }

    try {
      ResultExpr res = deserialize<ResultExpr>(
          line, {.def_perm_symm = Symmetry::Nonsymm,
                 .def_braket_symm = BraKetSymmetry::Nonsymm,
                 .def_col_symm = ColumnSymmetry::Nonsymm});
      groups.back().add(to_export_tree(res, retain_braket));
    } catch (...) {
      ExprPtr expr =
          deserialize(line, {.def_perm_symm = Symmetry::Nonsymm,
                             .def_braket_symm = BraKetSymmetry::Nonsymm,
                             .def_col_symm = ColumnSymmetry::Nonsymm});
      groups.back().add(to_export_tree(expr, retain_braket));
    }
  }

  return groups;
}

TEMPLATE_LIST_TEST_CASE("export_tests", "[export]", KnownGenerators) {
  using CurrentGen = TestType;
  using CurrentCtx = CurrentGen::Context;

  auto resetter = to_export_context();

  REQUIRE(Index(L"i_1") < Index(L"a_1"));

  // Safe-guard that template magic works
  const std::size_t n_generators = 7;

  const std::set<std::string> known_formats =
      known_format_names(KnownGenerators{});
  REQUIRE(known_formats.size() == n_generators);

  const std::vector<std::filesystem::path> test_files =
      enumerate_export_tests();
  REQUIRE(!test_files.empty());

  for (const std::filesystem::path &current : test_files) {
    const std::string section_name =
        current.filename().string() + " - " + CurrentGen{}.get_format_name();

    SECTION(section_name) {
      CurrentGen generator;
      CurrentCtx context;

      configure_context_defaults(context);

      // Parse test file
      std::optional<std::string> expected_output;
      std::string expression_spec;
      {
        std::ifstream in(current.native());
        bool finished_expr = false;
        bool inside_meta = false;
        bool set_format = false;
        std::string current_format;
        for (std::string line; std::getline(in, line);) {
          if (line.starts_with("=====")) {
            finished_expr = true;
            set_format = false;
            inside_meta = !inside_meta;
          } else if (!finished_expr) {
            expression_spec += line;
            expression_spec += "\n";
          } else if (inside_meta) {
            if (line.starts_with("#")) {
              // Comment
              continue;
            }
            if (!set_format) {
              current_format = boost::trim_copy(line);
              if (known_formats.find(current_format) == known_formats.end()) {
                FAIL("Unknown format '" + current_format + "'");
              }
              set_format = true;
            } else if (current_format == generator.get_format_name()) {
              // Context definition
              add_to_context(context, line);
            }
          } else if (current_format == generator.get_format_name()) {
            if (expected_output.has_value()) {
              expected_output.value() += "\n";
              expected_output.value() += line;
            } else {
              expected_output = line;
            }
          }
        }
      }

      if (!expected_output.has_value()) {
        // This is not a test case for the current generator
        continue;
      }

      REQUIRE(!expression_spec.empty());

      auto groups = parse_expression_spec(
          expression_spec, std::same_as<CurrentGen, JuliaTensorKitGenerator<>>);
      REQUIRE(!groups.empty());

      export_groups<>(groups, generator, context);

      REQUIRE_THAT(generator.get_generated_code(),
                   DiffedStringEquals(expected_output.value()));
    }
  }
}

TEST_CASE("export", "[export]") {
  auto resetter = to_export_context();

  SECTION("reordering_context") {
    REQUIRE(Index(L"i_1").space().dimension() >
            Index(L"u_1").space().dimension());
    REQUIRE(Index(L"a_1").space().dimension() >
            Index(L"i_1").space().dimension());

    std::vector<std::pair<std::wstring, std::array<std::string, 3>>> tests = {
        // Unchanged
        {L"t{a1;i1}:N-N-N",
         {"t{a1;i1}:N-N-N", "t{a1;i1}:N-N-N", "t{a1;i1}:N-N-N"}},
        {L"t{a1,a2;i1,i2}:N-N-S",
         {"t{a1,a2;i1,i2}:N-N-S", "t{a1,a2;i1,i2}:N-N-S",
          "t{a1,a2;i1,i2}:N-N-S"}},
        // Bra resorting
        {L"t{a1,i1}:S-N-S",
         {"t{;;i1,a1}:N-N-N", "t{a1,i1}:S-N-S", "t{a1,i1}:S-N-S"}},
        {L"t{i1,a1}:S-N-S",
         {"t{i1,a1}:S-N-S", "t{;;a1,i1}:N-N-N", "t{i1,a1}:S-N-S"}},
        // Ket resorting
        {L"t{;a1,i1}:S-N-S",
         {"t{;;i1,a1}:N-N-N", "t{;a1,i1}:S-N-S", "t{;a1,i1}:S-N-S"}},
        {L"t{;i1,a1}:S-N-S",
         {"t{;i1,a1}:S-N-S", "t{;;a1,i1}:N-N-N", "t{;i1,a1}:S-N-S"}},
        {L"t{;u1,a1}:S-N-S",
         {"t{;u1,a1}:S-N-S", "t{;;a1,u1}:N-N-N", "t{;u1,a1}:S-N-S"}},
        // BraKet swapping
        {L"t{a1;i1}:N-S-N",
         {"t{;;i1,a1}:N-N-N", "t{a1;i1}:N-S-N", "t{a1;i1}:N-S-N"}},
        {L"t{i1;a1}:N-S-N",
         {"t{i1;a1}:N-S-N", "t{;;a1,i1}:N-N-N", "t{i1;a1}:N-S-N"}},
        // Aux prioritization
        {L"t{i1;;a1}:N-N-N",
         {"t{i1;;a1}:N-N-N", "t{;;a1,i1}:N-N-N", "t{i1;;a1}:N-N-N"}},
        {L"t{a1;;i1}:N-N-N",
         {"t{;;i1,a1}:N-N-N", "t{a1;;i1}:N-N-N", "t{a1;;i1}:N-N-N"}},
        // Column-resorting (particle-symmetry)
        {L"t{i1,i2;u1,a1}:N-N-S",
         {"t{i1,i2;u1,a1}:N-N-S", "t{;;i2,i1,a1,u1}:N-N-N",
          "t{i1,i2;u1,a1}:N-N-S"}},
        {L"t{u1,a1;i1,i2}:N-N-S",
         {"t{u1,a1;i1,i2}:N-N-S", "t{;;a1,u1,i2,i1}:N-N-N",
          "t{u1,a1;i1,i2}:N-N-S"}},
        {L"t{u1,a1;i1,u2}:N-N-S",
         {"t{u1,a1;i1,u2}:N-N-S", "t{;;a1,u1,u2,i1}:N-N-N",
          "t{u1,a1;i1,u2}:N-N-S"}},
        {L"t{i1,u2;u1,a1}:N-N-S",
         {"t{i1,u2;u1,a1}:N-N-S", "t{;;u2,i1,a1,u1}:N-N-N",
          "t{i1,u2;u1,a1}:N-N-S"}},
        // combined
        {L"t{u1,a1;u2,i1}:N-S-S",
         {"t{;;u2,i1,u1,a1}:N-N-N", "t{;;a1,u1,i1,u2}:N-N-N",
          "t{u1,a1;u2,i1}:N-S-S"}},
    };

    ReorderingContext ctx(MemoryLayout::Unspecified);

    for (MemoryLayout layout :
         {MemoryLayout::RowMajor, MemoryLayout::ColumnMajor,
          MemoryLayout::Unspecified}) {
      CAPTURE(layout);

      ctx.set_memory_layout(layout);

      for (const auto &[input, candidates] : tests) {
        CAPTURE(toUtf8(input));

        const std::string &expected =
            candidates.at(static_cast<std::size_t>(layout));

        Tensor tensor =
            deserialize(input, {.def_perm_symm = Symmetry::Nonsymm,
                                .def_braket_symm = BraKetSymmetry::Nonsymm,
                                .def_col_symm = ColumnSymmetry::Nonsymm})
                ->as<Tensor>();
        bool rewritten = ctx.rewrite(tensor);
        REQUIRE_THAT(tensor, EquivalentTo(expected));
        REQUIRE(rewritten == (toUtf8(input) != expected));
      }
    }
  }

  SECTION("reordering_context with densities") {
    // a density is reordered into the aux-only layout like any other tensor
    ReorderingContext ctx(MemoryLayout::ColumnMajor);
    for (const auto &[density, like] :
         {std::pair{L"γ{u1,a1;u2,i1}", L"t{u1,a1;u2,i1}:A-H-S"},
          std::pair{L"Γ{u1,a1;u2,i1}", L"t{u1,a1;u2,i1}:N-H-S"}}) {
      CAPTURE(toUtf8(density));
      Tensor tensor = deserialize(density)->as<Tensor>();
      Tensor reference = deserialize(like)->as<Tensor>();
      bool rewritten = false;
      REQUIRE_NOTHROW(rewritten = ctx.rewrite(tensor));
      REQUIRE(rewritten == ctx.rewrite(reference));
      REQUIRE(rewritten);
      REQUIRE(tensor.bra_rank() == 0);
      REQUIRE(tensor.ket_rank() == 0);
      REQUIRE(tensor.aux() == reference.aux());
    }

    // ... also when exporting
    JuliaTensorOperationsGenerator<> generator;
    JuliaTensorOperationsGeneratorContext julia_ctx;
    for (const auto *space : {"i", "a", "u"})
      julia_ctx.set_tag(
          get_default_context().index_basis_registry()->retrieve(space), space);
    REQUIRE_NOTHROW(export_expression(
        to_export_tree(deserialize<ResultExpr>(
            L"R{u1,a1;u2,i1} = γ{u1,a1;u2,i1} + Γ{u1,a1;u2,i1}")),
        generator, julia_ctx));
  }

  SECTION("generation_optimizer") {
    TextGenerator<TextGeneratorContext> textgen;
    GenerationOptimizer<TextGenerator<TextGeneratorContext>> generator(textgen);
    TextGeneratorContext ctx;

    Variable v1{L"v1"};
    Variable v2{L"v2"};
    Variable v3{L"v3"};

    SECTION("unchanged") {
      export_expression(to_export_tree(deserialize<ResultExpr>(L"v1 = 2 v2")),
                        generator, ctx);

      REQUIRE_THAT(generator.get_generated_code(),
                   DiffedStringEquals("Declare variable v1\n"
                                      "Declare variable v2\n"
                                      "\n"
                                      "Create v1 and initialize to zero\n"
                                      "Load v2\n"
                                      "Compute v1 += 2 v2\n"
                                      "Unload v2\n"
                                      "Persist v1\n"));
    }
    SECTION("elided load/unload") {
      SECTION("single") {
        export_expression(
            to_export_tree(deserialize<ResultExpr>(L"v1 = 2 v2 + 4 v2 v3")),
            generator, ctx);

        REQUIRE_THAT(generator.get_generated_code(),
                     DiffedStringEquals("Declare variable v1\n"
                                        "Declare variable v2\n"
                                        "Declare variable v3\n"
                                        "\n"
                                        "Create v1 and initialize to zero\n"
                                        "Load v2\n"
                                        "Load v3\n"
                                        "Compute v1 += 4 v2 v3\n"
                                        "Unload v3\n"
                                        "Compute v1 += 2 v2\n"
                                        "Unload v2\n"
                                        "Persist v1\n"));
      }
      SECTION("multiple") {
        export_expression(to_export_tree(deserialize<ResultExpr>(
                              L"ECC = 2 g{i1,i2;a1,a2} t{a1,a2;i1,i2} "
                              "- g{i1,i2;a1,a2} t{a2,a1;i1,i2}")),
                          generator, ctx);

        REQUIRE_THAT(
            generator.get_generated_code(),
            DiffedStringEquals(
                "Declare index i_1\n"
                "Declare index i_2\n"
                "Declare index a_1\n"
                "Declare index a_2\n"
                "\n"
                "Declare variable ECC\n"
                "\n"
                "Declare tensor g[i_1, i_2, a_1, a_2]\n"
                "Declare tensor t[a_1, a_2, i_1, i_2]\n"
                "\n"
                "Create ECC and initialize to zero\n"
                "Load g[i_1, i_2, a_1, a_2]\n"
                "Load t[a_1, a_2, i_1, i_2]\n"
                "Compute ECC += 2 g[i_1, i_2, a_1, a_2] t[a_1, a_2, i_1, i_2]\n"
                "Compute ECC += -1 g[i_1, i_2, a_1, a_2] t[a_2, a_1, i_1, "
                "i_2]\n"
                "Unload t[a_2, a_1, i_1, i_2]\n"
                "Unload g[i_1, i_2, a_1, a_2]\n"
                "Persist ECC\n"

                ));
      }
      SECTION("index batching") {
        ctx.set_batch_indices(std::vector<Index>{"i_1", "i_2"});

        export_expression(
            to_export_tree(deserialize<ResultExpr>(
                L"R{a1,a2;i1,i2} = A{a1,a2;i1,i2} + A{a3,a4;i1,i2} "
                L"B{a1,a2;a3,a4} + A{a1,a2;i3,i4} B{a3,a4;i1,i2}")),
            generator, ctx);
        REQUIRE_THAT(
            generator.get_generated_code(),
            DiffedStringEquals(
                "Declare index i_1\n"
                "Declare index i_2\n"
                "Declare index i_3\n"
                "Declare index i_4\n"
                "Declare index a_1\n"
                "Declare index a_2\n"
                "Declare index a_3\n"
                "Declare index a_4\n"
                "\n"
                "Declare tensor A[a_3, a_4, i_1, i_2]\n"
                "Declare tensor B[a_3, a_4, i_1, i_2]\n"
                "Declare tensor B[a_1, a_2, a_3, a_4]\n"
                "Declare tensor R[a_1, a_2, i_1, i_2]\n"
                "\n"
                "Start batching over i_1, i_2\n"
                "Create R[a_1, a_2, i_1, i_2] and initialize to zero\n"
                "Load A[a_3, a_4, i_1, i_2]\n"
                "Load B[a_1, a_2, a_3, a_4]\n"
                "Compute R[a_1, a_2, i_1, i_2] += A[a_3, a_4, i_1, i_2] B[a_1, "
                "a_2, a_3, a_4]\n"
                "Unload B[a_1, a_2, a_3, a_4]\n"
                "Compute R[a_1, a_2, i_1, i_2] += A[a_1, a_2, i_1, i_2]\n"
                // Note: this pair of unload/load A must not be eliminated due
                // to differences in their use of batching vs. non-batching
                // indices Hence, the former refers to only a multidimensional
                // slice of A, whereas the latter refers to the full A tensor.
                "Unload A[a_1, a_2, i_1, i_2]\n"
                "Load A[a_1, a_2, i_3, i_4]\n"
                "Load B[a_3, a_4, i_1, i_2]\n"
                "Compute R[a_1, a_2, i_1, i_2] += A[a_1, a_2, i_3, i_4] B[a_3, "
                "a_4, i_1, i_2]\n"
                "Unload B[a_3, a_4, i_1, i_2]\n"
                "Unload A[a_1, a_2, i_3, i_4]\n"
                "Persist R[a_1, a_2, i_1, i_2]\n"
                "End batching\n"));
      }

      // The following test cases will only become relevant once the optimizer
      // is further improved
#if 0
      SECTION("multiple with reordering") {
        // When stripping redundant load/unload operations, load of tensor C
        // needs to be moved before the load of B in order to retain
        // compatibility to frameworks with stack-based memory models
        export_expression(
            to_export_tree(deserialize<ResultExpr>(
                L"R{a1;i1} = 2 B{a1;i1} + B{a1;i1} C - C D{a1;i1}")),
            generator, ctx);

        REQUIRE_THAT(
            generator.get_generated_code(),
            DiffedStringEquals("Declare index i_1\n"
                               "Declare index a_1\n"
                               "\n"
                               "Declare variable C\n"
                               "\n"
                               "Declare tensor B[a_1, i_1]\n"
                               "Declare tensor D[a_1, i_1]\n"
                               "Declare tensor R[a_1, i_1]\n"
                               "\n"
                               "Create R[a_1, i_1] and initialize to zero\n"
                               "Load C\n"
                               "Load B[a_1, i_1]\n"
                               "Compute R[a_1, i_1] += 2 B[a_1, i_1]\n"
                               "Compute R[a_1, i_1] += B[a_1, i_1] C\n"
                               "Unload B[a_1, i_1]\n"
                               "Load D[a_1, i_1]\n"
                               "Compute R[a_1, i_1] += -1 C D[a_1, i_1]\n"
                               "Unload D[a_1, i_1]\n"
                               "Unload C\n"
                               "Persist R[a_1, i_1]\n"));
      }
      SECTION("reused intermediates") {
        export_expression(to_export_tree(deserialize<ResultExpr>(
                              L"ECC = 2 K{i1,i2;a1,a2} t{a1;i1} t{a2;i2} - "
                              L"K{i1,i2;a1,a2} t{a1;i2} t{a2;i1}")),
                          generator, ctx);

        REQUIRE_THAT(
            generator.get_generated_code(),
            DiffedStringEquals(
                "Declare index i_1\n"
                "Declare index i_2\n"
                "Declare index a_1\n"
                "Declare index a_2\n"
                "\n"
                "Declare variable ECC\n"
                "\n"
                "Declare tensor I[i_2, a_2]\n"
                "Declare tensor K[i_1, i_2, a_1, a_2]\n"
                "Declare tensor t[a_1, i_1]\n"
                "\n"
                "Create ECC and initialize to zero\n"
                "Load t[a_1, i_1]\n"
                "Create I[i_2, a_2] and initialize to zero\n"
                "Load K[i_1, i_2, a_1, a_2]\n"
                "Compute I[i_2, a_2] += K[i_1, i_2, a_1, a_2] t[a_1, i_1]\n"
                "Unload K[i_1, i_2, a_1, a_2]\n"
                "Compute ECC += 2 I[i_2, a_2] t[a_2, i_2]\n"
                "Unload I[i_2, a_2]\n"
                "Load I[i_1, a_2] and set it to zero\n"
                "Load K[i_1, i_2, a_1, a_2]\n"
                "Compute I[i_1, a_2] += K[i_1, i_2, a_1, a_2] t[a_1, i_2]\n"
                "Unload K[i_1, i_2, a_1, a_2]\n"
                "Compute ECC += -1 I[i_1, a_2] t[a_2, i_1]\n"
                "Unload I[i_1, a_2]\n"
                "Unload t[a_2, i_1]\n"
                "Persist ECC\n"));
      }
      SECTION("requires caution") {
        export_expression(
            to_export_tree(deserialize<ResultExpr>(
                L"R{a1,a2;i1,i2} = g{i3,i4;a3,a4} t{a4;i2} t{a1,a3;i1,i4} "
                L"t{a2;i3} "
                "- 2 g{i3,i4;a3,a4} t{a4;i2} t{a1,a3;i1,i3} t{a2;i4} ")),
            generator, ctx);

        REQUIRE_THAT(generator.get_generated_code(),
                     DiffedStringEquals(
                         "Declare index i_1\n"
                         "Declare index i_2\n"
                         "Declare index i_3\n"
                         "Declare index i_4\n"
                         "Declare index a_1\n"
                         "Declare index a_2\n"
                         "Declare index a_3\n"
                         "Declare index a_4\n"
                         "\n"
                         "Declare tensor I[i_3, i_4, i_2, a_3]\n"
                         "Declare tensor I[i_4, a_1, i_1, i_2]\n"
                         "Declare tensor R[a_1, a_2, i_1, i_2]\n"
                         "Declare tensor g[i_3, i_4, a_3, a_4]\n"
                         "Declare tensor t[a_4, i_2]\n"
                         "Declare tensor t[a_1, a_3, i_1, i_3]\n"
                         "\n"
                         "Create R[a_1, a_2, i_1, i_2] and initialize to zero\n"
                         "Create I[i_4, a_1, i_1, i_2] and initialize to zero\n"
                         "Create I[i_3, i_4, i_2, a_3] and initialize to zero\n"
                         "Load g[i_3, i_4, a_3, a_4]\n"
                         "Load t[a_4, i_2]\n"
                         "Compute I[i_3, i_4, i_2, a_3] += g[i_3, i_4, a_3, "
                         "a_4] t[a_4, i_2]\n"
                         "Unload t[a_4, i_2]\n"
                         "Unload g[i_3, i_4, a_3, a_4]\n"
                         "Load t[a_1, a_3, i_1, i_3]\n"
                         "Compute I[i_4, a_1, i_1, i_2] += I[i_3, i_4, i_2, "
                         "a_3] t[a_1, a_3, i_1, i_3]\n"
                         "Unload t[a_1, a_3, i_1, i_3]\n"
                         "Unload I[i_3, i_4, i_2, a_3]\n"
                         "Load t[a_2, i_4]\n"
                         "Compute R[a_1, a_2, i_1, i_2] += -2 I[i_4, a_1, i_1, "
                         "i_2] t[a_2, i_4]\n"
                         "Unload t[a_2, i_4]\n"
                         "Unload I[i_4, a_1, i_1, i_2]\n"
                         "Load I[i_3, a_1, i_1, i_2] and set it to zero\n"
                         "Load I[i_3, i_4, i_2, a_3] and set it to zero\n"
                         "Load g[i_3, i_4, a_3, a_4]\n"
                         "Load t[a_4, i_2]\n"
                         "Compute I[i_3, i_4, i_2, a_3] += g[i_3, i_4, a_3, "
                         "a_4] t[a_4, i_2]\n"
                         "Unload t[a_4, i_2]\n"
                         "Unload g[i_3, i_4, a_3, a_4]\n"
                         "Load t[a_1, a_3, i_1, i_4]\n"
                         "Compute I[i_3, a_1, i_1, i_2] += I[i_3, i_4, i_2, "
                         "a_3] t[a_1, a_3, i_1, i_4]\n"
                         "Unload t[a_1, a_3, i_1, i_4]\n"
                         "Unload I[i_3, i_4, i_2, a_3]\n"
                         "Load t[a_2, i_3]\n"
                         "Compute R[a_1, a_2, i_1, i_2] += I[i_3, a_1, i_1, "
                         "i_2] t[a_2, i_3]\n"
                         "Unload t[a_2, i_3]\n"
                         "Unload I[i_3, a_1, i_1, i_2]\n"
                         "Persist R[a_1, a_2, i_1, i_2]\n"));
      }
      SECTION("tbd2") {
        export_expression(
            to_export_tree(deserialize<ResultExpr>(
                L"R1{u_1;i_1;} = "
                L"+ 2 g{u_2, i_2;a_1, a_2;} (GAM0{u_3, u_4;u_5, u_2;} T2g{a_2, "
                L"u_5;u_3, i_2;}) T2g{a_1, u_1;i_1, u_4;} "
                "+ -4 g{u_2, i_2;a_1, a_2;} (GAM0{u_3, u_4;u_5, u_2;} T2g{a_1, "
                "u_5;i_2, u_3;}) T2g{a_2, u_1;u_4, i_1;} ")),
            generator, ctx);
        REQUIRE_THAT(generator.get_generated_code(), DiffedStringEquals(""));
      }
#endif
    }
  }

  SECTION("itf") {
    std::wstring int_label = L"g";
    ItfContext ctx;
    ctx.set_two_electron_integral_label(int_label);

    SECTION("remap_integrals") {
      SECTION("Unchanged") {
        Tensor tensor = deserialize(L"t{i1;a1}:N-N-N")->as<Tensor>();
        bool rewritten = ctx.rewrite(tensor);
        REQUIRE_THAT(tensor, EquivalentTo("t{i1;a1}:N-N-N"));
        REQUIRE_FALSE(rewritten);

        tensor = deserialize(int_label + L"{a1;i1}:N-N-N")->as<Tensor>();
        rewritten = ctx.rewrite(tensor);
        REQUIRE_THAT(tensor, EquivalentTo("g{a1;i1}:N-N-N"));
      }

      SECTION("K") {
        SECTION("occ,occ,occ,occ") {
          std::vector<Index> indices = {L"i_1", L"i_2", L"i_3", L"i_4"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            CAPTURE(indexPerm);

            REQUIRE(indexPerm.size() == 4);

            Tensor integral(int_label,
                            bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                            ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            bool rewritten = ctx.rewrite(integral);
            REQUIRE_THAT(integral, EquivalentTo("K{i1,i2;i3,i4}"));
            REQUIRE(rewritten);
          }
        }

        SECTION("virt,virt,occ,occ") {
          std::vector<Index> indices = {L"a_1", L"a_2", L"i_1", L"i_2"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            CAPTURE(indexPerm);
            REQUIRE(indexPerm.size() == 4);

            Tensor integral(int_label,
                            bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                            ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            bool rewritten = ctx.rewrite(integral);
            REQUIRE_THAT(integral, EquivalentTo("K{a1,a2;i1,i2}"));
            REQUIRE(rewritten);
          }
        }

        SECTION("virt,virt,virt,virt") {
          std::vector<Index> indices = {L"a_1", L"a_2", L"a_3", L"a_4"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            CAPTURE(indexPerm);
            REQUIRE(indexPerm.size() == 4);

            Tensor integral(int_label,
                            bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                            ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            bool rewritten = ctx.rewrite(integral);
            REQUIRE_THAT(integral, EquivalentTo("K{a1,a2;a3,a4}"));
            REQUIRE(rewritten);
          }
        }
      }

      SECTION("J") {
        SECTION("virt,occ,virt,occ") {
          std::vector<Index> indices = {L"a_1", L"i_1", L"a_2", L"i_2"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            CAPTURE(indexPerm);
            REQUIRE(indexPerm.size() == 4);

            Tensor integral(int_label,
                            bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                            ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            bool rewritten = ctx.rewrite(integral);
            REQUIRE_THAT(integral, EquivalentTo("J{a1,a2;i1,i2}"));
            REQUIRE(rewritten);
          }
        }

        SECTION("virt,occ,virt,virt") {
          std::vector<Index> indices = {L"a_1", L"i_1", L"a_2", L"a_3"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            CAPTURE(indexPerm);
            REQUIRE(indexPerm.size() == 4);

            Tensor integral(int_label,
                            bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                            ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            bool rewritten = ctx.rewrite(integral);
            REQUIRE_THAT(integral, EquivalentTo("J{a1,a2;a3,i1}"));
            REQUIRE(rewritten);
          }
        }
      }
    }
  }
}

TEST_CASE("ExportExpr", "[export]") {
  SECTION("id & equality") {
    Variable v(L"V");

    ExportExpr e1(v);
    ExportExpr e2(v);

    REQUIRE(e1.id() != e2.id());
    REQUIRE(e1 != e2);

    ExportExpr e3 = e1;
    REQUIRE(e1.id() == e3.id());
    REQUIRE(e1 == e3);
  }
}

TEST_CASE("JuliaTensorOperationsGenerator", "[export]") {
  SECTION("load scalar and set to zero") {
    JuliaTensorOperationsGenerator<> generator;
    JuliaTensorOperationsGeneratorContext ctx;
    generator.set_to_zero(Variable(L"x"), ctx);
    generator.load(Variable(L"y"), /*set_to_zero=*/true, ctx);

    REQUIRE(generator.get_generated_code() == "x = 0.0\ny = 0.0\n");
  }
}

TEST_CASE("PythonEinsumGenerator", "[export]") {
  auto resetter = to_export_context();

  auto registry = get_default_context().index_basis_registry();
  IndexSpace occ = registry->retrieve("i");
  IndexSpace virt = registry->retrieve("a");

  SECTION("NumPy backend - simple contraction") {
    // Create tensors: T[a1, a2] = F[a1, i1] * t[i1, a2]
    auto F = ex<Tensor>(L"F", bra{L"a_1"}, ket{L"i_1"});
    auto t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_2"});
    Tensor T(L"T", bra{L"a_1"}, ket{L"a_2"});

    // Build the expression
    ResultExpr result_expr(T, F * t);

    // Convert to export tree
    auto export_tree = to_export_tree(result_expr);

    // Set up context
    NumPyEinsumGeneratorContext ctx;
    ctx.set_shape(occ, "nocc");
    ctx.set_shape(virt, "nvirt");
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    // Generate code
    NumPyEinsumGenerator generator;
    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    // Verify generated code contains expected elements
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("np.zeros"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("np.load"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("np.einsum"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("optimize=True"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("T_vv +="));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("np.save"));
  }

  SECTION("PyTorch backend - simple contraction") {
    // Create tensors: T[a1, a2] = F[a1, i1] * t[i1, a2]
    auto F = ex<Tensor>(L"F", bra{L"a_1"}, ket{L"i_1"});
    auto t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_2"});
    Tensor T(L"T", bra{L"a_1"}, ket{L"a_2"});

    ResultExpr result_expr(T, F * t);
    auto export_tree = to_export_tree(result_expr);

    // Set up PyTorch context
    PyTorchEinsumGeneratorContext ctx;
    ctx.set_shape(occ, "nocc");
    ctx.set_shape(virt, "nvirt");
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    PyTorchEinsumGenerator generator;
    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    // Verify PyTorch-specific code
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("torch.zeros"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("torch.load"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("torch.einsum"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("torch.save"));
    REQUIRE_THAT(code, !Catch::Matchers::ContainsSubstring("optimize=True"));
  }

  SECTION("Scalar factor") {
    // Test with scalar prefactor: R = 0.5 * t1 * t2
    Index i1(L"i_1", occ);
    Index i2(L"i_2", occ);
    Index a1(L"a_1", virt);
    Index a2(L"a_2", virt);

    auto t1 = ex<Tensor>(L"t1", bra{L"i_1"}, ket{L"a_1"});
    auto t2 = ex<Tensor>(L"t2", bra{L"i_2"}, ket{L"a_2"});
    Tensor R(L"R", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});

    ResultExpr result_with_scalar(R, rational(1, 2) * t1 * t2);
    auto export_tree_scalar = to_export_tree(result_with_scalar);

    PyTorchEinsumGeneratorContext ctx;
    ctx.set_shape(occ, "nocc");
    ctx.set_shape(virt, "nvirt");
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    PyTorchEinsumGenerator generator;
    export_expression(export_tree_scalar, generator, ctx);

    std::string code = generator.get_generated_code();

    // Verify scalar factor is included
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("1/2"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("*"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring(".einsum('"));
  }

  SECTION("Scalar result - energy expression") {
    auto g = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});
    auto t = ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"});
    Variable E(L"E");

    ResultExpr result_expr(E, g * t);
    auto export_tree = to_export_tree(result_expr);

    PyTorchEinsumGeneratorContext ctx;
    ctx.set_shape(occ, "nocc");
    ctx.set_shape(virt, "nvirt");
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");

    PyTorchEinsumGenerator generator;
    export_expression(export_tree, generator, ctx);

    std::string code = generator.get_generated_code();

    // Verify scalar result (einsum ending with "->")
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("E +="));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("->'"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring(".einsum('"));
  }
}

TEST_CASE("export-basis-instance", "[export][basis]") {
  using Catch::Matchers::ContainsSubstring;
  auto resetter = to_export_context();
  const Index y = Index(L"a_7").replace_basis_instance(1);

  auto registry = get_default_context().index_basis_registry();
  const IndexSpace occ = registry->retrieve("i");
  const IndexSpace virt = registry->retrieve("a");
  const Tensor t(L"t", bra{y}, ket{L"i_1"});
  const Tensor t_null(L"t", bra{L"a_7"}, ket{L"i_1"});

  // a mode in a basis instance appends the instance to its block tag; the
  // generic block is named as before
  SECTION("ITF") {
    ItfContext ctx;
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");
    const ItfGenerator<ItfContext> generator;
    CHECK_THAT(generator.represent(t, ctx), ContainsSubstring("t:v1o["));
    CHECK_THAT(generator.represent(t_null, ctx), ContainsSubstring("t:vo["));
  }
  SECTION("Julia") {
    JuliaTensorOperationsGeneratorContext ctx;
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");
    const JuliaTensorOperationsGenerator<> generator;
    CHECK_THAT(generator.represent(t, ctx),
               ContainsSubstring("t_v1o[") && ContainsSubstring("a_7_1"));
    CHECK_THAT(generator.represent(t_null, ctx),
               ContainsSubstring("t_vo[") && ContainsSubstring("a_7"));
  }
  SECTION("Julia TensorKit") {
    // the domain shows the extent the tensor is allocated with
    JuliaTensorKitGeneratorContext ctx;
    ctx.set_tag(occ, "o");
    ctx.set_tag(virt, "v");
    JuliaTensorKitGenerator<> generator;
    generator.create(t, true, ctx);
    const auto code = generator.get_generated_code();
    const auto nv1 = ctx.get_dim(virt) + "_1";
    CHECK_THAT(code, ContainsSubstring("zeros(Float64, " + nv1) &&
                         ContainsSubstring("ℝ^" + nv1));
  }
  SECTION("text") {
    const TextGenerator<TextGeneratorContext> generator;
    CHECK(generator.represent(t, TextGeneratorContext{}) == "t[a_7<;1>, i_1]");
    CHECK(generator.represent(t_null, TextGeneratorContext{}) == "t[a_7, i_1]");
  }
  // the einsum exporters name blocks and shapes per mode, too
  SECTION("Python einsum") {
    auto exported = [&](auto generator, auto ctx) {
      ctx.set_shape(occ, "nocc");
      ctx.set_shape(virt, "nvirt");
      ctx.set_tag(occ, "o");
      ctx.set_tag(virt, "v");
      // the result is in basis 1, so its allocation shows the basis' extent
      ResultExpr result(Tensor(L"R", bra{y}, ket{L"i_1"}),
                        ex<Tensor>(L"f", bra{y}, ket{L"a_1"}) *
                            ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}));
      export_expression(to_export_tree(result), generator, ctx);
      return generator.get_generated_code();
    };
    for (auto const &code :
         {exported(NumPyEinsumGenerator{}, NumPyEinsumGeneratorContext{}),
          exported(PyTorchEinsumGenerator{}, PyTorchEinsumGeneratorContext{})})
      CHECK_THAT(code, ContainsSubstring("f_v1v") &&
                           ContainsSubstring("t_vo") &&
                           ContainsSubstring("R_v1o") &&
                           ContainsSubstring("(nvirt_1, nocc)"));
  }
  // a tag ending in a digit or `_`, or a dimension name ending in `_`, would
  // make an appended instance ambiguous (v1 + "" vs v + 1, nvirt_ + _1 vs
  // nvirt + _ + 1), so an instance-bearing mode in such a space is refused
  SECTION("ambiguous space tag") {
    const Index z = Index(L"a_7").replace_basis_instance(1);
    const Tensor tz(L"t", bra{z}, ket{L"i_1"});
    const Tensor t0(L"t", bra{L"a_7"}, ket{L"i_1"});
    ItfContext itf_ctx;
    itf_ctx.set_tag(occ, "o");
    itf_ctx.set_tag(virt, "v1");
    const ItfGenerator<ItfContext> itf;
    CHECK_NOTHROW(itf.represent(t0, itf_ctx));
    CHECK_THROWS_AS(itf.represent(tz, itf_ctx), Exception);
    itf_ctx.set_tag(virt, "v_");
    CHECK_THROWS_AS(itf.represent(tz, itf_ctx), Exception);
    JuliaTensorOperationsGeneratorContext julia_ctx;
    julia_ctx.set_tag(occ, "o");
    julia_ctx.set_tag(virt, "v");
    julia_ctx.set_dim(occ, "nocc");
    julia_ctx.set_dim(virt, "nvirt_");
    JuliaTensorOperationsGenerator<JuliaTensorOperationsGeneratorContext> julia;
    CHECK_NOTHROW(julia.create(t0, true, julia_ctx));
    CHECK_THROWS_AS(julia.create(tz, true, julia_ctx), Exception);
    NumPyEinsumGeneratorContext py_ctx;
    py_ctx.set_shape(occ, "nocc");
    py_ctx.set_shape(virt, "nvirt_");
    CHECK_NOTHROW(py_ctx.get_shape_tuple(t0));
    CHECK_THROWS_AS(py_ctx.get_shape_tuple(tz), Exception);
  }

  // a negative instance is spelled with an `_` in place of the minus sign,
  // which no identifier may contain; a letter would read as the tag of a space
  SECTION("negative instance") {
    const Index z = Index(L"a_7").replace_basis_instance(-12);
    const Tensor tz(L"t", bra{z}, ket{L"i_1"});
    ItfContext itf_ctx;
    itf_ctx.set_tag(occ, "o");
    itf_ctx.set_tag(virt, "v");
    CHECK_THAT(ItfGenerator<ItfContext>{}.represent(tz, itf_ctx),
               ContainsSubstring("t:v_12o["));
    // with a space tagged `m`, a letter for the sign would make t{a<;-1>;i}
    // and t{a,u<;1>;i} one block
    const IndexSpace u = registry->retrieve("u");
    itf_ctx.set_tag(u, "m");
    const Tensor t_neg(L"t", bra{Index(L"a_1").replace_basis_instance(-1)},
                       ket{L"i_1"});
    const Tensor t_m(
        L"t", bra{Index(L"a_1"), Index(L"u_1").replace_basis_instance(1)},
        ket{L"i_1"});
    const ItfGenerator<ItfContext> itf;
    CHECK(itf.get_name(t_neg, itf_ctx) != itf.get_name(t_m, itf_ctx));
    JuliaTensorOperationsGeneratorContext julia_ctx;
    julia_ctx.set_tag(occ, "o");
    julia_ctx.set_tag(virt, "v");
    CHECK_THAT(JuliaTensorOperationsGenerator<>{}.represent(tz, julia_ctx),
               ContainsSubstring("t_v_12o[") && ContainsSubstring("a_7__12"));
    NumPyEinsumGeneratorContext py_ctx;
    py_ctx.set_shape(occ, "nocc");
    py_ctx.set_shape(virt, "nvirt");
    py_ctx.set_tag(occ, "o");
    py_ctx.set_tag(virt, "v");
    NumPyEinsumGenerator py;
    ResultExpr result(Tensor(L"R", bra{z}, ket{L"i_1"}),
                      ex<Tensor>(L"f", bra{z}, ket{L"a_1"}) *
                          ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}));
    export_expression(to_export_tree(result), py, py_ctx);
    CHECK_THAT(py.get_generated_code(),
               ContainsSubstring("R_v_12o") && ContainsSubstring("f_v_12v"));
    CHECK_THAT(py.get_generated_code(), !ContainsSubstring("-12"));
    CHECK(py_ctx.get_shape_tuple(tz) == "(nvirt__12, nocc)");
  }
}
