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
#include <SeQuant/core/index_space_registry.hpp>
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
  auto registry = get_default_context().index_space_registry();
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
  auto registry = get_default_context().index_space_registry();
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
  auto registry = get_default_context().index_space_registry();
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
  auto registry = get_default_context().index_space_registry();
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
        get_default_context().index_space_registry()->retrieve(space),
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
  // the generated code is real arithmetic: the test tensors live over a real
  // basis, where the `S` braket letter is derivable
  auto real_basis = tests::scoped_real_basis();

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
    // the `S` braket letters below are derivable only over a real basis
    auto real_basis = tests::scoped_real_basis();
    REQUIRE(Index(L"i_1").space().approximate_size() >
            Index(L"u_1").space().approximate_size());
    REQUIRE(Index(L"a_1").space().approximate_size() >
            Index(L"i_1").space().approximate_size());

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
                              "- g{i1,i2;a1,a2} t{a2,a1;i1,i2}",
                              {.def_braket_symm = Hermiticity::NonHermitian})),
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

TEST_CASE("PythonEinsumGenerator", "[export]") {
  auto resetter = to_export_context();

  auto registry = get_default_context().index_space_registry();
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

TEST_CASE("exported names of marked tensors", "[export]") {
  using namespace sequant;
  auto resetter = to_export_context();
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
           TensorSymmetries{.conjugation_parity = ConjugationParity::None});
  Tensor ta = t;
  REQUIRE(ta.adjoint() == 1);
  Tensor tk = t;
  REQUIRE(tk.kconjugate() == 1);
  TextGenerator<TextGeneratorContext> gen;
  TextGeneratorContext ctx;
  REQUIRE(gen.represent(t, ctx) == "t[a_1, i_1]");
  REQUIRE(gen.represent(ta, ctx) == "t_adj[i_1, a_1]");
  REQUIRE(gen.represent(tk, ctx) == "t_conj[a_1, i_1]");
}

TEST_CASE("exported names of reordered marked tensors", "[export]") {
  using namespace sequant;
  auto resetter = to_export_context();

  ItfContext ctx;
  configure_context_defaults(ctx);
  ItfGenerator<ItfContext> gen;

  // the array name that the generator produces, without the index-space tags
  auto array_name = [&](const Tensor &tensor) {
    const std::string name = gen.get_name(tensor, ctx);
    return name.substr(0, name.find(':'));
  };

  const Tensor t(
      L"t", bra{L"i_1", L"a_1"}, ket{L"i_2", L"a_2"},
      TensorSymmetries{.perm = Symmetry::Symm,
                       .conjugation_parity = ConjugationParity::None});

  Tensor adjointed = t;
  REQUIRE(adjointed.adjoint() == 1);
  REQUIRE(adjointed.adjointed());
  REQUIRE(array_name(adjointed) == "t_adj");

  Tensor kconjugated = t;
  REQUIRE(kconjugated.kconjugate() == 1);
  REQUIRE(kconjugated.kconjugated());
  REQUIRE(array_name(kconjugated) == "t_conj");

  // the index reordering must leave the array name of a marked tensor alone
  Tensor reordered = t;
  REQUIRE(ctx.rewrite(reordered));
  REQUIRE(array_name(reordered) == "t");

  Tensor reordered_adjointed = adjointed;
  REQUIRE(ctx.rewrite(reordered_adjointed));
  REQUIRE(array_name(reordered_adjointed) == "t_adj");

  Tensor reordered_kconjugated = kconjugated;
  REQUIRE(ctx.rewrite(reordered_kconjugated));
  REQUIRE(array_name(reordered_kconjugated) == "t_conj");

  // the import-name map tells the marked arrays apart from the bare one
  ctx.set_import_name(reordered, "T");
  ctx.set_import_name(reordered_adjointed, "TADJ");
  ctx.set_import_name(reordered_kconjugated, "TCONJ");
  REQUIRE(ctx.import_name(reordered).value() == "T");
  REQUIRE(ctx.import_name(reordered_adjointed).value() == "TADJ");
  REQUIRE(ctx.import_name(reordered_kconjugated).value() == "TCONJ");
}

TEST_CASE("full export of marked tensors", "[export]") {
  using namespace sequant;
  auto resetter = to_export_context();

  const auto t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_2"});
  const auto t_adj = ex<Tensor>(L"t⁺", bra{L"i_1"}, ket{L"a_2"});
  REQUIRE(t_adj->as<Tensor>().adjointed());
  const auto f = ex<Tensor>(L"f", bra{L"a_2"}, ket{L"a_1"});
  const Tensor R(L"R", bra{L"i_1"}, ket{L"a_1"});
  const ResultExpr result(R, t * f + t_adj * f);

  SECTION("itf tells a marked tensor apart from its bare twin") {
    // the reordering leaves these Nonsymm, aux-less tensors alone, so only the
    // marks can keep the two arrays apart
    for (bool rewriting : {true, false}) {
      CAPTURE(rewriting);

      ItfContext ctx;
      configure_context_defaults(ctx);
      ctx.enable_rewriting(rewriting);
      ItfGenerator<ItfContext> gen;
      export_expression(to_export_tree(result), gen, ctx);
      const std::string code = gen.get_generated_code();
      CAPTURE(code);

      // two declarations under two names
      REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("tensor: t:ce["));
      REQUIRE_THAT(code,
                   Catch::Matchers::ContainsSubstring("tensor: t_adj:ce["));
      // and two terminals in the load-strategy bookkeeping: the marked array
      // is one the host code supplies, not one the generated code builds
      REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("load t:ce["));
      REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("load t_adj:ce["));
    }
  }

  SECTION("the text generator tells them apart as well") {
    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    export_expression(to_export_tree(result), gen, ctx);
    const std::string code = gen.get_generated_code();
    CAPTURE(code);

    REQUIRE_THAT(
        code, Catch::Matchers::ContainsSubstring("Declare tensor t[i_1, a_2]"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring(
                           "Declare tensor t_adj[i_1, a_2]"));
    // the marked array is contracted under its own name, in the slot order
    // the `⁺` denotes
    REQUIRE_THAT(code,
                 Catch::Matchers::ContainsSubstring(
                     "Compute R[i_1, a_1] += t_adj[i_1, a_2] f[a_2, a_1]"));
  }

  SECTION("an import name set on the marked tensor as written is honoured") {
    // the name is registered on the tensor that still carries its `⁺`, before
    // any rewrite; the map is keyed on the folded label, which is the name the
    // rest of the pipeline looks up
    ItfContext ctx;
    configure_context_defaults(ctx);
    ctx.set_import_name(t_adj->as<Tensor>(), "TADJ");
    REQUIRE(ctx.import_name(t_adj->as<Tensor>()).value() == "TADJ");

    ItfGenerator<ItfContext> gen;
    export_expression(to_export_tree(ResultExpr(R, t_adj * f)), gen, ctx);
    const std::string code = gen.get_generated_code();
    CAPTURE(code);

    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("TADJ"));
  }

  SECTION("a K-conjugated terminal is imported under its own name") {
    // over a complex basis a `꙳` is a leaf of its own, so both arrays are
    // terminals and both are imported
    const TensorSymmetries syms{.conjugation_parity = ConjugationParity::None};
    const auto tp = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_2"}, syms);
    const auto tc = ex<Tensor>(L"t꙳", bra{L"i_1"}, ket{L"a_2"}, syms);
    REQUIRE(tc->as<Tensor>().kconjugated());

    ItfContext ctx;
    configure_context_defaults(ctx);
    ItfGenerator<ItfContext> gen;
    export_expression(to_export_tree(ResultExpr(R, tp * f + tc * f)), gen, ctx);
    const std::string code = gen.get_generated_code();
    CAPTURE(code);

    REQUIRE_THAT(code,
                 Catch::Matchers::ContainsSubstring("tensor: t:ce[jc], t:ce"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring(
                           "tensor: t_conj:ce[jc], t_conj:ce"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("load t:ce[jc]"));
    REQUIRE_THAT(code,
                 Catch::Matchers::ContainsSubstring("load t_conj:ce[jc]"));
  }

  SECTION("the integral remap does not see a marked integral") {
    // the marks are folded into the label before any context rewrite, so a
    // marked integral reaches the g->J/K remap under its own name and is not
    // recognized by it
    const auto g = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});
    const auto g_adj =
        ex<Tensor>(L"g⁺", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});
    REQUIRE(g_adj->as<Tensor>().adjointed());
    const auto t2 = ex<Tensor>(L"t2", bra{L"a_1", L"a_2"}, ket{L"i_2", L"a_3"});
    const Tensor R2(L"R2", bra{L"i_1"}, ket{L"a_3"});

    ItfContext ctx;
    configure_context_defaults(ctx);
    ctx.set_two_electron_integral_label(L"g");
    ItfGenerator<ItfContext> gen;
    export_expression(to_export_tree(ResultExpr(R2, g * t2 + g_adj * t2)), gen,
                      ctx);
    const std::string code = gen.get_generated_code();
    CAPTURE(code);

    // the bare integral is still remapped
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("tensor: K:eecc["));
    // the marked one is not
    REQUIRE_THAT(code,
                 Catch::Matchers::ContainsSubstring("tensor: g_adj:ccee["));
    REQUIRE_THAT(code, !Catch::Matchers::ContainsSubstring("K_adj"));
    REQUIRE_THAT(code, !Catch::Matchers::ContainsSubstring("J_adj"));
  }
}

TEST_CASE("a marked leaf reaches the generator as the array it denotes",
          "[export]") {
  using namespace sequant;
  auto resetter = to_export_context();

  auto generate = [](const ResultExpr &result) {
    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    export_expression(to_export_tree(result), gen, ctx);
    return gen.get_generated_code();
  };

  const auto f = ex<Tensor>(L"f", bra{L"a_1"}, ket{L"a_2"});
  const Tensor R(L"R", bra{L"i_1"}, ket{L"a_2"});

  SECTION("an adjointed leaf is named _adj") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
    REQUIRE(t.adjoint() == 1);
    REQUIRE(t.adjointed());

    const std::string code = generate(ResultExpr(R, ex<Tensor>(t) * f));
    CAPTURE(code);

    // the `⁺` is spelled by the array name, over the slots it denotes
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("t_adj[i_1, a_1]"));
    // and the array the leaf stores is nowhere in the generated code
    REQUIRE_THAT(code, !Catch::Matchers::ContainsSubstring(" t["));
  }

  SECTION("a K-conjugated leaf is named _conj") {
    // over this complex basis a `꙳` whose parity leaves it unresolved names
    // an array of its own, which the leaf stores as written
    const TensorSymmetries no_parity{.conjugation_parity =
                                         ConjugationParity::None};
    Tensor r(L"r", bra{L"a_1"}, ket{L"i_1"}, no_parity);
    REQUIRE(r.kconjugate() == 1);
    REQUIRE(r.kconjugated());

    // the dummy pairs r's bra with a ket, as the network requires
    const auto w = ex<Tensor>(L"w", bra{L"a_2"}, ket{L"a_1"});
    const std::string code = generate(ResultExpr(R, ex<Tensor>(r) * w));
    CAPTURE(code);

    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("r_conj[a_1, i_1]"));
    REQUIRE_THAT(code, !Catch::Matchers::ContainsSubstring(" r["));
  }
}

TEST_CASE("a pruned scalar prefactor keeps its conjugation", "[export]") {
  using namespace sequant;
  auto resetter = to_export_context();

  // the text generator prunes every scalar it can
  REQUIRE(TextGenerator<TextGeneratorContext>{}.prunable_scalars() ==
          PrunableScalars::All);

  auto generate = [](const ResultExpr &result) {
    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    export_expression(to_export_tree(result), gen, ctx);
    return gen.get_generated_code();
  };

  const auto t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"});
  const auto f = ex<Tensor>(L"f", bra{L"a_1"}, ket{L"a_2"});

  // a scalar leaf stores the unmarked spelling, its conjugation riding the
  // node's transform
  auto conjugated_power = []() {
    auto p = ex<Power>(L"x", rational(2));
    p->as<Power>().conjugate();
    return p;
  };

  SECTION("pruned out of the tree") {
    // two tensor factors leave the product a subtree to hold the pruned
    // scalar, so the prefactor is taken out of the tree
    const Tensor R(L"R", bra{L"i_1"}, ket{L"a_2"});
    const std::string code =
        generate(ResultExpr(R, conjugated_power() * t * f));
    CAPTURE(code);

    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("conj(x^2)"));
  }

  SECTION("kept in the tree") {
    // one tensor factor: pruning the scalar would make the tree vanish, so it
    // reaches the generator through the computation instead
    const Tensor R(L"R", bra{L"i_1"}, ket{L"a_1"});
    const std::string code = generate(ResultExpr(R, conjugated_power() * t));
    CAPTURE(code);

    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("conj(x^2)"));
  }

  SECTION("a conjugated power base is spelled as written") {
    // the base's own marker is part of the stored spelling, not of the
    // transform, and survives the pruning unchanged
    auto x = ex<Variable>(L"x");
    x->as<Variable>().conjugate();
    const auto p = ex<Power>(std::move(x), rational(2));

    const Tensor R(L"R", bra{L"i_1"}, ket{L"a_2"});
    const std::string code = generate(ResultExpr(R, p * t * f));
    CAPTURE(code);

    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("conj(x)^2"));
  }
}

TEST_CASE("a folded name must come from one tensor", "[export]") {
  using namespace sequant;
  auto resetter = to_export_context();

  // `t⁺` is exported as `t_adj`, so a tensor already written `t_adj` would
  // share every label-keyed map with it and the generated code would read one
  // buffer under two readings
  const auto t_adj_written = ex<Tensor>(L"t_adj", bra{L"i_1"}, ket{L"a_2"});
  REQUIRE_FALSE(t_adj_written->as<Tensor>().adjointed());
  const auto t_marked = ex<Tensor>(L"t⁺", bra{L"i_1"}, ket{L"a_2"});
  REQUIRE(t_marked->as<Tensor>().adjointed());
  REQUIRE(export_label(t_marked->as<Tensor>()) ==
          export_label(t_adj_written->as<Tensor>()));
  const auto f = ex<Tensor>(L"f", bra{L"a_2"}, ket{L"a_1"});
  const Tensor R(L"R", bra{L"i_1"}, ket{L"a_1"});

  SECTION("the collision is refused") {
    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    REQUIRE_THROWS_AS(
        export_expression(
            to_export_tree(ResultExpr(R, t_adj_written * f + t_marked * f)),
            gen, ctx),
        Exception);
  }

  SECTION("distinct folded names are fine") {
    const auto t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_2"});
    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    REQUIRE_NOTHROW(export_expression(
        to_export_tree(ResultExpr(R, t * f + t_marked * f)), gen, ctx));
    const std::string code = gen.get_generated_code();
    CAPTURE(code);
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring(
                           "Declare tensor t_adj[i_1, a_2]"));
  }

  SECTION("one tensor used twice is not a collision") {
    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    REQUIRE_NOTHROW(export_expression(
        to_export_tree(ResultExpr(R, t_marked * f + t_marked * f)), gen, ctx));
  }

  SECTION("the guard spans every tree of one export") {
    // the declarations of all trees are merged into one global block, so two
    // trees are as much a collision as two factors of one tree
    const Tensor S(L"S", bra{L"i_1"}, ket{L"a_1"});
    std::vector<ExpressionGroup<>> groups;
    groups.emplace_back();
    groups.back().add(to_export_tree(ResultExpr(R, t_adj_written * f)));
    groups.back().add(to_export_tree(ResultExpr(S, t_marked * f)));

    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    REQUIRE_THROWS_AS(export_groups<>(std::move(groups), gen, ctx), Exception);
  }
}

TEST_CASE("a context rewrite folds a tensor's marks into its label",
          "[export]") {
  using namespace sequant;
  auto resetter = to_export_context();

  ItfContext ctx;
  configure_context_defaults(ctx);
  ctx.set_two_electron_integral_label(L"g");

  // a plain Nonsymm, aux-less tensor gives the reordering nothing to do, so
  // only the fold can report a change
  Tensor bare(L"t", bra{L"i_1"}, ket{L"a_1"});
  REQUIRE_FALSE(ctx.rewrite(bare));
  REQUIRE(bare.label() == L"t");

  Tensor marked(L"t⁺", bra{L"i_1"}, ket{L"a_1"});
  REQUIRE(marked.adjointed());
  REQUIRE(ctx.rewrite(marked));
  REQUIRE(marked.label() == L"t_adj");
  REQUIRE_FALSE(marked.adjointed());
  REQUIRE_FALSE(marked.kconjugated());
  REQUIRE_THAT(marked, EquivalentTo("t_adj{i_1;a_1}"));

  // the fold precedes the two-electron integral remap, which matches on the
  // bare label
  Tensor bare_integral(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});
  REQUIRE(ctx.rewrite(bare_integral));
  REQUIRE((bare_integral.label() == L"J" || bare_integral.label() == L"K"));

  Tensor marked_integral(L"g⁺", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});
  REQUIRE(marked_integral.adjointed());
  REQUIRE(ctx.rewrite(marked_integral));
  REQUIRE(marked_integral.label() == L"g_adj");
  REQUIRE_FALSE(marked_integral.adjointed());
}
