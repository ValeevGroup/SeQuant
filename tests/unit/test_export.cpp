#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/export_expr.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/export/expression_group.hpp>
#include <SeQuant/core/export/generation_optimizer.hpp>
#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/export/julia_itensor.hpp>
#include <SeQuant/core/export/julia_tensor_kit.hpp>
#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/export/reordering_context.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include "catch2_sequant.hpp"

#include <boost/algorithm/string.hpp>

#include <filesystem>
#include <fstream>
#include <optional>
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
    throw std::runtime_error("Invalid base dir for export tests");
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
    ItfGenerator<ItfContext>
>;
// clang-format on

template <typename Generator>
std::string get_format_name() {
  Generator g;

  return g.get_format_name();
}

template <typename... Generator>
std::set<std::string> known_format_names(std::tuple<Generator...> generators) {
  std::set<std::string> names;

  (names.insert(get_format_name<Generator>()), ...);

  return names;
}

void configure_context_defaults(TextGeneratorContext &ctx) {}

void configure_context_defaults(ItfContext &ctx) {
  auto registry = get_default_context().index_space_registry();
  IndexSpace occ = registry->retrieve("i");
  IndexSpace virt = registry->retrieve("a");

  ctx.set_tag(occ, "c");
  ctx.set_tag(virt, "e");

  ctx.set_name(occ, "Closed");
  ctx.set_name(virt, "External");
}

void configure_context_defaults(JuliaTensorOperationsGeneratorContext &ctx) {
  auto registry = get_default_context().index_space_registry();
  IndexSpace occ = registry->retrieve("i");
  IndexSpace virt = registry->retrieve("a");

  ctx.set_dim(occ, "nocc");
  ctx.set_dim(virt, "nv");

  ctx.set_tag(occ, "o");
  ctx.set_tag(virt, "v");
}

void add_to_context(TextGeneratorContext &ctx, std::string_view key,
                    std::string_view value) {
  throw std::runtime_error(
      "TextGeneratorContext doesn't support specifications");
}
void add_to_context(ItfContext &ctx, std::string_view key,
                    std::string_view value) {}
void add_to_context(JuliaTensorOperationsGeneratorContext &ctx,
                    std::string_view key, std::string_view value) {
  auto parse_space_map = [](std::string_view spec) {
    auto pos = spec.find("->");
    if (pos == std::string_view::npos) {
      throw std::runtime_error("Malformed space map");
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
    throw std::runtime_error(
        "Unsupported key in Julia context specification '" + std::string(key) +
        "'");
  }
}

template <typename Context>
void add_to_context(Context &ctx, const std::string &line) {
  auto pos = line.find(":");
  if (pos == std::string::npos) {
    throw std::runtime_error(
        "Malformed context specification: missing ':' in '" + line + "'");
  }

  std::string key = line.substr(0, pos);
  boost::trim(key);
  std::string value = line.substr(pos + 1);
  boost::trim(value);

  if (key.empty()) {
    throw std::runtime_error("Malformed context specification: Empty key");
  }
  if (value.empty()) {
    throw std::runtime_error("Malformed context specification: Empty value");
  }

  add_to_context(ctx, key, value);
}

std::vector<ExpressionGroup<>> parse_expression_spec(const std::string &spec) {
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
      assert(!name.empty());
      groups.emplace_back(std::move(name));
      continue;
    }

    if (groups.empty()) {
      groups.emplace_back();
    }

    try {
      ResultExpr res =
          parse_result_expr(to_wstring(line), Symmetry::nonsymm,
                            BraKetSymmetry::nonsymm, ParticleSymmetry::nonsymm);
      groups.back().add(to_export_tree(res));
    } catch (...) {
      ExprPtr expr =
          parse_expr(to_wstring(line), Symmetry::nonsymm,
                     BraKetSymmetry::nonsymm, ParticleSymmetry::nonsymm);
      groups.back().add(to_export_tree(expr));
    }
  }

  return groups;
}

[[nodiscard]] auto to_export_context() {
  auto reg = std::make_shared<IndexSpaceRegistry>();
  reg->add(L"i", 0b001, is_particle, 10);
  reg->add(L"a", 0b010, is_vacuum_occupied, is_reference_occupied, is_hole,
           100);
  reg->add(L"u", 0b100, is_vacuum_occupied, is_reference_occupied, is_hole,
           is_particle, 5);

  return set_scoped_default_context(
      {.index_space_registry_shared_ptr = std::move(reg)});
}

TEMPLATE_LIST_TEST_CASE("export_tests", "[export]", KnownGenerators) {
  using CurrentGen = TestType;
  using CurrentCtx = CurrentGen::Context;

  auto resetter = to_export_context();

  REQUIRE(Index(L"i_1") < Index(L"a_1"));

  // Safe-guard that template magic works
  const std::size_t n_generators = 5;

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

      auto groups = parse_expression_spec(expression_spec);
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
            parse_expr(input, Symmetry::nonsymm, BraKetSymmetry::nonsymm,
                       ParticleSymmetry::nonsymm)
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
      export_expression(to_export_tree(parse_result_expr(L"v1 = 2 v2")),
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
            to_export_tree(parse_result_expr(L"v1 = 2 v2 + 4 v2 v3")),
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
        export_expression(to_export_tree(parse_result_expr(
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

      // The following test cases will only become relevant once the optimizer
      // is further improved
#if 0
      SECTION("multiple with reordering") {
        // When stripping redundant load/unload operations, load of tensor C
        // needs to be moved before the load of B in order to retain
        // compatibility to frameworks with stack-based memory models
        export_expression(
            to_export_tree(parse_result_expr(
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
        export_expression(to_export_tree(parse_result_expr(
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
            to_export_tree(parse_result_expr(
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
            to_export_tree(parse_result_expr(
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
        Tensor tensor = parse_expr(L"t{i1;a1}:N-N-N")->as<Tensor>();
        bool rewritten = ctx.rewrite(tensor);
        REQUIRE_THAT(tensor, EquivalentTo("t{i1;a1}:N-N-N"));
        REQUIRE_FALSE(rewritten);

        tensor = parse_expr(int_label + L"{a1;i1}:N-N-N")->as<Tensor>();
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
