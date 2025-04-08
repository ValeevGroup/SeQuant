#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/expression_group.hpp>
#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/export/itf_generator.hpp>
#include <SeQuant/core/export/julia_itensor.hpp>
#include <SeQuant/core/export/julia_tensor_kit.hpp>
#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/utility/string.hpp>
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

class ItfContext : public sequant::itf::Context {
 public:
  ItfContext() = default;

  int compare(const sequant::Index &lhs, const sequant::Index &rhs) const {
    return rhs.space().type().to_int32() - lhs.space().type().to_int32();
  }

  std::wstring get_base_label(const sequant::IndexSpace &space) const {
    return L"a";
  }
  std::wstring get_tag(const sequant::IndexSpace &space) const { return L"b"; }
  std::wstring get_name(const sequant::IndexSpace &space) const {
    return L"tenshi";
  }
};

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
	ItfGenerator<ItfGeneratorContext>
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

void configure_context_defaults(ItfGeneratorContext &ctx) {
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
void add_to_context(ItfGeneratorContext &ctx, std::string_view key,
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

std::vector<ExpressionGroup<EvalExpr>> parse_expression_spec(
    const std::string &spec) {
  std::vector<ExpressionGroup<EvalExpr>> groups;

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
      ResultExpr res = parse_result_expr(to_wstring(line));
      groups.back().add(binarize(res));
    } catch (...) {
      ExprPtr expr = parse_expr(to_wstring(line));
      groups.back().add(binarize(expr));
    }
  }

  return groups;
}

TEMPLATE_LIST_TEST_CASE("export_tests", "[export]", KnownGenerators) {
  using CurrentGen = TestType;
  using CurrentCtx = CurrentGen::Context;

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

      export_groups<EvalExpr>(groups, generator, context);

      REQUIRE_THAT(generator.get_generated_code(),
                   DiffedStringEquals(expected_output.value()));
    }
  }
}

TEST_CASE("export", "[export]") {
  SECTION("itf") {
    const ItfContext ctx;

    SECTION("remap_integrals") {
      using namespace sequant::itf::detail;
      SECTION("Unchanged") {
        auto expr = parse_expr(L"t{i1;a1}");
        auto remapped = expr;
        remap_integrals(expr, ctx);
        REQUIRE(remapped == expr);

        expr = parse_expr(L"t{i1;a1} f{a1;i1} + first{a1;i1} second{i1;a1}");
        remapped = expr;
        remap_integrals(expr, ctx);
        REQUIRE(remapped == expr);
      }

      SECTION("K") {
        SECTION("occ,occ,occ,occ") {
          std::vector<Index> indices = {L"i_1", L"i_2", L"i_3", L"i_4"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            REQUIRE(indexPerm.size() == 4);

            ExprPtr integralExpr = ex<Tensor>(
                L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            auto transformed = integralExpr;
            remap_integrals(transformed, ctx);

            CAPTURE(indexPerm);
            REQUIRE_THAT(transformed, EquivalentTo("K{i1,i2;i3,i4}"));
          }
        }

        SECTION("virt,virt,occ,occ") {
          std::vector<Index> indices = {L"a_1", L"a_2", L"i_1", L"i_2"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            REQUIRE(indexPerm.size() == 4);

            ExprPtr integralExpr = ex<Tensor>(
                L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            auto transformed = integralExpr;
            remap_integrals(transformed, ctx);

            CAPTURE(indexPerm);

            REQUIRE_THAT(transformed, EquivalentTo("K{a1,a2;i1,i2}"));
          }
        }

        SECTION("virt,virt,virt,virt") {
          std::vector<Index> indices = {L"a_1", L"a_2", L"a_3", L"a_4"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            REQUIRE(indexPerm.size() == 4);

            ExprPtr integralExpr = ex<Tensor>(
                L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            auto transformed = integralExpr;
            remap_integrals(transformed, ctx);

            CAPTURE(indexPerm);

            REQUIRE_THAT(transformed, EquivalentTo("K{a1,a2;a3,a4}"));
          }
        }
      }

      SECTION("J") {
        SECTION("virt,occ,virt,occ") {
          std::vector<Index> indices = {L"a_1", L"i_1", L"a_2", L"i_2"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            REQUIRE(indexPerm.size() == 4);

            ExprPtr integralExpr = ex<Tensor>(
                L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            auto transformed = integralExpr;
            remap_integrals(transformed, ctx);

            CAPTURE(indexPerm);

            REQUIRE_THAT(transformed, EquivalentTo("J{a1,a2;i1,i2}"));
          }
        }

        SECTION("virt,occ,virt,virt") {
          std::vector<Index> indices = {L"a_1", L"i_1", L"a_2", L"a_3"};
          REQUIRE(indices.size() == 4);

          for (const std::vector<std::size_t> &indexPerm :
               twoElectronIntegralSymmetries()) {
            REQUIRE(indexPerm.size() == 4);

            ExprPtr integralExpr = ex<Tensor>(
                L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
                ket{indices[indexPerm[2]], indices[indexPerm[3]]});

            auto transformed = integralExpr;
            remap_integrals(transformed, ctx);

            CAPTURE(indexPerm);

            REQUIRE_THAT(transformed, EquivalentTo("J{a1,a2;a3,i1}"));
          }
        }
      }

      SECTION("f") {
        SECTION("same_space") {
          ExprPtr expr = parse_expr(L"f{i2;i1} + f{a1;a2}");

          remap_integrals(expr, ctx);

          REQUIRE_THAT(expr, EquivalentTo("f{i1;i2} + f{a1;a2}"));
        }
        SECTION("different_space") {
          ExprPtr expr = parse_expr(L"f{a1;i1} + f{i1;a1}");

          remap_integrals(expr, ctx);

          REQUIRE_THAT(expr, EquivalentTo("f{a1;i1} + f{a1;i1}"));
        }
      }
    }
  }
}
