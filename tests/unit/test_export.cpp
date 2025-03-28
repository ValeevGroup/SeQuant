#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/export/julia_itensor.hpp>
#include <SeQuant/core/export/julia_tensor_kit.hpp>
#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include "catch2_sequant.hpp"

#include <boost/algorithm/string.hpp>

#include <algorithm>
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
    JuliaTensorOperationsGenerator<JuliaTensorOperationsGeneratorContext>
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

TEMPLATE_LIST_TEST_CASE("export_tests", "[export]", KnownGenerators) {
  using CurrentGen = TestType;
  using CurrentCtx = CurrentGen::Context;

  // Safe-guard that template magic works
  const std::size_t n_generators = 4;

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
      std::string expression;
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
            expression += line;
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

      REQUIRE(!expression.empty());

      ExprPtr input = parse_expr(to_wstring(expression));
      // CAPTURE(*input);

      auto tree = binarize(input);

      export_expression(tree, generator, context);

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

  SECTION("Julia") {
    using namespace sequant;
    Index a_1 = Index("a_1");
    Index i_1 = Index("i_1");
    IndexSpace a = a_1.space();
    IndexSpace i = i_1.space();
    std::map<IndexSpace, std::string> index_tags;
    std::map<IndexSpace, std::string> index_dims;
    index_tags[a] = "v";
    index_tags[i] = "o";
    index_dims[a] = "nv";
    index_dims[i] = "nocc";

    SECTION("TensorOperations") {
      JuliaTensorOperationsGeneratorContext ctx(index_tags, index_dims);
      JuliaTensorOperationsGenerator generator;

      SECTION("represent complex") {
        Complex<rational> z1(1.0, -2.0);
        auto c = ex<Constant>(z1);
        auto tree = binarize(c * ex<Variable>(L"Dummy"));
        export_expression(tree, generator);
        std::string expected =
            "\n"
            "Z = 0.0\n"
            "Dummy = deserialize(\"Dummy.jlbin\")\n"
            "@tensor Z += (1-2im) * Dummy\n"
            "Dummy = nothing\n"
            "return Z\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("complex scalar addition") {
        // z1 and z2 are complex scalars to be read from disc
        ExprPtr expr = parse_expr(L"z1 + z2");
        auto tree = binarize(expr);
        export_expression(tree, generator);
        std::string expected =
            "\n"
            "Z = 0.0\n"
            "z1 = deserialize(\"z1.jlbin\")\n"
            "z2 = deserialize(\"z2.jlbin\")\n"
            "@tensor Z += z1 + z2\n"
            "z2 = nothing\n"
            "z1 = nothing\n"
            "return Z\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("proto indices") {
        ExprPtr expr = parse_expr(L"g{i1<a1>;}+u{i1<a1>;}");
        auto tree = binarize(expr);
        REQUIRE_THROWS_WITH(export_expression(tree, generator, ctx),
                            "Proto Indices are not (yet) supported!");
      }

      SECTION("binary contraction") {
        ExprPtr expr = parse_expr(L"g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_3,i_4}");
        auto tree = binarize(expr);

        export_expression(tree, generator, ctx);

        std::string expected =
            "\n"
            "\n"
            "I_oooo = zeros(Float64, nocc, nocc, nocc, nocc)\n"
            "g_oovv = deserialize(\"g_oovv.jlbin\")\n"
            "T2_vvoo = deserialize(\"T2_vvoo.jlbin\")\n"
            "@tensor I_oooo[ i_1, i_2, i_3, i_4 ] += g_oovv[ i_1, i_2, a_1, "
            "a_2 ] * T2_vvoo[ a_1, a_2, i_3, i_4 ]\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "return I_oooo\n";
        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("ternary") {
        ExprPtr expr = parse_expr(L"A{a2;i2} B{i2;a1} C{i1;a2}");
        auto tree = binarize(expr);

        export_expression(tree, generator, ctx);

        std::string expected =
            "\n"
            "\n"
            "I_ov = zeros(Float64, nocc, nv)\n"
            "I_vv = zeros(Float64, nv, nv)\n"
            "A_vo = deserialize(\"A_vo.jlbin\")\n"
            "B_ov = deserialize(\"B_ov.jlbin\")\n"
            "@tensor I_vv[ a_2, a_1 ] += A_vo[ a_2, i_2 ] * B_ov[ i_2, a_1 ]\n"
            "B_ov = nothing\n"
            "A_vo = nothing\n"
            "C_ov = deserialize(\"C_ov.jlbin\")\n"
            "@tensor I_ov[ i_1, a_1 ] += I_vv[ a_2, a_1 ] * C_ov[ i_1, a_2 ]\n"
            "C_ov = nothing\n"
            "I_vv = nothing\n"
            "return I_ov\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("Binary Dot , Scalar Multiplication and Sum") {
        ExprPtr expr = parse_expr(
            L"2 g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_1,i_2} - 1 "
            L"g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_2,i_1}");
        auto tree = binarize(expr);
        export_expression(tree, generator, ctx);
        std::string expected =
            "\n"
            "\n"
            "\n"
            "Z = 0.0\n"
            "g_oovv = deserialize(\"g_oovv.jlbin\")\n"
            "T2_vvoo = deserialize(\"T2_vvoo.jlbin\")\n"
            "@tensor Z += 2 * g_oovv[ i_1, i_2, a_1, a_2 ] * T2_vvoo[ a_1, "
            "a_2, i_1, i_2 ]\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "g_oovv = deserialize(\"g_oovv.jlbin\")\n"
            "T2_vvoo = deserialize(\"T2_vvoo.jlbin\")\n"
            "@tensor Z += -1 * g_oovv[ i_1, i_2, a_1, a_2 ] * T2_vvoo[ a_1, "
            "a_2, i_2, i_1 ]\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "return Z\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("binary + ternary") {
        auto tree = binarize(
            parse_expr(L"A{a1;i1} B{i1;i2} + A{a1;i1} B{i1;i3} C{i3;i2}"));

        export_expression(tree, generator, ctx);

        std::string expected =
            "\n"
            "\n"
            "I_vo = zeros(Float64, nv, nocc)\n"
            "I2_vo = zeros(Float64, nv, nocc)\n"
            "A_vo = deserialize(\"A_vo.jlbin\")\n"
            "B_oo = deserialize(\"B_oo.jlbin\")\n"
            "@tensor I2_vo[ a_1, i_3 ] += A_vo[ a_1, i_1 ] * B_oo[ i_1, i_3 ]\n"
            "B_oo = nothing\n"
            "A_vo = nothing\n"
            "C_oo = deserialize(\"C_oo.jlbin\")\n"
            "@tensor I_vo[ a_1, i_2 ] += I2_vo[ a_1, i_3 ] * C_oo[ i_3, i_2 ]\n"
            "C_oo = nothing\n"
            "I2_vo = nothing\n"
            "A_vo = deserialize(\"A_vo.jlbin\")\n"
            "B_oo = deserialize(\"B_oo.jlbin\")\n"
            "@tensor I_vo[ a_1, i_2 ] += A_vo[ a_1, i_1 ] * B_oo[ i_1, i_2 ]\n"
            "B_oo = nothing\n"
            "A_vo = nothing\n"
            "return I_vo\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
    }

    SECTION("TensorKit") {
      JuliaTensorKitGeneratorContext ctx(index_tags, index_dims);
      JuliaTensorKitGenerator generator;

      SECTION("binary contraction") {
        ExprPtr expr = parse_expr(L"g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_3,i_4}");
        auto tree = binarize(expr);

        export_expression(tree, generator, ctx);

        std::string expected =
            "\n"
            "\n"
            "I_oooo = TensorMap(zeros(Float64, nocc, nocc, nocc, nocc), ℝ^nocc "
            "⊗ ℝ^nocc, ℝ^nocc ⊗ ℝ^nocc)\n"
            "g_oovv = TensorMap(deserialize(\"g_oovv.jlbin\"), ℝ^nocc ⊗ "
            "ℝ^nocc, ℝ^nv ⊗ ℝ^nv)\n"
            "T2_vvoo = TensorMap(deserialize(\"T2_vvoo.jlbin\"), ℝ^nv ⊗ "
            "ℝ^nv, ℝ^nocc ⊗ ℝ^nocc)\n"
            "@tensor I_oooo[ i_1, i_2, i_3, i_4 ] += g_oovv[ i_1, i_2, a_1, "
            "a_2 ] * T2_vvoo[ a_1, a_2, i_3, i_4 ]\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "return I_oooo\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("ternary") {
        ExprPtr expr = parse_expr(L"A{a2;i2} B{i2;a1} C{i1;a2}");
        auto tree = binarize(expr);

        export_expression(tree, generator, ctx);

        std::string expected =
            "\n"
            "\n"
            "I_ov = TensorMap(zeros(Float64, nocc, nv), ℝ^nocc, ℝ^nv)\n"
            "I_vv = TensorMap(zeros(Float64, nv, nv), ℝ^nv, ℝ^nv)\n"
            "A_vo = TensorMap(deserialize(\"A_vo.jlbin\"), ℝ^nv, ℝ^nocc)\n"
            "B_ov = TensorMap(deserialize(\"B_ov.jlbin\"), ℝ^nocc, ℝ^nv)\n"
            "@tensor I_vv[ a_2, a_1 ] += A_vo[ a_2, i_2 ] * B_ov[ i_2, a_1 ]\n"
            "B_ov = nothing\n"
            "A_vo = nothing\n"
            "C_ov = TensorMap(deserialize(\"C_ov.jlbin\"), ℝ^nocc, ℝ^nv)\n"
            "@tensor I_ov[ i_1, a_1 ] += I_vv[ a_2, a_1 ] * C_ov[ i_1, a_2 ]\n"
            "C_ov = nothing\n"
            "I_vv = nothing\n"
            "return I_ov\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("Binary Dot , Scalar Multiplication and Sum") {
        ExprPtr expr = parse_expr(
            L"2 g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_1,i_2} - 1 "
            L"g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_2,i_1}");
        auto tree = binarize(expr);
        export_expression(tree, generator, ctx);
        std::string expected =
            "\n"
            "\n"
            "\n"
            "Z = 0.0\n"
            "g_oovv = TensorMap(deserialize(\"g_oovv.jlbin\"), ℝ^nocc ⊗ "
            "ℝ^nocc, ℝ^nv ⊗ ℝ^nv)\n"
            "T2_vvoo = TensorMap(deserialize(\"T2_vvoo.jlbin\"), ℝ^nv ⊗ "
            "ℝ^nv, ℝ^nocc ⊗ ℝ^nocc)\n"
            "@tensor Z += 2 * g_oovv[ i_1, i_2, a_1, a_2 ] * T2_vvoo[ a_1, "
            "a_2, i_1, i_2 ]\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "g_oovv = TensorMap(deserialize(\"g_oovv.jlbin\"), ℝ^nocc ⊗ "
            "ℝ^nocc, ℝ^nv ⊗ ℝ^nv)\n"
            "T2_vvoo = TensorMap(deserialize(\"T2_vvoo.jlbin\"), ℝ^nv ⊗ "
            "ℝ^nv, ℝ^nocc ⊗ ℝ^nocc)\n"
            "@tensor Z += -1 * g_oovv[ i_1, i_2, a_1, a_2 ] * T2_vvoo[ a_1, "
            "a_2, i_2, i_1 ]\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "return Z\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("binary + ternary") {
        auto tree = binarize(
            parse_expr(L"A{a1;i1} B{i1;i2} + A{a1;i1} B{i1;i3} C{i3;i2}"));

        export_expression(tree, generator, ctx);

        std::string expected =
            "\n"
            "\n"
            "I_vo = TensorMap(zeros(Float64, nv, nocc), ℝ^nv, ℝ^nocc)\n"
            "I2_vo = TensorMap(zeros(Float64, nv, nocc), ℝ^nv, ℝ^nocc)\n"
            "A_vo = TensorMap(deserialize(\"A_vo.jlbin\"), ℝ^nv, ℝ^nocc)\n"
            "B_oo = TensorMap(deserialize(\"B_oo.jlbin\"), ℝ^nocc, ℝ^nocc)\n"
            "@tensor I2_vo[ a_1, i_3 ] += A_vo[ a_1, i_1 ] * B_oo[ i_1, i_3 ]\n"
            "B_oo = nothing\n"
            "A_vo = nothing\n"
            "C_oo = TensorMap(deserialize(\"C_oo.jlbin\"), ℝ^nocc, ℝ^nocc)\n"
            "@tensor I_vo[ a_1, i_2 ] += I2_vo[ a_1, i_3 ] * C_oo[ i_3, i_2 ]\n"
            "C_oo = nothing\n"
            "I2_vo = nothing\n"
            "A_vo = TensorMap(deserialize(\"A_vo.jlbin\"), ℝ^nv, ℝ^nocc)\n"
            "B_oo = TensorMap(deserialize(\"B_oo.jlbin\"), ℝ^nocc, ℝ^nocc)\n"
            "@tensor I_vo[ a_1, i_2 ] += A_vo[ a_1, i_1 ] * B_oo[ i_1, i_2 ]\n"
            "B_oo = nothing\n"
            "A_vo = nothing\n"
            "return I_vo\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
    }

    SECTION("ITensor") {
      JuliaITensorGeneratorContext ctx(index_tags, index_dims);
      JuliaITensorGenerator generator;

      SECTION("binary contraction") {
        ExprPtr expr = parse_expr(L"g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_3,i_4}");
        auto tree = binarize(expr);

        export_expression(tree, generator, ctx);

        std::string expected =
            "i_1 = Index(nocc, \"i_1\")\n"
            "i_2 = Index(nocc, \"i_2\")\n"
            "i_3 = Index(nocc, \"i_3\")\n"
            "i_4 = Index(nocc, \"i_4\")\n"
            "a_1 = Index(nv, \"a_1\")\n"
            "a_2 = Index(nv, \"a_2\")\n"
            "\n"
            "\n"
            "I_oooo = ITensor(zeros(Float64, nocc, nocc, nocc, nocc), i_1, "
            "i_2, i_3, i_4)\n"
            "g_oovv = ITensor(deserialize(\"g_oovv.jlbin\"), i_1, i_2, a_1, "
            "a_2)\n"
            "T2_vvoo = ITensor(deserialize(\"T2_vvoo.jlbin\"), a_1, a_2, i_3, "
            "i_4)\n"
            "I_oooo += g_oovv * T2_vvoo\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "return I_oooo\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("ternary") {
        ExprPtr expr = parse_expr(L"A{a2;i2} B{i2;a1} C{i1;a2}");
        auto tree = binarize(expr);

        export_expression(tree, generator, ctx);

        std::string expected =
            "i_1 = Index(nocc, \"i_1\")\n"
            "i_2 = Index(nocc, \"i_2\")\n"
            "a_1 = Index(nv, \"a_1\")\n"
            "a_2 = Index(nv, \"a_2\")\n"
            "\n"
            "\n"
            "I_ov = ITensor(zeros(Float64, nocc, nv), i_1, a_1)\n"
            "I_vv = ITensor(zeros(Float64, nv, nv), a_2, a_1)\n"
            "A_vo = ITensor(deserialize(\"A_vo.jlbin\"), a_2, i_2)\n"
            "B_ov = ITensor(deserialize(\"B_ov.jlbin\"), i_2, a_1)\n"
            "I_vv += A_vo * B_ov\n"
            "B_ov = nothing\n"
            "A_vo = nothing\n"
            "C_ov = ITensor(deserialize(\"C_ov.jlbin\"), i_1, a_2)\n"
            "I_ov += I_vv * C_ov\n"
            "C_ov = nothing\n"
            "I_vv = nothing\n"
            "return I_ov\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("Binary Dot , Scalar Multiplication and Sum") {
        ExprPtr expr = parse_expr(
            L"2 g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_1,i_2} - 1 "
            L"g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_2,i_1}");
        auto tree = binarize(expr);
        export_expression(tree, generator, ctx);
        std::string expected =
            "i_1 = Index(nocc, \"i_1\")\n"
            "i_2 = Index(nocc, \"i_2\")\n"
            "a_1 = Index(nv, \"a_1\")\n"
            "a_2 = Index(nv, \"a_2\")\n"
            "\n"
            "\n"
            "\n"
            "tmpvar = 0.0\n"
            "Z = ITensor(tmpvar)\n"
            "g_oovv = ITensor(deserialize(\"g_oovv.jlbin\"), i_1, i_2, a_1, "
            "a_2)\n"
            "T2_vvoo = ITensor(deserialize(\"T2_vvoo.jlbin\"), a_1, a_2, i_1, "
            "i_2)\n"
            "Z += 2 * g_oovv * T2_vvoo\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "g_oovv = ITensor(deserialize(\"g_oovv.jlbin\"), i_1, i_2, a_1, "
            "a_2)\n"
            "T2_vvoo = ITensor(deserialize(\"T2_vvoo.jlbin\"), a_1, a_2, i_2, "
            "i_1)\n"
            "Z += -1 * g_oovv * T2_vvoo\n"
            "T2_vvoo = nothing\n"
            "g_oovv = nothing\n"
            "return Z\n";
        REQUIRE(generator.get_generated_code() == expected);
      }

      SECTION("binary + ternary") {
        auto tree = binarize(
            parse_expr(L"A{a1;i1} B{i1;i2} + A{a1;i1} B{i1;i3} C{i3;i2}"));

        export_expression(tree, generator, ctx);

        std::string expected =
            "i_1 = Index(nocc, \"i_1\")\n"
            "i_2 = Index(nocc, \"i_2\")\n"
            "i_3 = Index(nocc, \"i_3\")\n"
            "a_1 = Index(nv, \"a_1\")\n"
            "\n"
            "\n"
            "I_vo = ITensor(zeros(Float64, nv, nocc), a_1, i_2)\n"
            "I2_vo = ITensor(zeros(Float64, nv, nocc), a_1, i_3)\n"
            "A_vo = ITensor(deserialize(\"A_vo.jlbin\"), a_1, i_1)\n"
            "B_oo = ITensor(deserialize(\"B_oo.jlbin\"), i_1, i_3)\n"
            "I2_vo += A_vo * B_oo\n"
            "B_oo = nothing\n"
            "A_vo = nothing\n"
            "C_oo = ITensor(deserialize(\"C_oo.jlbin\"), i_3, i_2)\n"
            "I_vo += I2_vo * C_oo\n"
            "C_oo = nothing\n"
            "I2_vo = nothing\n"
            "A_vo = ITensor(deserialize(\"A_vo.jlbin\"), a_1, i_1)\n"
            "B_oo = ITensor(deserialize(\"B_oo.jlbin\"), i_1, i_2)\n"
            "I_vo += A_vo * B_oo\n"
            "B_oo = nothing\n"
            "A_vo = nothing\n"
            "return I_vo\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
    }
  }
}
