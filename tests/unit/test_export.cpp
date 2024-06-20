#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/parse.hpp>

#include <string>
#include <vector>

namespace Catch {

// Note: Again, template specialization doesn't seem to be used from inside
// ::Catch::Details::stringify for some reason
template <>
struct StringMaker<sequant::ExprPtr> {
  static std::string convert(const sequant::ExprPtr &expr) {
    return sequant::toUtf8(sequant::deparse(expr, false));
  }
};

}  // namespace Catch

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

#define CAPTURE_EXPR(expr) \
  INFO(#expr " := " << ::Catch::StringMaker<sequant::ExprPtr>::convert(expr))

TEST_CASE("Export capabilities", "[exports]") {
  using namespace sequant;

  SECTION("Tree processing") {
    TextGenerator generator;
    SECTION("basics") {
      SECTION("binary_dot") {
        ExprPtr expr = parse_expr(L"A{a1;i1} B{i1;a1}");
        auto tree = eval_node<EvalExpr>(expr);

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index a_1\n"
            "\n"
            "Declare variable Z\n"
            "\n"
            "Declare tensor A[a_1, i_1]\n"
            "Declare tensor B[i_1, a_1]\n"
            "\n"
            "Load Z\n"
            "Load A[a_1, i_1]\n"
            "Load B[i_1, a_1]\n"
            "Compute Z += A[a_1, i_1] B[i_1, a_1]\n"
            "Unload B[i_1, a_1]\n"
            "Unload A[a_1, i_1]\n"
            "Persist Z\n";
        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("binary_contract") {
        ExprPtr expr = parse_expr(L"A{a1;i2} B{i1;a1}");
        auto tree = eval_node<EvalExpr>(expr);

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index i_2\n"
            "Declare index a_1\n"
            "\n"
            "Declare tensor A[a_1, i_2]\n"
            "Declare tensor B[i_1, a_1]\n"
            "Declare tensor I[i_1, i_2]\n"
            "\n"
            "Create I[i_1, i_2] and initialize to zero\n"
            "Load A[a_1, i_2]\n"
            "Load B[i_1, a_1]\n"
            "Compute I[i_1, i_2] += A[a_1, i_2] B[i_1, a_1]\n"
            "Unload B[i_1, a_1]\n"
            "Unload A[a_1, i_2]\n"
            "Persist I[i_1, i_2]\n";
        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("ternary") {
        ExprPtr expr = parse_expr(L"A{a2;i2} B{i2;a1} C{i1;a2}");
        auto tree = eval_node<EvalExpr>(expr);

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index i_2\n"
            "Declare index a_1\n"
            "Declare index a_2\n"
            "\n"
            "Declare tensor A[a_2, i_2]\n"
            "Declare tensor B[i_2, a_1]\n"
            "Declare tensor C[i_1, a_2]\n"
            "Declare tensor I[i_1, a_1]\n"
            "Declare tensor I[a_2, a_1]\n"
            "\n"
            "Create I[i_1, a_1] and initialize to zero\n"
            "Create I[a_2, a_1] and initialize to zero\n"
            "Load A[a_2, i_2]\n"
            "Load B[i_2, a_1]\n"
            "Compute I[a_2, a_1] += A[a_2, i_2] B[i_2, a_1]\n"
            "Unload B[i_2, a_1]\n"
            "Unload A[a_2, i_2]\n"
            "Load C[i_1, a_2]\n"
            "Compute I[i_1, a_1] += I[a_2, a_1] C[i_1, a_2]\n"
            "Unload C[i_1, a_2]\n"
            "Unload I[a_2, a_1]\n"
            "Persist I[i_1, a_1]\n";
        REQUIRE(generator.get_generated_code() == expected);
      }
    }
    SECTION("scalar * tensor") {
      ExprPtr expr = parse_expr(L"42 * A{a1;i1}");
      auto tree = eval_node<EvalExpr>(expr);

      export_expression(tree, generator);

      std::string expected =
          "Declare index i_1\n"
          "Declare index a_1\n"
          "\n"
          "Declare tensor A[a_1, i_1]\n"
          "Declare tensor I[a_1, i_1]\n"
          "\n"
          "Create I[a_1, i_1] and initialize to zero\n"
          "Load A[a_1, i_1]\n"
          "Compute I[a_1, i_1] += 42 A[a_1, i_1]\n"
          "Unload A[a_1, i_1]\n"
          "Persist I[a_1, i_1]\n";

      REQUIRE(generator.get_generated_code() == expected);
    }
    SECTION("constant scalar factor") {
      ExprPtr expr = parse_expr(L"5 * A{a2;i2} B{i2;a1}");
      auto tree = eval_node<EvalExpr>(expr);

      export_expression(tree, generator);

      std::string expected =
          "Declare index i_2\n"
          "Declare index a_1\n"
          "Declare index a_2\n"
          "\n"
          "Declare tensor A[a_2, i_2]\n"
          "Declare tensor B[i_2, a_1]\n"
          "Declare tensor I[a_2, a_1]\n"
          "\n"
          "Create I[a_2, a_1] and initialize to zero\n"
          "Load A[a_2, i_2]\n"
          "Load B[i_2, a_1]\n"
          "Compute I[a_2, a_1] += 5 A[a_2, i_2] B[i_2, a_1]\n"
          "Unload B[i_2, a_1]\n"
          "Unload A[a_2, i_2]\n"
          "Persist I[a_2, a_1]\n";

      REQUIRE(generator.get_generated_code() == expected);
    }
    SECTION("variable scalar factor") {
      ExprPtr expr = parse_expr(L"myVar * A{a2;i2} B{i2;a1}");
      auto tree = eval_node<EvalExpr>(expr);

      export_expression(tree, generator);

      std::string expected =
          "Declare index i_2\n"
          "Declare index a_1\n"
          "Declare index a_2\n"
          "\n"
          "Declare variable myVar\n"
          "\n"
          "Declare tensor A[a_2, i_2]\n"
          "Declare tensor B[i_2, a_1]\n"
          "Declare tensor I[a_2, a_1]\n"
          "\n"
          "Create I[a_2, a_1] and initialize to zero\n"
          "Load A[a_2, i_2]\n"
          "Load B[i_2, a_1]\n"
          "Load myVar\n"
          "Compute I[a_2, a_1] += myVar A[a_2, i_2] B[i_2, a_1]\n"
          "Unload myVar\n"
          "Unload B[i_2, a_1]\n"
          "Unload A[a_2, i_2]\n"
          "Persist I[a_2, a_1]\n";

      REQUIRE(generator.get_generated_code() == expected);
    }
    SECTION("tree rebalancing") {
      ExprPtr first = parse_expr(L"A{a1;i1} B{a2;a1}");
      ExprPtr second = parse_expr(L"C{i1;a3}");

      /*
       * Note that we are creating the expression tree
       *
       *                    R
       *                  /   \
       *                 C     I
       *                     /   \
       *                    A     B
       * which is suboptimal for evaluation as this requires (for a backend that
       * can only deal with stack-like memory allocations) to have R, C, I, A
       * and B in memory all at the same time. Thus, this tree has to be
       * rearranged to
       *                    R
       *                  /   \
       *                 I     C
       *               /   \
       *              A     B
       * which allows to not have C in memory while computing I from A and B,
       * even in with stack-like memory model.
       */
      auto tree = eval_node<EvalExpr>(
          ex<Product>(ExprPtrList{second, first}, Product::Flatten::No));

      export_expression(tree, generator);

      std::string expected =
          "Declare index i_1\n"
          "Declare index a_1\n"
          "Declare index a_2\n"
          "Declare index a_3\n"
          "\n"
          "Declare tensor A[a_1, i_1]\n"
          "Declare tensor B[a_2, a_1]\n"
          "Declare tensor C[i_1, a_3]\n"
          "Declare tensor I[a_2, i_1]\n"
          "Declare tensor I[a_2, a_3]\n"
          "\n"
          "Create I[a_2, a_3] and initialize to zero\n"
          "Create I[a_2, i_1] and initialize to zero\n"
          "Load A[a_1, i_1]\n"
          "Load B[a_2, a_1]\n"
          "Compute I[a_2, i_1] += A[a_1, i_1] B[a_2, a_1]\n"
          "Unload B[a_2, a_1]\n"
          "Unload A[a_1, i_1]\n"
          "Load C[i_1, a_3]\n"
          "Compute I[a_2, a_3] += I[a_2, i_1] C[i_1, a_3]\n"
          "Unload C[i_1, a_3]\n"
          "Unload I[a_2, i_1]\n"
          "Persist I[a_2, a_3]\n";

      REQUIRE(generator.get_generated_code() == expected);
    }
    SECTION("duplicate tensor blocks") {
      SECTION("non-overlapping") {
        // Computation of first and second leads to an intermediate of the same
        // shape/dimension. Since these two results are not needed
        // simultaneously, we can allocate a single tensor for one and reuse it
        // for the other (after having reset it to zero before)
        ExprPtr first = parse_expr(L"A{a1;i2} B{i1;a1}");
        ExprPtr second = parse_expr(L"C{a2;i4} D{i3;a2}");
        ExprPtr third = parse_expr(L"E{i2;a3}");

        auto tree = eval_node<EvalExpr>(ex<Product>(
            ExprPtrList{
                ex<Product>(ExprPtrList{first, third}, Product::Flatten::No),
                second},
            Product::Flatten::No));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index i_2\n"
            "Declare index i_3\n"
            "Declare index i_4\n"
            "Declare index a_1\n"
            "Declare index a_2\n"
            "Declare index a_3\n"
            "\n"
            "Declare tensor A[a_1, i_2]\n"
            "Declare tensor B[i_1, a_1]\n"
            "Declare tensor C[a_2, i_4]\n"
            "Declare tensor D[i_3, a_2]\n"
            "Declare tensor E[i_2, a_3]\n"
            "Declare tensor I[i_1, i_2]\n"
            "Declare tensor I[i_1, a_3]\n"
            "Declare tensor I[i_1, i_3, a_3, i_4]\n"
            "\n"
            "Create I[i_1, i_3, a_3, i_4] and initialize to zero\n"
            "Create I[i_1, a_3] and initialize to zero\n"
            "Create I[i_1, i_2] and initialize to zero\n"
            "Load A[a_1, i_2]\n"
            "Load B[i_1, a_1]\n"
            "Compute I[i_1, i_2] += A[a_1, i_2] B[i_1, a_1]\n"
            "Unload B[i_1, a_1]\n"
            "Unload A[a_1, i_2]\n"
            "Load E[i_2, a_3]\n"
            "Compute I[i_1, a_3] += I[i_1, i_2] E[i_2, a_3]\n"
            "Unload E[i_2, a_3]\n"
            "Unload I[i_1, i_2]\n"
            // TODO: Instead of unloading I and then loading
            // it again, we could simply keep it loaded instead.
            "Load I[i_3, i_4] and set it to zero\n"
            "Load C[a_2, i_4]\n"
            "Load D[i_3, a_2]\n"
            "Compute I[i_3, i_4] += C[a_2, i_4] D[i_3, a_2]\n"
            "Unload D[i_3, a_2]\n"
            "Unload C[a_2, i_4]\n"
            "Compute I[i_1, i_3, a_3, i_4] += I[i_1, a_3] I[i_3, i_4]\n"
            "Unload I[i_3, i_4]\n"
            "Unload I[i_1, a_3]\n"
            "Persist I[i_1, i_3, a_3, i_4]\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("overlapping") {
        ExprPtr first = parse_expr(L"A{a1;i1} B{i2;a1}");
        ExprPtr second = parse_expr(L"C{i1;i3}");

        auto tree = eval_node<EvalExpr>(
            ex<Product>(ExprPtrList{first, second}, Product::Flatten::No));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index i_2\n"
            "Declare index i_3\n"
            "Declare index a_1\n"
            "\n"
            "Declare tensor A[a_1, i_1]\n"
            "Declare tensor B[i_2, a_1]\n"
            "Declare tensor C[i_1, i_3]\n"
            "Declare tensor I[i_2, i_3]\n"
            // Note that in this case both I tensors have to be used
            // at the same time (at some point) and therefore we have
            // to  have different tensors, even though they have the
            // same shape.
            "Declare tensor I2[i_2, i_1]\n"
            "\n"
            "Create I[i_2, i_3] and initialize to zero\n"
            "Create I2[i_2, i_1] and initialize to zero\n"
            "Load A[a_1, i_1]\n"
            "Load B[i_2, a_1]\n"
            "Compute I2[i_2, i_1] += A[a_1, i_1] B[i_2, a_1]\n"
            "Unload B[i_2, a_1]\n"
            "Unload A[a_1, i_1]\n"
            "Load C[i_1, i_3]\n"
            "Compute I[i_2, i_3] += I2[i_2, i_1] C[i_1, i_3]\n"
            "Unload C[i_1, i_3]\n"
            "Unload I2[i_2, i_1]\n"
            "Persist I[i_2, i_3]\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("Duplicate leaf tensors") {
        auto tree = eval_node<EvalExpr>(parse_expr(L"A{a1;i1} A{a1;i1}"));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index a_1\n"
            "\n"
            "Declare tensor A[a_1, i_1]\n"
            "Declare tensor I[a_1, a_1, i_1, i_1]\n"
            "\n"
            "Create I[a_1, a_1, i_1, i_1] and initialize to zero\n"
            // Note that we only load A once, even though we use it twice
            "Load A[a_1, i_1]\n"
            "Compute I[a_1, a_1, i_1, i_1] += A[a_1, i_1] A[a_1, i_1]\n"
            "Unload A[a_1, i_1]\n"
            "Persist I[a_1, a_1, i_1, i_1]\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
    }
    SECTION("sum") {
      SECTION("unary + unary") {
        auto tree = eval_node<EvalExpr>(parse_expr(L"A{a1;i1} + B{a1;i1}"));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index a_1\n"
            "\n"
            "Declare tensor A[a_1, i_1]\n"
            "Declare tensor B[a_1, i_1]\n"
            "Declare tensor I[a_1, i_1]\n"
            "\n"
            "Create I[a_1, i_1] and initialize to zero\n"
            "Load A[a_1, i_1]\n"
            "Load B[a_1, i_1]\n"
            "Compute I[a_1, i_1] += A[a_1, i_1] + B[a_1, i_1]\n"
            "Unload B[a_1, i_1]\n"
            "Unload A[a_1, i_1]\n"
            "Persist I[a_1, i_1]\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("binary + binary") {
        auto tree = eval_node<EvalExpr>(
            parse_expr(L"A{a1;i1} B{a2;i2} + A{a1;i1} B{a2;i2}"));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index i_2\n"
            "Declare index a_1\n"
            "Declare index a_2\n"
            "\n"
            "Declare tensor A[a_1, i_1]\n"
            "Declare tensor B[a_2, i_2]\n"
            "Declare tensor I[a_2, a_1, i_2, i_1]\n"
            "\n"
            "Create I[a_2, a_1, i_2, i_1] and initialize to zero\n"
            "Load A[a_1, i_1]\n"
            "Load B[a_2, i_2]\n"
            "Compute I[a_2, a_1, i_2, i_1] += A[a_1, i_1] B[a_2, i_2]\n"
            "Unload B[a_2, i_2]\n"
            "Unload A[a_1, i_1]\n"
            // TODO: there is some unload/load pairs that could be removed
            "Load A[a_1, i_1]\n"
            "Load B[a_2, i_2]\n"
            "Compute I[a_2, a_1, i_2, i_1] += A[a_1, i_1] B[a_2, i_2]\n"
            "Unload B[a_2, i_2]\n"
            "Unload A[a_1, i_1]\n"
            "Persist I[a_2, a_1, i_2, i_1]\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("binary + unary") {
        auto tree =
            eval_node<EvalExpr>(parse_expr(L"A{a1;i2} B{i2;i1} + C{a1;i1}"));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index i_2\n"
            "Declare index a_1\n"
            "\n"
            "Declare tensor A[a_1, i_2]\n"
            "Declare tensor B[i_2, i_1]\n"
            "Declare tensor C[a_1, i_1]\n"
            "Declare tensor I[a_1, i_1]\n"
            "\n"
            "Create I[a_1, i_1] and initialize to zero\n"
            "Load A[a_1, i_2]\n"
            "Load B[i_2, i_1]\n"
            "Compute I[a_1, i_1] += A[a_1, i_2] B[i_2, i_1]\n"
            "Unload B[i_2, i_1]\n"
            "Unload A[a_1, i_2]\n"
            "Load C[a_1, i_1]\n"
            "Compute I[a_1, i_1] += C[a_1, i_1]\n"
            "Unload C[a_1, i_1]\n"
            "Persist I[a_1, i_1]\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("unary + unary + unary") {
        auto tree =
            eval_node<EvalExpr>(parse_expr(L"A{a1;i1} + B{a1;i1} + C{a1;i1}"));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index a_1\n"
            "\n"
            "Declare tensor A[a_1, i_1]\n"
            "Declare tensor B[a_1, i_1]\n"
            "Declare tensor C[a_1, i_1]\n"
            "Declare tensor I[a_1, i_1]\n"
            "\n"
            "Create I[a_1, i_1] and initialize to zero\n"
            "Load A[a_1, i_1]\n"
            "Load B[a_1, i_1]\n"
            // TODO: Doing I += A and I += B separately would reduce
            // the amount of tensors that have to be loaded simultaneously
            "Compute I[a_1, i_1] += A[a_1, i_1] + B[a_1, i_1]\n"
            "Unload B[a_1, i_1]\n"
            "Unload A[a_1, i_1]\n"
            "Load C[a_1, i_1]\n"
            "Compute I[a_1, i_1] += C[a_1, i_1]\n"
            "Unload C[a_1, i_1]\n"
            "Persist I[a_1, i_1]\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
    }
  }

  SECTION("remap_integrals") {
    using namespace sequant::itf::detail;
    SECTION("Unchanged") {
      auto expr = parse_expr(L"t{i1;a1}");
      auto remapped = expr;
      remap_integrals(expr);
      REQUIRE(remapped == expr);

      expr = parse_expr(L"t{i1;a1} f{a1;i1} + first{a1;i1} second{i1;a1}");
      remapped = expr;
      remap_integrals(expr);
      REQUIRE(remapped == expr);
    }

    SECTION("K") {
      SECTION("occ,occ,occ,occ") {
        std::vector<Index> indices = {L"i_1", L"i_2", L"i_3", L"i_4"};
        REQUIRE(indices.size() == 4);
        const ExprPtr expected = parse_expr(L"K{i1,i2;i3,i4}");

        for (const std::vector<std::size_t> &indexPerm :
             twoElectronIntegralSymmetries()) {
          REQUIRE(indexPerm.size() == 4);

          ExprPtr integralExpr = ex<Tensor>(
              L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
              ket{indices[indexPerm[2]], indices[indexPerm[3]]});

          auto transformed = integralExpr;
          remap_integrals(transformed);

          CAPTURE(indexPerm);
          CAPTURE_EXPR(integralExpr);
          CAPTURE_EXPR(transformed);
          CAPTURE_EXPR(expected);

          REQUIRE(transformed == expected);
        }
      }

      SECTION("virt,virt,occ,occ") {
        std::vector<Index> indices = {L"a_1", L"a_2", L"i_1", L"i_2"};
        REQUIRE(indices.size() == 4);
        const ExprPtr expected = parse_expr(L"K{a1,a2;i1,i2}");

        for (const std::vector<std::size_t> &indexPerm :
             twoElectronIntegralSymmetries()) {
          REQUIRE(indexPerm.size() == 4);

          ExprPtr integralExpr = ex<Tensor>(
              L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
              ket{indices[indexPerm[2]], indices[indexPerm[3]]});

          auto transformed = integralExpr;
          remap_integrals(transformed);

          CAPTURE(indexPerm);
          CAPTURE_EXPR(integralExpr);
          CAPTURE_EXPR(transformed);
          CAPTURE_EXPR(expected);

          REQUIRE(transformed == expected);
        }
      }

      SECTION("virt,virt,virt,virt") {
        std::vector<Index> indices = {L"a_1", L"a_2", L"a_3", L"a_4"};
        REQUIRE(indices.size() == 4);
        const ExprPtr expected = parse_expr(L"K{a1,a2;a3,a4}");

        for (const std::vector<std::size_t> &indexPerm :
             twoElectronIntegralSymmetries()) {
          REQUIRE(indexPerm.size() == 4);

          ExprPtr integralExpr = ex<Tensor>(
              L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
              ket{indices[indexPerm[2]], indices[indexPerm[3]]});

          auto transformed = integralExpr;
          remap_integrals(transformed);

          CAPTURE(indexPerm);
          CAPTURE_EXPR(integralExpr);
          CAPTURE_EXPR(transformed);
          CAPTURE_EXPR(expected);

          REQUIRE(transformed == expected);
        }
      }
    }

    SECTION("J") {
      SECTION("virt,occ,virt,occ") {
        std::vector<Index> indices = {L"a_1", L"i_1", L"a_2", L"i_2"};
        REQUIRE(indices.size() == 4);
        const ExprPtr expected = parse_expr(L"J{a1,a2;i1,i2}");

        for (const std::vector<std::size_t> &indexPerm :
             twoElectronIntegralSymmetries()) {
          REQUIRE(indexPerm.size() == 4);

          ExprPtr integralExpr = ex<Tensor>(
              L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
              ket{indices[indexPerm[2]], indices[indexPerm[3]]});

          auto transformed = integralExpr;
          remap_integrals(transformed);

          CAPTURE(indexPerm);
          CAPTURE_EXPR(integralExpr);
          CAPTURE_EXPR(transformed);
          CAPTURE_EXPR(expected);

          REQUIRE(transformed == expected);
        }
      }

      SECTION("virt,occ,virt,virt") {
        std::vector<Index> indices = {L"a_1", L"i_1", L"a_2", L"a_3"};
        REQUIRE(indices.size() == 4);
        const ExprPtr expected = parse_expr(L"J{a1,a2;a3,i1}");

        for (const std::vector<std::size_t> &indexPerm :
             twoElectronIntegralSymmetries()) {
          REQUIRE(indexPerm.size() == 4);

          ExprPtr integralExpr = ex<Tensor>(
              L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
              ket{indices[indexPerm[2]], indices[indexPerm[3]]});

          auto transformed = integralExpr;
          remap_integrals(transformed);

          CAPTURE(indexPerm);
          CAPTURE_EXPR(integralExpr);
          CAPTURE_EXPR(transformed);
          CAPTURE_EXPR(expected);

          REQUIRE(transformed == expected);
        }
      }
    }

    SECTION("f") {
      SECTION("same_space") {
        ExprPtr expr = parse_expr(L"f{i2;i1} + f{a1;a2}");
        const ExprPtr expected = parse_expr(L"f{i1;i2} + f{a1;a2}");

        remap_integrals(expr);

        CAPTURE_EXPR(expr);
        CAPTURE_EXPR(expected);

        REQUIRE(expr == expected);
      }
      SECTION("different_space") {
        ExprPtr expr = parse_expr(L"f{a1;i1} + f{i1;a1}");
        const ExprPtr expected = parse_expr(L"f{a1;i1} + f{a1;i1}");

        remap_integrals(expr);

        CAPTURE_EXPR(expr);
        CAPTURE_EXPR(expected);

        REQUIRE(expr == expected);
      }
    }
  }
}

#undef CAPTURE_EXPR
