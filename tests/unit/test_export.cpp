#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/julia_itensor.hpp>
#include <SeQuant/core/export/julia_tensor_kit.hpp>
#include <SeQuant/core/export/julia_tensor_operations.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>

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
            "Create Z and initialize to zero\n"
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

    SECTION("binary + ternary") {
      auto tree = eval_node<EvalExpr>(
          parse_expr(L"A{a1;i1} B{i1;i2} + A{a1;i1} B{i1;i3} C{i3;i2}"));

      export_expression(tree, generator);

      std::string expected =
          "Declare index i_1\n"
          "Declare index i_2\n"
          "Declare index i_3\n"
          "Declare index a_1\n"
          "\n"
          "Declare tensor A[a_1, i_1]\n"
          "Declare tensor B[i_1, i_3]\n"
          "Declare tensor C[i_3, i_2]\n"
          "Declare tensor I[a_1, i_2]\n"
          "Declare tensor I2[a_1, i_3]\n"
          "\n"
          "Create I[a_1, i_2] and initialize to zero\n"
          "Create I2[a_1, i_3] and initialize to zero\n"
          "Load A[a_1, i_1]\n"
          "Load B[i_1, i_3]\n"
          "Compute I2[a_1, i_3] += A[a_1, i_1] B[i_1, i_3]\n"
          "Unload B[i_1, i_3]\n"
          "Unload A[a_1, i_1]\n"
          "Load C[i_3, i_2]\n"
          "Compute I[a_1, i_2] += I2[a_1, i_3] C[i_3, i_2]\n"
          "Unload C[i_3, i_2]\n"
          "Unload I2[a_1, i_3]\n"
          "Load A[a_1, i_1]\n"
          "Load B[i_1, i_2]\n"
          "Compute I[a_1, i_2] += A[a_1, i_1] B[i_1, i_2]\n"
          "Unload B[i_1, i_2]\n"
          "Unload A[a_1, i_1]\n"
          "Persist I[a_1, i_2]\n";

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
            "Declare tensor I[i_1, i_3, i_4, a_3]\n"
            "\n"
            "Create I[i_1, i_3, i_4, a_3] and initialize to zero\n"
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
            "Compute I[i_1, i_3, i_4, a_3] += I[i_1, a_3] I[i_3, i_4]\n"
            "Unload I[i_3, i_4]\n"
            "Unload I[i_1, a_3]\n"
            "Persist I[i_1, i_3, i_4, a_3]\n";

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
            "Declare tensor I[a_1, a_2, i_1, i_2]\n"
            "\n"
            "Create I[a_1, a_2, i_1, i_2] and initialize to zero\n"
            "Load A[a_1, i_1]\n"
            "Load B[a_2, i_2]\n"
            "Compute I[a_1, a_2, i_1, i_2] += A[a_1, i_1] B[a_2, i_2]\n"
            "Unload B[a_2, i_2]\n"
            "Unload A[a_1, i_1]\n"
            // TODO: there is some unload/load pairs that could be removed
            "Load A[a_1, i_1]\n"
            "Load B[a_2, i_2]\n"
            "Compute I[a_1, a_2, i_1, i_2] += A[a_1, i_1] B[a_2, i_2]\n"
            "Unload B[a_2, i_2]\n"
            "Unload A[a_1, i_1]\n"
            "Persist I[a_1, a_2, i_1, i_2]\n";

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
      SECTION("sum of index perms") {
        auto tree = eval_node<EvalExpr>(
            parse_expr(L"-1 g{a_2,i_3;a_3,i_1} * t{a_1,a_3;i_3,i_2} + "
                       L"2 g{a_1,i_3;i_1,a_3} * t{a_2,a_3;i_2,i_3}"));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index i_2\n"
            "Declare index i_3\n"
            "Declare index a_1\n"
            "Declare index a_2\n"
            "Declare index a_3\n"
            "\n"
            "Declare tensor I[a_1, a_2, i_1, i_2]\n"
            "Declare tensor g[a_1, i_3, i_1, a_3]\n"
            "Declare tensor g[a_2, i_3, a_3, i_1]\n"
            "Declare tensor t[a_1, a_3, i_3, i_2]\n"
            "\n"
            "Create I[a_1, a_2, i_1, i_2] and initialize to zero\n"
            "Load g[a_2, i_3, a_3, i_1]\n"
            "Load t[a_1, a_3, i_3, i_2]\n"
            "Compute I[a_1, a_2, i_1, i_2] += -1 g[a_2, i_3, a_3, i_1] "
            "t[a_1, a_3, i_3, i_2]\n"
            "Unload t[a_1, a_3, i_3, i_2]\n"
            "Unload g[a_2, i_3, a_3, i_1]\n"
            "Load g[a_1, i_3, i_1, a_3]\n"
            "Load t[a_2, a_3, i_2, i_3]\n"
            "Compute I[a_1, a_2, i_1, i_2] += 2 g[a_1, i_3, i_1, a_3] "
            "t[a_2, a_3, i_2, i_3]\n"
            "Unload t[a_2, a_3, i_2, i_3]\n"
            "Unload g[a_1, i_3, i_1, a_3]\n"
            "Persist I[a_1, a_2, i_1, i_2]\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
      SECTION("dot + dot") {
        auto tree = eval_node<EvalExpr>(
            parse_expr(L"2 g{i_1,i_2;a_1,a_2} * t{a_1,a_2;i_1,i_2} - "
                       L"g{i_1,i_2;a_1,a_2} * t{a_1,a_2;i_2,i_1}"));

        export_expression(tree, generator);

        std::string expected =
            "Declare index i_1\n"
            "Declare index i_2\n"
            "Declare index a_1\n"
            "Declare index a_2\n"
            "\n"
            "Declare variable Z\n"
            "\n"
            "Declare tensor g[i_1, i_2, a_1, a_2]\n"
            "Declare tensor t[a_1, a_2, i_1, i_2]\n"
            "\n"
            "Create Z and initialize to zero\n"
            "Load g[i_1, i_2, a_1, a_2]\n"
            "Load t[a_1, a_2, i_1, i_2]\n"
            "Compute Z += 2 g[i_1, i_2, a_1, a_2] t[a_1, a_2, i_1, i_2]\n"
            "Unload t[a_1, a_2, i_1, i_2]\n"
            "Unload g[i_1, i_2, a_1, a_2]\n"
            "Load g[i_1, i_2, a_1, a_2]\n"
            "Load t[a_1, a_2, i_2, i_1]\n"
            "Compute Z += -1 g[i_1, i_2, a_1, a_2] t[a_1, a_2, i_2, i_1]\n"
            "Unload t[a_1, a_2, i_2, i_1]\n"
            "Unload g[i_1, i_2, a_1, a_2]\n"
            "Persist Z\n";

        REQUIRE(generator.get_generated_code() == expected);
      }
    }
    SECTION("Intermediate with same shape as summation result") {
      auto tree = eval_node<EvalExpr>(
          parse_expr(L"A{a1;i1} B{i1;i2} + ( A{a1;i1} B{i1;i3} ) C{i3;i2}"));

      export_expression(tree, generator);

      std::string expected =
          "Declare index i_1\n"
          "Declare index i_2\n"
          "Declare index i_3\n"
          "Declare index a_1\n"
          "\n"
          "Declare tensor A[a_1, i_1]\n"
          "Declare tensor B[i_1, i_3]\n"
          "Declare tensor C[i_3, i_2]\n"
          "Declare tensor I[a_1, i_2]\n"
          "Declare tensor I2[a_1, i_3]\n"
          "\n"
          "Create I[a_1, i_2] and initialize to zero\n"
          "Create I2[a_1, i_3] and initialize to zero\n"
          "Load A[a_1, i_1]\n"
          "Load B[i_1, i_3]\n"
          "Compute I2[a_1, i_3] += A[a_1, i_1] B[i_1, i_3]\n"
          "Unload B[i_1, i_3]\n"
          "Unload A[a_1, i_1]\n"
          "Load C[i_3, i_2]\n"
          "Compute I[a_1, i_2] += I2[a_1, i_3] C[i_3, i_2]\n"
          "Unload C[i_3, i_2]\n"
          "Unload I2[a_1, i_3]\n"
          "Load A[a_1, i_1]\n"
          "Load B[i_1, i_2]\n"
          "Compute I[a_1, i_2] += A[a_1, i_1] B[i_1, i_2]\n"
          "Unload B[i_1, i_2]\n"
          "Unload A[a_1, i_1]\n"
          "Persist I[a_1, i_2]\n";

      REQUIRE(generator.get_generated_code() == expected);
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
        auto tree = eval_node<EvalExpr>(c * ex<Variable>(L"Dummy"));
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
        auto tree = eval_node<EvalExpr>(expr);
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
        auto tree = eval_node<EvalExpr>(expr);
        REQUIRE_THROWS_WITH(export_expression(tree, generator, ctx),
                            "Proto Indices are not (yet) supported!");
      }

      SECTION("binary contraction") {
        ExprPtr expr = parse_expr(L"g{i_1,i_2;a_1,a_2} * T2{a_1,a_2;i_3,i_4}");
        auto tree = eval_node<EvalExpr>(expr);

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
        auto tree = eval_node<EvalExpr>(expr);

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
        auto tree = eval_node<EvalExpr>(expr);
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
        auto tree = eval_node<EvalExpr>(
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
        auto tree = eval_node<EvalExpr>(expr);

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
        auto tree = eval_node<EvalExpr>(expr);

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
        auto tree = eval_node<EvalExpr>(expr);
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
        auto tree = eval_node<EvalExpr>(
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
        auto tree = eval_node<EvalExpr>(expr);

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
        auto tree = eval_node<EvalExpr>(expr);

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
        auto tree = eval_node<EvalExpr>(expr);
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
        auto tree = eval_node<EvalExpr>(
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

#undef CAPTURE_EXPR
