#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/parse.hpp>

#include <codecvt>
#include <locale>
#include <vector>

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

TEST_CASE("Itf export", "[exports]") {
  using namespace sequant;
  ItfContext ctx;

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
        const ExprPtr expected = parse_expr(L"K{i1,i2;i3,i4}");

        for (const std::vector<std::size_t> &indexPerm :
             twoElectronIntegralSymmetries()) {
          REQUIRE(indexPerm.size() == 4);

          ExprPtr integralExpr = ex<Tensor>(
              L"g", bra{indices[indexPerm[0]], indices[indexPerm[1]]},
              ket{indices[indexPerm[2]], indices[indexPerm[3]]});

          auto transformed = integralExpr;
          remap_integrals(transformed, ctx);

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
          remap_integrals(transformed, ctx);

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
          remap_integrals(transformed, ctx);

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
          remap_integrals(transformed, ctx);

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
          remap_integrals(transformed, ctx);

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

        remap_integrals(expr, ctx);

        CAPTURE_EXPR(expr);
        CAPTURE_EXPR(expected);

        REQUIRE(expr == expected);
      }
      SECTION("different_space") {
        ExprPtr expr = parse_expr(L"f{a1;i1} + f{i1;a1}");
        const ExprPtr expected = parse_expr(L"f{a1;i1} + f{a1;i1}");

        remap_integrals(expr, ctx);

        CAPTURE_EXPR(expr);
        CAPTURE_EXPR(expected);

        REQUIRE(expr == expected);
      }
    }
  }
}
