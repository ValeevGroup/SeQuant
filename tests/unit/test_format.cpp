#include <SeQuant/core/io/format.hpp>

#include <SeQuant/core/expr.hpp>

#include <catch2/catch_test_macros.hpp>

#include <format>
#include <initializer_list>
#include <sstream>
#include <string>

namespace {

struct LatexOnlyExpr : sequant::Expr {
  std::wstring to_latex() const override { return L"\\mathcal{X}"; }
  type_id_type type_id() const override { return get_type_id<LatexOnlyExpr>(); }
  sequant::ExprPtr clone() const override {
    return sequant::ex<LatexOnlyExpr>();
  }
  void adjoint() override {}
  bool static_equal(const sequant::Expr&) const override { return true; }
};

}  // namespace

TEST_CASE("expression formatting defaults", "[format]") {
  using namespace sequant;

  const Variable variable(L"x");
  const Expr& expr = variable;
  const auto ptr = ex<Power>(L"x", 2);
  const ExprPtr null;

  REQUIRE(std::format("{}", variable) == "{x}");
  REQUIRE(std::format("{}", expr) == "{x}");
  REQUIRE(std::format("{}", ptr) == "{x}^{2}");
  REQUIRE(std::format("{}", null) == "NULL");
  REQUIRE(std::format(L"{}", variable) == L"{x}");
  REQUIRE(std::format(L"{}", expr) == L"{x}");
  REQUIRE(std::format(L"{}", ptr) == L"{x}^{2}");
  REQUIRE(std::format(L"{}", null) == L"NULL");
}

TEST_CASE("expression stream insertion", "[format]") {
  using namespace sequant;

  const Variable variable(L"Ж");
  const Expr& expr = variable;
  const auto ptr = ex<Power>(L"x", 2);
  const ExprPtr null;
  std::ostringstream narrow;
  narrow << variable << '|' << expr << '|' << ptr << '|' << null;
  REQUIRE(narrow.str() == "{Ж}|{Ж}|{x}^{2}|NULL");

  std::wostringstream wide;
  wide << variable << L'|' << expr << L'|' << ptr << L'|' << null;
  REQUIRE(wide.str() == L"{Ж}|{Ж}|{x}^{2}|NULL");
  REQUIRE(std::format("{}", variable) == "{Ж}");
  REQUIRE(std::format(L"{}", variable) == L"{Ж}");
  REQUIRE(std::format("{:s}", variable) == "Ж");
  REQUIRE(std::format(L"{:s}", variable) == L"Ж");
}

TEST_CASE("expression format selectors", "[format]") {
  using namespace sequant;

  const Variable variable(L"x");
  const Expr& expr = variable;
  const auto ptr = ex<Power>(L"x", 2);
  const ExprPtr null;

  REQUIRE(std::format("{:l}|{:latex}", variable, expr) == "{x}|{x}");
  REQUIRE(std::format("{:s}|{:serialize}", ptr, ptr) == "x^(2)|x^(2)");
  REQUIRE(std::format(L"{:l}|{:latex}", ptr, ptr) == L"{x}^{2}|{x}^{2}");
  REQUIRE(std::format(L"{:s}|{:serialize}", variable, expr) == L"x|x");
  REQUIRE(std::format("{:s}|{:serialize}", null, null) == "NULL|NULL");
  REQUIRE(std::format(L"{:s}|{:serialize}", null, null) == L"NULL|NULL");

  for (const auto spec : {"{:bogus}", "{:latexx}", "{:ls}", "{:>20}"}) {
    REQUIRE_THROWS_AS(std::vformat(spec, std::make_format_args(ptr)),
                      Exception);
  }
  REQUIRE_THROWS_AS(std::vformat(L"{:bogus}", std::make_wformat_args(ptr)),
                    Exception);
}

TEST_CASE("format custom expressions", "[format]") {
  const LatexOnlyExpr custom;
  const sequant::Expr& expr = custom;
  const auto ptr = custom.clone();
  REQUIRE(std::format("{}|{}|{}", custom, expr, ptr) ==
          "\\mathcal{X}|\\mathcal{X}|\\mathcal{X}");
  REQUIRE(std::format(L"{}", ptr) == L"\\mathcal{X}");
  REQUIRE_THROWS_AS(std::format("{:s}", ptr), sequant::Exception);
  REQUIRE_THROWS_AS(std::format(L"{:serialize}", custom), sequant::Exception);
}
