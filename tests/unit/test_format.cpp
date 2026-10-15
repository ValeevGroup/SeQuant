#include <SeQuant/core/io/format.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/string.hpp>

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
  // stateless; operator== has already matched the type id
  bool static_equal(const sequant::Expr&) const override { return true; }
};

}  // namespace

#if defined(__cpp_lib_format_ranges)
// expressions are ranges; they must not be formatted as such
static_assert(std::formattable<sequant::Expr, char>);
static_assert(std::formattable<sequant::Sum, char>);
static_assert(std::formattable<sequant::Product, wchar_t>);
static_assert(std::formattable<sequant::Variable, wchar_t>);
static_assert(std::formattable<LatexOnlyExpr, char>);
#endif

TEST_CASE("expression formatting defaults", "[format]") {
  using namespace sequant;

  const Variable variable(L"x");
  const Expr& expr = variable;
  const auto ptr = ex<Power>(L"x", 2);
  const ExprPtr null;

  REQUIRE(std::format("{}", variable) == toUtf8(serialize(expr)));
  REQUIRE(std::format("{}", expr) == toUtf8(serialize(expr)));
  REQUIRE(std::format("{}", ptr) == toUtf8(serialize(ptr)));
  REQUIRE(std::format("{}", null).empty());
  REQUIRE(std::format(L"{}", variable) == serialize(expr));
  REQUIRE(std::format(L"{}", expr) == serialize(expr));
  REQUIRE(std::format(L"{}", ptr) == serialize(ptr));
  REQUIRE(std::format(L"{}", null).empty());
}

TEST_CASE("tensor formatting", "[format]") {
  using namespace sequant;

  const auto t = ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"});
  REQUIRE(std::format("{}", t) == toUtf8(serialize(t)));
  REQUIRE(std::format("{:l}", t) == toUtf8(to_latex(t)));
  std::wostringstream wide;
  wide << t;
  REQUIRE(wide.str() == serialize(t));
}

TEST_CASE("result expression formatting", "[format]") {
  using namespace sequant;

  const Tensor r(L"R", bra{L"a_1"}, ket{L"i_1"});
  const ResultExpr result(r, ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}));
  REQUIRE(std::format("{}", result) == toUtf8(serialize(result)));
  REQUIRE(std::format(L"{:l}", result) == to_latex(result));
  std::wostringstream wide;
  wide << result;
  REQUIRE(wide.str() == serialize(result));
}

TEST_CASE("expression stream insertion", "[format]") {
  using namespace sequant;

  const Variable variable(L"Ж");
  const Expr& expr = variable;
  const auto ptr = ex<Power>(L"x", 2);
  const ExprPtr null;
  const auto var_text = serialize(expr);
  const auto ptr_text = serialize(ptr);
  std::ostringstream narrow;
  narrow << variable << '|' << expr << '|' << ptr << '|' << null;
  REQUIRE(narrow.str() ==
          toUtf8(var_text + L"|" + var_text + L"|" + ptr_text + L"|"));

  std::wostringstream wide;
  wide << variable << L'|' << expr << L'|' << ptr << L'|' << null;
  REQUIRE(wide.str() == var_text + L"|" + var_text + L"|" + ptr_text + L"|");
  REQUIRE(std::format("{}", variable) == toUtf8(var_text));
  REQUIRE(std::format(L"{}", variable) == var_text);
  REQUIRE(std::format("{:s}", variable) == toUtf8(var_text));
  REQUIRE(std::format(L"{:s}", variable) == var_text);
}

TEST_CASE("expression format selectors", "[format]") {
  using namespace sequant;

  const Variable variable(L"x");
  const Expr& expr = variable;
  const auto ptr = ex<Power>(L"x", 2);
  const ExprPtr null;

  const auto var_latex = toUtf8(to_latex(expr));
  const auto ptr_text = toUtf8(serialize(ptr));
  REQUIRE(std::format("{:l}|{:latex}", variable, expr) ==
          var_latex + "|" + var_latex);
  REQUIRE(std::format("{:s}|{:serialize}", ptr, ptr) ==
          ptr_text + "|" + ptr_text);
  REQUIRE(std::format(L"{:l}|{:latex}", ptr, ptr) ==
          to_latex(ptr) + L"|" + to_latex(ptr));
  REQUIRE(std::format(L"{:s}|{:serialize}", variable, expr) ==
          serialize(expr) + L"|" + serialize(expr));
  REQUIRE(std::format("{:l}|{:latex}", null, null) == "|");
  REQUIRE(std::format(L"{:l}|{:latex}", null, null) == L"|");
  REQUIRE(std::format("{:s}|{:serialize}", null, null) == "|");
  REQUIRE(std::format(L"{:s}|{:serialize}", null, null) == L"|");
  REQUIRE(io::serialization::from_string<ExprPtr>(std::format("{:s}", null)) ==
          nullptr);
  REQUIRE(io::serialization::from_string<ExprPtr>(
              std::format(L"{:serialize}", null)) == nullptr);

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
  REQUIRE(std::format("{:l}|{:latex}|{:l}", custom, expr, ptr) ==
          "\\mathcal{X}|\\mathcal{X}|\\mathcal{X}");
  REQUIRE(std::format(L"{:latex}", ptr) == L"\\mathcal{X}");
  REQUIRE_THROWS_AS(std::format("{}", ptr), sequant::Exception);
  REQUIRE_THROWS_AS(std::format(L"{}", custom), sequant::Exception);
}
