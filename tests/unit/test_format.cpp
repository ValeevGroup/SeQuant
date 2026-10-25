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
  std::string to_latex() const override { return "\\mathcal{X}"; }
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
static_assert(std::formattable<sequant::Product, char>);
static_assert(std::formattable<sequant::Variable, char>);
static_assert(std::formattable<LatexOnlyExpr, char>);
#endif

TEST_CASE("expression formatting defaults", "[format]") {
  using namespace sequant;

  const Variable variable("x");
  const Expr& expr = variable;
  const auto ptr = ex<Power>("x", 2);
  const ExprPtr null;

  REQUIRE(std::format("{}", variable) == serialize(expr));
  REQUIRE(std::format("{}", expr) == serialize(expr));
  REQUIRE(std::format("{}", ptr) == serialize(ptr));
  REQUIRE(std::format("{}", null).empty());
}

TEST_CASE("tensor formatting", "[format]") {
  using namespace sequant;

  const auto t = ex<Tensor>("t", bra{"a_1"}, ket{"i_1"});
  REQUIRE(std::format("{}", t) == serialize(t));
  REQUIRE(std::format("{:l}", t) == to_latex(t));
  std::ostringstream stream;
  stream << t;
  REQUIRE(stream.str() == serialize(t));
}

TEST_CASE("result expression formatting", "[format]") {
  using namespace sequant;

  const Tensor r("R", bra{"a_1"}, ket{"i_1"});
  const ResultExpr result(r, ex<Tensor>("t", bra{"a_1"}, ket{"i_1"}));
  REQUIRE(std::format("{}", result) == serialize(result));
  REQUIRE(std::format("{:l}", result) == to_latex(result));
  std::ostringstream stream;
  stream << result;
  REQUIRE(stream.str() == serialize(result));
}

TEST_CASE("expression stream insertion", "[format]") {
  using namespace sequant;

  const Variable variable("Ж");
  const Expr& expr = variable;
  const auto ptr = ex<Power>("x", 2);
  const ExprPtr null;
  const auto var_text = serialize(expr);
  const auto ptr_text = serialize(ptr);
  std::ostringstream stream;
  stream << variable << '|' << expr << '|' << ptr << '|' << null;
  REQUIRE(stream.str() == var_text + "|" + var_text + "|" + ptr_text + "|");
  REQUIRE(std::format("{}", variable) == var_text);
  REQUIRE(std::format("{:s}", variable) == var_text);
}

TEST_CASE("expression format selectors", "[format]") {
  using namespace sequant;

  const Variable variable("x");
  const Expr& expr = variable;
  const auto ptr = ex<Power>("x", 2);
  const ExprPtr null;

  const auto var_latex = to_latex(expr);
  const auto ptr_text = serialize(ptr);
  REQUIRE(std::format("{:l}|{:latex}", variable, expr) ==
          var_latex + "|" + var_latex);
  REQUIRE(std::format("{:s}|{:serialize}", ptr, ptr) ==
          ptr_text + "|" + ptr_text);
  REQUIRE(std::format("{:l}|{:latex}", ptr, ptr) ==
          to_latex(ptr) + "|" + to_latex(ptr));
  REQUIRE(std::format("{:s}|{:serialize}", variable, expr) ==
          serialize(expr) + "|" + serialize(expr));
  REQUIRE(std::format("{:l}|{:latex}", null, null) == "|");
  REQUIRE(std::format("{:s}|{:serialize}", null, null) == "|");
  REQUIRE(io::serialization::from_string<ExprPtr>(std::format("{:s}", null)) ==
          nullptr);
  REQUIRE(io::serialization::from_string<ExprPtr>(
              std::format("{:serialize}", null)) == nullptr);

  for (const auto spec : {"{:bogus}", "{:latexx}", "{:ls}", "{:>20}"}) {
    REQUIRE_THROWS_AS(std::vformat(spec, std::make_format_args(ptr)),
                      Exception);
  }
}

TEST_CASE("format custom expressions", "[format]") {
  const LatexOnlyExpr custom;
  const sequant::Expr& expr = custom;
  const auto ptr = custom.clone();
  REQUIRE(std::format("{:l}|{:latex}|{:l}", custom, expr, ptr) ==
          "\\mathcal{X}|\\mathcal{X}|\\mathcal{X}");
  REQUIRE(std::format("{:latex}", ptr) == "\\mathcal{X}");
  REQUIRE_THROWS_AS(std::format("{}", ptr), sequant::Exception);
  REQUIRE_THROWS_AS(std::format("{}", custom), sequant::Exception);
}
