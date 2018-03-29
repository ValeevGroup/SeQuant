//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <iostream>
#include "../../src/SeQuant2/wick.hpp"

struct Dummy : public sequant2::Expr {
  virtual ~Dummy() = default;
  bool operator==(const sequant2::Expr &expr) const override {
    const auto &expr_cast_ptr = dynamic_cast<const Dummy &>(expr);
    return true;
  }
  std::wstring to_latex() const override {
    return L"{\\text{Dummy}}";
  }
};

template <typename T>
struct VecExpr : public std::vector<T>, public sequant2::Expr {
  using base_type = std::vector<T>;
  using base_type::begin;
  using base_type::end;

  VecExpr() = default;
  template <typename U> VecExpr(std::initializer_list<U> elements) : std::vector<T>(elements) {}
  virtual ~VecExpr() = default;
  bool operator==(const sequant2::Expr &expr) const override {
    const auto &expr_cast_ptr = dynamic_cast<const VecExpr &>(expr);
    return true;
  }
  std::wstring to_latex() const override {
    std::wstring result = L"{\\text{VecExpr}\\{";
    for(const auto& e: *this) {
      result += std::to_wstring(e) + L" ";
    }
    result += L"\\}}";
    return result;
  }
};

TEST_CASE("Expr", "[elements]") {

  using namespace sequant2;

  SECTION("expr") {
    REQUIRE_NOTHROW(Expr{});
    Expr ex0{};
    REQUIRE_NOTHROW(std::make_shared<Dummy>());
    const auto ex1 = std::make_shared<Dummy>();
    REQUIRE_NOTHROW(std::make_shared<VecExpr<double>>());
    const auto ex2 = std::make_shared<VecExpr<double>>();
    REQUIRE_NOTHROW(std::make_shared<VecExpr<double>>(std::initializer_list<double>{1.0, 2.0, 3.0}));
    const auto ex3 = std::make_shared<VecExpr<double>>(std::initializer_list<double>{1.0, 2.0, 3.0});
  }

  SECTION("iteration") {
    Expr ex0{};
    using ranges::size;
    REQUIRE(begin(ex0) == end(ex0));
    REQUIRE(size(ex0) == 0);

    const auto ex1 = std::make_shared<Dummy>();
    REQUIRE(begin(*ex1) == end(*ex1));
    REQUIRE(size(*ex1) == 0);

    const auto ex2 = std::make_shared<VecExpr<double>>(std::initializer_list<double>{1.0, 2.0, 3.0});
    REQUIRE(begin(*ex2) != end(*ex2));
    REQUIRE(size(*ex2) == 3);
    REQUIRE(begin(ex2->expr()) == end(ex2->expr()));
    REQUIRE(size(ex2->expr()) == 0);
  }

  SECTION("scaled_product") {
    REQUIRE_NOTHROW(ScaledProduct{});
    ScaledProduct sp0{};
    REQUIRE(sp0.scalar() == 1.0);
    REQUIRE(sp0.factors().empty());

    REQUIRE_NOTHROW(sp0.append(2.0, std::make_shared<Dummy>()));
    REQUIRE(sp0.scalar() == 2.0);
    REQUIRE(sp0.factors().size() == 1);
    REQUIRE(*(sp0.factors()[0]) == Dummy{});
  }

  SECTION("latex") {
    ScaledProduct sp0{};
    sp0.append(2.0, std::make_shared<Dummy>());
    REQUIRE(to_latex(sp0) == L"{2.000000 \\times {\\text{Dummy}}}");
  }

}