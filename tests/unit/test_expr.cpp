//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <iostream>
#include "../../src/SeQuant2/wick.hpp"

TEST_CASE("Expr", "[elements]") {

  using namespace sequant2;

  struct Dummy : public Expr {
    virtual ~Dummy() = default;
    bool operator==(const Expr &expr) const override {
      const auto &expr_cast_ptr = dynamic_cast<const Dummy &>(expr);
      return true;
    }
    std::wstring to_latex() const override {
      return L"{\\text{Dummy}}";
    }
  };

  SECTION("expr") {
    REQUIRE_NOTHROW(Expr{});
    Expr ex0{};
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
    REQUIRE(sp0.to_latex() == L"{2.000000 \\times {\\text{Dummy}}}");
  }

}