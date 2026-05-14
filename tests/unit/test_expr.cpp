//
// Created by Eduard Valeyev on 3/23/18.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/primitives.hpp>

struct Dummy : public sequant::Expr {
  virtual ~Dummy() = default;
  std::wstring to_latex() const override { return L"{\\text{Dummy}}"; }
  type_id_type type_id() const override { return get_type_id<Dummy>(); };
  sequant::ExprPtr clone() const override { return sequant::ex<Dummy>(); }
  bool static_equal(const sequant::Expr &) const override { return true; }
};

template <typename T>
struct VecExpr : public std::vector<T>, public sequant::Expr {
  using base_type = std::vector<T>;
  using base_type::begin;
  using base_type::end;
  using base_type::size;

  VecExpr() = default;
  template <typename U>
  VecExpr(std::initializer_list<U> elements) : std::vector<T>(elements) {}
  template <typename Iter>
  VecExpr(Iter begin, Iter end) : std::vector<T>(begin, end) {}
  virtual ~VecExpr() = default;
  std::wstring to_latex() const override {
    std::wstring result = L"{\\text{VecExpr}\\{";
    for (const auto &e : *this) {
      if constexpr (sequant::Expr::is_shared_ptr_of_expr_or_derived<T>::value) {
        result += e->to_latex() + L" ";
      } else {
        result += std::to_wstring(e) + L" ";
      }
    }
    result += L"\\}}";
    return result;
  }

  type_id_type type_id() const override { return get_type_id<VecExpr<T>>(); };

 private:
  cursor begin_cursor() const override {
    if constexpr (sequant::Expr::is_shared_ptr_of_expr<T>::value) {
      return base_type::empty() ? Expr::begin_cursor()
                                : cursor{&base_type::at(0)};
    } else {
      return Expr::begin_cursor();
    }
  };
  cursor end_cursor() const override {
    if constexpr (sequant::Expr::is_shared_ptr_of_expr<T>::value) {
      return base_type::empty() ? Expr::end_cursor()
                                : cursor{&base_type::at(0) + base_type::size()};
    } else {
      return Expr::end_cursor();
    }
  };
  cursor begin_cursor() override {
    return const_cast<const VecExpr &>(*this).begin_cursor();
  };
  cursor end_cursor() override {
    return const_cast<const VecExpr &>(*this).end_cursor();
  };

  bool static_equal(const sequant::Expr &that) const override {
    return static_cast<const base_type &>(*this) ==
           static_cast<const base_type &>(static_cast<const VecExpr &>(that));
  }

  sequant::ExprPtr clone() const override {
    return sequant::ex<VecExpr>(this->begin(), this->end());
  }
};

struct Adjointable : public sequant::Expr {
  Adjointable() = default;
  Adjointable(int v) : v(v) {}
  virtual ~Adjointable() = default;
  std::wstring to_latex() const override {
    return L"{\\text{Adjointable}{" + std::to_wstring(v) + L"}}";
  }
  type_id_type type_id() const override { return get_type_id<Adjointable>(); };
  sequant::ExprPtr clone() const override {
    return sequant::ex<Adjointable>(v);
  }
  bool static_equal(const sequant::Expr &that) const override {
    return v == that.as<Adjointable>().v;
  }
  void adjoint() override { v = -v; };

  int v = 1;
};

struct latex_visitor {
  void operator()(const std::shared_ptr<sequant::Expr> &expr) {
    result += expr->to_latex();
  }
  std::wstring result{};
};

TEST_CASE("expr", "[elements]") {
  using namespace sequant;
  SECTION("constructors") {
    REQUIRE_NOTHROW(std::make_shared<Constant>(2));
    const auto ex2 = std::make_shared<Constant>(2);
    REQUIRE_NOTHROW(std::make_shared<Variable>(L"q"));
    const auto ex_q = std::make_shared<Variable>(L"q");
    REQUIRE_NOTHROW(std::make_shared<VecExpr<double>>());
    const auto ex3 = std::make_shared<VecExpr<double>>();
    REQUIRE_NOTHROW(std::make_shared<VecExpr<double>>(
        std::initializer_list<double>{1.0, 2.0, 3.0}));
    const auto ex4 = std::make_shared<VecExpr<double>>(
        std::initializer_list<double>{1.0, 2.0, 3.0});
    REQUIRE_NOTHROW(std::make_shared<VecExpr<std::shared_ptr<Constant>>>(
        std::initializer_list<std::shared_ptr<Constant>>{
            std::make_shared<Constant>(1), std::make_shared<Constant>(2),
            std::make_shared<Constant>(3)}));
    const auto ex5 = std::make_shared<VecExpr<std::shared_ptr<Constant>>>(
        std::initializer_list<std::shared_ptr<Constant>>{
            std::make_shared<Constant>(1), std::make_shared<Constant>(2),
            std::make_shared<Constant>(3)});
    REQUIRE_NOTHROW(std::make_shared<Dummy>());
    const auto ex1 = std::make_shared<Dummy>();
  }

  SECTION("accessors") {
    {
      const auto ex = std::make_shared<Constant>(2);
      REQUIRE(ex->is_atom());
    }
    {
      const auto ex = std::make_shared<Variable>(L"q");
      REQUIRE(ex->is_atom());
    }
    {
      const auto ex = std::make_shared<VecExpr<double>>(
          std::initializer_list<double>{1.0, 2.0, 3.0});
      REQUIRE(ex->is_atom());
    }
    {
      const auto e = ex<Constant>(1) + ex<Constant>(2);
      REQUIRE(!e->is_atom());
    }
  }

  SECTION("comparison") {
    {
      const auto ex1 = std::make_shared<Constant>(1);
      const auto ex2 = std::make_shared<Constant>(2);
      const auto ex3 = std::make_shared<Constant>(1);
      const auto ex4 = std::make_shared<VecExpr<double>>();
      const auto ex5 = std::make_shared<VecExpr<ExprPtr>>(ExprPtrList{
          std::make_shared<Constant>(1), std::make_shared<Constant>(2),
          std::make_shared<Constant>(3)});
      const auto ex0 = std::make_shared<Dummy>();

      // type ids get assigned in the order of use, which is program dependent,
      // only check basic relations here
      REQUIRE(ex0->type_id() == Expr::get_type_id<Dummy>());
      REQUIRE(ex1->type_id() == Expr::get_type_id<Constant>());
      REQUIRE(ex4->type_id() == Expr::get_type_id<VecExpr<double>>());
      REQUIRE(ex4->type_id() <
              Expr::get_type_id<VecExpr<float>>());  // VecExpr<float> had not
                                                     // been used yet

      REQUIRE(*ex0 == *ex0);
      REQUIRE(*ex1 == *ex1);
      REQUIRE(*ex2 == *ex2);
      REQUIRE(*ex1 == *ex3);
      REQUIRE(*ex4 == *ex4);
      REQUIRE(*ex5 == *ex5);
      REQUIRE(*ex0 != *ex1);
    }
  }

  SECTION("iteration") {
    const auto ex1 = std::make_shared<Dummy>();
    REQUIRE(begin(*ex1) == end(*ex1));
    REQUIRE(size(*ex1) == 0);

    const auto ex2 = std::make_shared<Constant>(2);
    REQUIRE(begin(*ex2) == end(*ex2));
    REQUIRE(size(*ex2) == 0);

    const auto ex3 = std::make_shared<VecExpr<double>>();
    REQUIRE(begin(*ex3) == end(*ex3));
    REQUIRE(size(*ex3) == 0);
    REQUIRE(begin(ex3->expr()) == end(ex3->expr()));
    REQUIRE(size(ex3->expr()) == 0);

    const auto ex4 = std::make_shared<VecExpr<double>>(
        std::initializer_list<double>{1.0, 2.0, 3.0});
    REQUIRE(begin(*ex4) != end(*ex4));  // uses VecExpr::{begin,end}
    REQUIRE(size(*ex4) == 3);           // uses VecExpr::size
    REQUIRE(begin(ex4->expr()) == end(ex4->expr()));  // uses Expr::{begin,end}
    REQUIRE(size(ex4->expr()) == 0);                  // uses Expr::{begin,end}

    const auto ex5_init = std::vector<std::shared_ptr<Constant>>{
        std::make_shared<Constant>(1), std::make_shared<Constant>(2),
        std::make_shared<Constant>(3)};
    const auto ex5 = std::make_shared<VecExpr<std::shared_ptr<Constant>>>(
        begin(ex5_init), end(ex5_init));
    REQUIRE(begin(*ex5) != end(*ex5));
    REQUIRE(size(*ex5) == 3);
    REQUIRE(begin(ex5->expr()) == end(ex5->expr()));
    REQUIRE(size(ex5->expr()) == 0);

    {
      auto ex6 =
          std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));
      REQUIRE(begin(*ex6) != end(*ex6));
      REQUIRE(size(*ex6) == 3);
      REQUIRE(begin(ex6->expr()) != end(ex6->expr()));
      REQUIRE(size(ex6->expr()) == 3);
      const auto &front_ptr = *begin(ex6->expr());
      auto front_ptr_cast = std::dynamic_pointer_cast<Constant>(front_ptr);
      REQUIRE(front_ptr_cast);
      REQUIRE(front_ptr_cast->value() == 1);
      {
        *begin(ex6->expr()) = std::make_shared<Constant>(4);
        auto &front_ptr = *begin(ex6->expr());
        auto front_ptr_cast = std::dynamic_pointer_cast<Constant>(front_ptr);
        REQUIRE(front_ptr_cast);
        REQUIRE(front_ptr_cast->value() == 4);
      }
    }
  }

  SECTION("constant") {
    const auto ex = std::make_shared<Constant>(2);
    REQUIRE(ex->value() == 2);
    REQUIRE(ex->value<int>() == 2);
    REQUIRE(ex->value<std::complex<int>>() == std::complex<int>{2, 0});
    REQUIRE_THROWS_AS(ex->value<Dummy>(), Exception);
    // sequant::rational is convertible to bool
    REQUIRE_NOTHROW(ex->value<bool>());
    REQUIRE_THROWS_AS(std::make_shared<Constant>(-2)->value<unsigned int>(),
                      std::range_error);
  }

  SECTION("power") {
    const auto c2 = ex<Constant>(rational{1, 2});
    const auto vx = ex<Variable>(L"x");

    {  // constructors
      REQUIRE_NOTHROW(Power(c2, rational{1, 2}));
      REQUIRE_NOTHROW(Power(vx, rational{3, 1}));

      // convenience ctors:
      REQUIRE(Power(L"x", 2) == Power(vx, rational{2}));
      REQUIRE(Power(L"x", rational{1, 2}) == Power(vx, rational{1, 2}));
      REQUIRE(Power(2, 3) == Power(ex<Constant>(2), rational{3}));
      REQUIRE(Power(rational{2, 3}, 2) ==
              Power(ex<Constant>(rational{2, 3}), rational{2}));
      if constexpr (sequant::assert_behavior() ==
                    sequant::AssertBehavior::Throw) {
        // base must be a Constant or Variable; Power-of-Power is not allowed
        auto inner = ex<Power>(c2, rational{1, 2});
        REQUIRE_THROWS(Power(inner, rational{2, 3}));

        // 0^n is defined only for n >= 0
        REQUIRE_THROWS(Power(ex<Constant>(0), rational{-1}));
      }
    }

    {  // accessors
      Power p(c2, rational{1, 2});
      REQUIRE(p.base() == ex<Constant>(rational{1, 2}));
      REQUIRE(p.exponent() == rational{1, 2});

      // is_zero: base == 0, exponent > 0
      Power pz(ex<Constant>(0), rational{2});
      REQUIRE(pz.is_zero());
      // 0^0 is not zero by our convention
      Power pz2(ex<Constant>(0), rational{0});
      REQUIRE(!pz2.is_zero());
      REQUIRE(!p.is_zero());
    }

    {  // comparison
      Power p1(c2, rational{1, 2});
      Power p2(c2, rational{1, 2});
      REQUIRE(p1 == p2);

      Power p3(c2, rational{1, 3});
      REQUIRE(!(p1 == p3));

      Power p4(vx, rational{1, 2});
      REQUIRE(!(p1 == p4));

      // static_less_than: same base, compare by exponent
      Power plt_a(vx, rational{1, 2});
      Power plt_b(vx, rational{3, 4});
      REQUIRE(plt_a < plt_b);
      Power plt_c(vx, rational{1, 2});
      REQUIRE(!(plt_a < plt_c));
    }

    {  // operator*=
      // b^e1 *= b^e2 -> b^(e1+e2)
      Power pa(vx, rational{1, 2});
      Power pb(vx, rational{1, 3});
      pa *= pb;
      REQUIRE(pa.exponent() == rational{5, 6});  // 1/2 + 1/3

      // bare base: b^e *= b -> b^(e+1)
      Power pc(vx, rational{1, 2});
      pc *= vx.as<Variable>();
      REQUIRE(pc.exponent() == rational{3, 2});

      // (b^e)* *= b* -> ((b^(e+1))*) (Power-level conjugation; bare base
      // also conjugated)
      Power pc_conj(vx, rational{1, 2});
      pc_conj.conjugate();
      Variable vx_conj(L"x");
      vx_conj.conjugate();
      pc_conj *= vx_conj;
      REQUIRE(pc_conj.conjugated());
      REQUIRE(pc_conj.exponent() == rational{3, 2});

      // mirror case: (b*^e)* *= b -> ((b*^(e+1))*) (base Variable is conj'd,
      // Power conj'd, matches)
      Power pc_conj2(ex<Variable>(vx_conj), rational{1, 2});
      pc_conj2.conjugate();
      pc_conj2 *= vx.as<Variable>();
      REQUIRE(pc_conj2.exponent() == rational{3, 2});

      // 2^{1/2} * 2^{1/2} = 2
      Power pe(ex<Constant>(2), rational{1, 2});
      Power pf(ex<Constant>(2), rational{1, 2});
      pe *= pf;
      REQUIRE(pe.exponent() == rational{1});
      REQUIRE(to_latex(pe) == Constant(2).to_latex());
    }

    {  // Power should NOT be absorbed into Product::scalar_
      auto p = ex<Power>(vx, rational{1, 2});
      auto prod = ex<Product>(Product{});
      prod->as<Product>().append(1, p, Product::Flatten::Yes);
      REQUIRE(prod->as<Product>().factors().size() == 1);
      REQUIRE(prod->as<Product>().scalar() == 1);

      // simplify folds a foldable Power-in-Product into the Product scalar
      auto pwf = ex<Power>(2, 2) * ex<Variable>(L"x");
      simplify(pwf);
      REQUIRE(pwf->is<Product>());
      REQUIRE(pwf->as<Product>().scalar() == 4);

      // non-foldable Power stays as a factor; scalar is 1
      auto pwnf = ex<Power>(2, rational{1, 2}) * ex<Variable>(L"x");
      simplify(pwnf);
      REQUIRE(pwnf->is<Product>());
      REQUIRE(pwnf->as<Product>().scalar() == 1);
    }
  }

  SECTION("scaled_product") {
    REQUIRE_NOTHROW(Product{});
    Product sp0{};
    REQUIRE(sp0.scalar() == 1);
    REQUIRE(sp0.factors().empty());

    REQUIRE_NOTHROW(sp0.append(2, std::make_shared<Dummy>()));
    REQUIRE(sp0.scalar() == 2);
    REQUIRE(sp0.factors().size() == 1);
    REQUIRE(*(sp0.factors()[0]) == Dummy{});

    REQUIRE(begin(sp0.expr()) != end(sp0.expr()));
    REQUIRE(ranges::size(sp0.expr()) == 1);

    REQUIRE_NOTHROW(sp0.scale(2));
    REQUIRE(sp0.scalar() == 4);
  }

  SECTION("adjoint") {
    {  // not implemented by default
      const auto e = std::make_shared<Dummy>();
      REQUIRE_THROWS_AS(e->adjoint(), Exception);
    }
    {  // implemented in Adjointable
      const auto e = std::make_shared<Adjointable>();
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE_NOTHROW(adjoint(e));  // check free-function adjoint
    }
    {  // Constant
      const auto e = std::make_shared<Constant>(Constant::scalar_type{1, 2});
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->value() == Constant::scalar_type{1, -2});
    }
    {  // Variable
      const auto e = std::make_shared<Variable>(L"q");
      REQUIRE(e->conjugated() == false);
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->label() == L"q");
      REQUIRE(e->conjugated() == true);
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->conjugated() == false);
    }
    {  // Product
      const auto e = std::make_shared<Product>();
      e->append(Constant::scalar_type{2, -1}, ex<Adjointable>());
      e->append(1, ex<Adjointable>(-2));
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->scalar() == Constant::scalar_type{2, 1});
      REQUIRE(e->factors()[0]->as<Adjointable>().v == 2);
      REQUIRE(e->factors()[1]->as<Adjointable>().v == -1);
    }
    {  // CProduct
      const auto e = std::make_shared<CProduct>();
      e->append(Constant::scalar_type{2, -1}, ex<Adjointable>());
      e->append(1, ex<Adjointable>(-2));
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->scalar() == Constant::scalar_type{2, 1});
      REQUIRE(e->factors()[0]->as<Adjointable>().v == -1);
      REQUIRE(e->factors()[1]->as<Adjointable>().v == 2);
    }
    {  // NCProduct
      const auto e = std::make_shared<NCProduct>();
      e->append(Constant::scalar_type{2, -1}, ex<Adjointable>());
      e->append(1, ex<Adjointable>(-2));
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->scalar() == Constant::scalar_type{2, 1});
      REQUIRE(e->factors()[0]->as<Adjointable>().v == 2);
      REQUIRE(e->factors()[1]->as<Adjointable>().v == -1);
    }
    {  // Sum
      const auto e = std::make_shared<Sum>();
      e->append(ex<Adjointable>());
      e->append(ex<Adjointable>(-2));
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->summands()[0]->as<Adjointable>().v == -1);
      REQUIRE(e->summands()[1]->as<Adjointable>().v == 2);
    }
    {  // Power: adjoint flips the conjugation flag; base/exponent unchanged
      Power pv(ex<Variable>(L"z"), rational{1, 2});
      REQUIRE(!pv.conjugated());
      pv.adjoint();
      REQUIRE(pv.conjugated());
      REQUIRE(!pv.base()->as<Variable>().conjugated());
      REQUIRE(pv.exponent() == rational{1, 2});
      // double adjoint is identity
      pv.adjoint();
      REQUIRE(!pv.conjugated());

      using scalar_type = Constant::scalar_type;
      // i^{2} = -1
      auto i_sq = ex<Power>(ex<Constant>(scalar_type{0, 1}), 2);
      Power::flatten(i_sq);
      REQUIRE(i_sq->is<Constant>());
      REQUIRE(i_sq->as<Constant>().value() == scalar_type{-1});

      // (1+i)^{2} = 2i; ((1+i)^{2})* = -2i
      auto one_plus_i = ex<Constant>(scalar_type{1, 1});  // (1+i)
      auto square = ex<Power>(one_plus_i, 2);             // (1+i)^{2}
      square->as<Power>().adjoint();
      REQUIRE(square->as<Power>().conjugated());
      Power::flatten(square);
      REQUIRE(square->is<Constant>());
      REQUIRE(square->as<Constant>().value() ==
              scalar_type{0, -2});  // (2i)^{*} = -2i
    }
  }

  SECTION("clone") {
    {  // Variable
      const auto e = std::make_shared<Variable>(L"q");
      REQUIRE_NOTHROW(e->adjoint());
      const auto e_clone = e->clone();
      REQUIRE(e_clone.is<Variable>());
      REQUIRE(e->label() == e_clone.as<Variable>().label());
      REQUIRE(e->conjugated() == e_clone.as<Variable>().conjugated());
    }
  }  // SECTION("clone")

  SECTION("latex") {
    {  // Variable
      const auto e = std::make_shared<Variable>(L"q");
      REQUIRE(e->to_latex() == L"{q}");
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->to_latex() == L"{{q}^*}");
      REQUIRE_NOTHROW(e->adjoint());
      REQUIRE(e->to_latex() == L"{q}");
    }

    {  // Power
      const auto c2 = ex<Constant>(rational{1, 2});

      Power p1(c2, rational{1, 2});
      REQUIRE(to_latex(p1) == L"{\\frac{1}{2^{\\frac{1}{2}}}}");

      Power p_exp1(c2, rational{1});
      REQUIRE(to_latex(p_exp1) == L"{{{\\frac{1}{2}}}}");

      Power pv(ex<Variable>(L"x"), rational{2, 1});
      REQUIRE(to_latex(pv) == L"{x}^{2}");
    }

    Product sp0{};
    sp0.append(2, std::make_shared<Dummy>());
    REQUIRE(to_latex(sp0) == L"{{{2}}{\\text{Dummy}}}");

    // VecExpr<ExprPtr>
    {
      const auto ex5_init = std::vector<std::shared_ptr<Constant>>{
          std::make_shared<Constant>(1), std::make_shared<Constant>(2),
          std::make_shared<Constant>(3)};
      auto ex6 =
          std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));
      REQUIRE(ex6->to_latex() ==
              L"{\\text{VecExpr}\\{{{{1}}} {{{2}}} {{{3}}} \\}}");
    }

    // to_latex_align
    {
      const auto e = std::make_shared<Sum>();
      e->append(ex<Adjointable>(1));
      e->append(ex<Adjointable>(2));
      e->append(ex<Adjointable>(3));
      e->append(ex<Adjointable>(4));
      // std::wcout << "to_latex(e) = " << to_latex(e) << std::endl;
      REQUIRE(to_latex(*e) ==
              L"{ \\bigl({\\text{Adjointable}{1}} + {\\text{Adjointable}{2}} + "
              L"{\\text{Adjointable}{3}} + {\\text{Adjointable}{4}}\\bigr) }");
      // std::wcout << "to_latex_align(e) = " << to_latex_align(e) << std::endl;
      REQUIRE(to_latex_align(e) ==
              L"\\begin{align}\n"
              "& ({\\text{Adjointable}{1}} \\\\\n"
              "& + {\\text{Adjointable}{2}} \\\\\n"
              "& + {\\text{Adjointable}{3}} \\\\\n"
              "& + {\\text{Adjointable}{4}})\n"
              "\\end{align}");
      // std::wcout << "to_latex_align(e,5,2) = " << to_latex_align(e,5,2) <<
      // std::endl;
      REQUIRE(to_latex_align(e, 5, 2) ==
              L"\\begin{align}\n"
              "& ({\\text{Adjointable}{1}} + {\\text{Adjointable}{2}} \\\\\n"
              "& + {\\text{Adjointable}{3}} + {\\text{Adjointable}{4}})\n"
              "\\end{align}");
      // std::wcout << "to_latex_align(e,1,2) = " << to_latex_align(e,1,2) <<
      // std::endl;
      REQUIRE(to_latex_align(e, 1, 2) ==
              L"\\begin{align}\n"
              "& ({\\text{Adjointable}{1}} + {\\text{Adjointable}{2}} \\\\\n"
              "& + {\\text{Adjointable}{3}} \n"
              "\\end{align}\n"
              "\\begin{align}\n"
              "& + {\\text{Adjointable}{4}})\n"
              "\\end{align}");
    }
  }  // SECTION("latex")

  SECTION("visitor") {
    // read-only visitor
    {
      const auto ex5_init = std::vector<std::shared_ptr<Constant>>{
          std::make_shared<Constant>(1), std::make_shared<Constant>(2),
          std::make_shared<Constant>(3)};
      ExprPtr ex6 =
          std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));

      auto ex = ex6 + ex6;

      latex_visitor v1{};
      ex->visit(v1);

      //      std::wcout << "v1.result = " << v1.result << std::endl;
      REQUIRE(
          v1.result ==
          L"{{{1}}}{{{2}}}{{{3}}}{\\text{VecExpr}\\{{{{1}}} {{{2}}} {{{3}}} "
          L"\\}}{{{1}}}{{{2}}}{{{3}}}{\\text{VecExpr}\\{{{{1}}} {{{2}}} "
          L"{{{3}}} \\}}{ \\bigl({\\text{VecExpr}\\{{{{1}}} {{{2}}} {{{3}}} "
          L"\\}} + {\\text{VecExpr}\\{{{{1}}} {{{2}}} {{{3}}} \\}}\\bigr) }");

      latex_visitor v2{};
      ex->visit(v2, /* atoms_only = */ true);
      //      std::wcout << "v2.result = " << v2.result << std::endl;
      REQUIRE(v2.result ==
              L"{{{1}}}{{{2}}}{{{3}}}{{{1}}}{{{"
              L"2}}}{{{3}}}");
    }

    // mutating visitor
    {
      auto x = (ex<Constant>(1) + ex<Constant>(2)) *
               (ex<Constant>(3) + ex<Constant>(4));
      x->visit([](ExprPtr &expr) {
        if (expr->is<Constant>()) {
          expr = ex<Constant>(2 * expr->as<Constant>().value());
        }
      });
      CHECK_NOTHROW(simplify(x));
      CHECK(x->is<Constant>());
      CHECK(x->as<Constant>().value() == 84);
    }
  }

  SECTION("range") {
    {
      REQUIRE_NOTHROW(expr_range{});
      expr_range exrng{};
      REQUIRE(ranges::begin(exrng) == ranges::begin(exrng));
      REQUIRE(ranges::begin(exrng) == ranges::end(exrng));
    }

    // compares indices in address provided by cursor::address() to a list of
    // indices
    auto compare =
        [](const container::svector<std::pair<ExprPtr *, int64_t>> &address1,
           std::initializer_list<int> address2) {
          return address1.size() == address2.size() &&
                 std::equal(
                     begin(address1), end(address1), begin(address2),
                     [](const auto &parent_and_index1, const auto &index2) {
                       return parent_and_index1.second == index2;
                     });
        };

    {
      auto x = (ex<Constant>(1) + ex<Constant>(2)) *
               (ex<Constant>(3) + ex<Constant>(4));
      REQUIRE_NOTHROW(expr_range{x});
      expr_range exrng{x};
      REQUIRE(ranges::begin(exrng) == ranges::begin(exrng));
      REQUIRE(ranges::begin(exrng) != ranges::end(exrng));
      REQUIRE(std::distance(ranges::begin(exrng), ranges::end(exrng)) == 2);

      auto i = 0;
      for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
        switch (i) {
          case 0:
            REQUIRE(to_latex(*it) == L"{{{3}}}");
            REQUIRE(compare(ranges::get_cursor(it).address(), {0, 0}));
            REQUIRE(ranges::get_cursor(it).ordinal() == 0);
            break;
          case 1:
            REQUIRE(to_latex(*it) == L"{{{7}}}");
            REQUIRE(compare(ranges::get_cursor(it).address(), {1, 0}));
            REQUIRE(ranges::get_cursor(it).ordinal() == 1);
            break;
        }
        ++i;
      }
    }

    {
      auto x =
          (ex<Constant>(1) +
           ex<Constant>(2) * (ex<Constant>(3) - ex<Dummy>())) *
          (ex<Constant>(5) + (ex<Constant>(6) + ex<Dummy>()) * ex<Constant>(8));
      REQUIRE_NOTHROW(expr_range{x});
      expr_range exrng{x};
      REQUIRE(std::distance(ranges::begin(exrng), ranges::end(exrng)) == 6);

      auto i = 0;
      for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
        switch (i) {
          case 0:
            REQUIRE(to_latex(*it) == L"{{{1}}}");
            REQUIRE(compare(ranges::get_cursor(it).address(), {0, 0}));
            REQUIRE(ranges::get_cursor(it).ordinal() == 0);
            break;
          case 1:
            REQUIRE(to_latex(*it) == L"{{{3}}}");
            REQUIRE(compare(ranges::get_cursor(it).address(), {0, 1, 0, 0}));
            REQUIRE(ranges::get_cursor(it).ordinal() == 1);
            break;
          case 2:
            REQUIRE(to_latex(*it) == L"{\\text{Dummy}}");
            REQUIRE(compare(
                ranges::get_cursor(it).address(),
                {0, 1, 0, 1, 0}));  // -1 x Dummy = Product(-1){Dummy}!!!
            REQUIRE(ranges::get_cursor(it).ordinal() == 2);
            break;
          case 3:
            REQUIRE(to_latex(*it) == L"{{{5}}}");
            REQUIRE(compare(ranges::get_cursor(it).address(), {1, 0}));
            REQUIRE(ranges::get_cursor(it).ordinal() == 3);
            break;
          case 4:
            REQUIRE(to_latex(*it) == L"{{{6}}}");
            REQUIRE(compare(ranges::get_cursor(it).address(), {1, 1, 0, 0}));
            REQUIRE(ranges::get_cursor(it).ordinal() == 4);
            break;
          case 5:
            REQUIRE(to_latex(*it) == L"{\\text{Dummy}}");
            REQUIRE(compare(ranges::get_cursor(it).address(), {1, 1, 0, 1}));
            REQUIRE(ranges::get_cursor(it).ordinal() == 5);
            break;
        }
        ++i;
      }
    }
  }

  SECTION("expand") {
    {
      auto x =
          (ex<Constant>(1) + ex<Dummy>()) * (ex<Constant>(3) + ex<Dummy>());
      REQUIRE(to_latex(x) ==
              L"{{ \\bigl({{{1}}} + {\\text{Dummy}}\\bigr) }{ \\bigl({{{3}}} "
              L"+ {\\text{Dummy}}\\bigr) }}");
      expand(x);
      REQUIRE(to_latex(x) ==
              L"{ \\bigl({{{3}}} + {{\\text{Dummy}}} + {{{3}}"
              L"{\\text{Dummy}}} + {{\\text{Dummy}}{\\text{Dummy}}}\\bigr) }");
      rapid_simplify(x);
    }
    {
      auto x =
          (ex<Constant>(1) +
           ex<Constant>(2) * (ex<Constant>(3) - ex<Dummy>())) *
          (ex<Constant>(5) * (ex<Constant>(6) + ex<Dummy>()) + ex<Dummy>());
      REQUIRE(to_latex(x) ==
              L"{{ \\bigl({{{1}}} + {{{2}}{ \\bigl({{{3}}} - {"
              L"{\\text{Dummy}}}\\bigr) }}\\bigr) }{ \\bigl({{{5}}"
              L"{ \\bigl({{{6}}} + {\\text{Dummy}}\\bigr) }} + "
              L"{\\text{Dummy}}\\bigr) }}");
      expand(x);
      REQUIRE(to_latex(x) ==
              L"{ \\bigl({{{30}}} + {{{5}}{\\text{Dummy}}} + "
              L"{{\\text{Dummy}}} + {{{180}}} + {{{30}}"
              L"{\\text{Dummy}}} + {{{6}}{\\text{Dummy}}} - {{{60}}"
              L"{\\text{Dummy}}} - {{{10}}"
              L"{\\text{Dummy}}{\\text{Dummy}}} - {{{2}}"
              L"{\\text{Dummy}}{\\text{Dummy}}}\\bigr) }");
    }
  }

  SECTION("flatten") {
    {  // sums of sums
      auto x = ex<Constant>(1) + (ex<Dummy>() + ex<Constant>(3)) + ex<Dummy>();
      REQUIRE(to_latex(x) ==
              L"{ \\bigl({\\text{Dummy}} + {{{4}}} + {\\text{Dummy}}\\bigr) }");
      // make nested sums by visitation
      x->visit([](ExprPtr &e) {
        if (e.is<Dummy>()) {
          e = ex<Dummy>() + ex<Dummy>();
        }
      });
      REQUIRE(to_latex(x) ==
              L"{ \\bigl({ \\bigl({\\text{Dummy}} + {\\text{Dummy}}\\bigr) } + "
              L"{{{4}}} + { \\bigl({\\text{Dummy}} + {\\text{Dummy}}\\bigr) "
              L"}\\bigr) }");
      flatten(x);
      REQUIRE(to_latex(x) ==
              L"{ \\bigl({\\text{Dummy}} + {\\text{Dummy}} + {{{4}}} + "
              L"{\\text{Dummy}} + {\\text{Dummy}}\\bigr) }");
    }
    {  // products of products
      auto x = ex<Constant>(1) * (ex<Dummy>() * ex<Constant>(3)) * ex<Dummy>();
      REQUIRE(to_latex(x) == L"{{{3}}{\\text{Dummy}}{\\text{Dummy}}}");
      // make nested sums by visitation
      x->visit([](ExprPtr &e) {
        if (e.is<Dummy>()) {
          e = ex<Dummy>() * ex<Dummy>();
        }
      });
      REQUIRE(to_latex(x) ==
              L"{{{3}}\\bigl({{\\text{Dummy}}{\\text{Dummy}}}\\bigr)\\bigl({{"
              L"\\text{Dummy}}{\\text{Dummy}}}\\bigr)}");
      flatten(x);
      REQUIRE(to_latex(x) ==
              L"{{{3}}{\\text{Dummy}}{\\text{Dummy}}{\\text{Dummy}}{\\text{"
              L"Dummy}}}");
    }
    {  // products of products and sums
      auto x = ex<Constant>(1) * (ex<Dummy>() * ex<Constant>(3)) * ex<Dummy>();
      REQUIRE(to_latex(x) == L"{{{3}}{\\text{Dummy}}{\\text{Dummy}}}");
      // make nested sums by visitation
      x->visit([](ExprPtr &e) {
        if (e.is<Dummy>()) {
          e = ex<Dummy>() * (ex<Dummy>() + ex<Dummy>());
        }
      });
      REQUIRE(to_latex(x) ==
              L"{{{3}}\\bigl({{\\text{Dummy}}{ \\bigl({\\text{Dummy}} + "
              L"{\\text{Dummy}}\\bigr) }}\\bigr)\\bigl({{\\text{Dummy}}{ "
              L"\\bigl({\\text{Dummy}} + {\\text{Dummy}}\\bigr) }}\\bigr)}");
      flatten(x);
      REQUIRE(to_latex(x) ==
              L"{{{3}}{\\text{Dummy}}{ \\bigl({\\text{Dummy}} + "
              L"{\\text{Dummy}}\\bigr) }{\\text{Dummy}}{ "
              L"\\bigl({\\text{Dummy}} + {\\text{Dummy}}\\bigr) }}");
    }

    {  // Power::flatten

      // Constant base + integer exponent folds in place
      auto pf1 = ex<Power>(2, 3);
      Power::flatten(pf1);
      REQUIRE(pf1 == ex<Constant>(rational{8}));

      // non-integer exponent without an exact root is a no-op
      auto pf2 = ex<Power>(2, rational{1, 2});
      Power::flatten(pf2);
      REQUIRE(pf2->is<Power>());

      // Variable base is a no-op
      auto pf3 = ex<Power>(L"x", 2);
      Power::flatten(pf3);
      REQUIRE(pf3->is<Power>());

      // Variable base, zero exponent: x^0 = 1
      auto pf4 = ex<Power>(L"x", 0);
      Power::flatten(pf4);
      REQUIRE(pf4 == ex<Constant>(rational{1}));

      // 2^(-20) = 1/1048576
      auto pf5 = ex<Power>(2, -20);
      Power::flatten(pf5);
      REQUIRE(pf5 == ex<Constant>(rational{1, 1048576}));

      // square-root exponent with perfect-square base folds
      // 4^(1/2) = 2
      auto pf6 = ex<Power>(4, rational{1, 2});
      Power::flatten(pf6);
      REQUIRE(pf6 == ex<Constant>(rational{2}));

      // (1/4)^(1/2) = 1/2
      auto pf7 = ex<Power>(rational{1, 4}, rational{1, 2});
      Power::flatten(pf7);
      REQUIRE(pf7 == ex<Constant>(rational{1, 2}));

      // (1/4)^(-1/2) = 2
      auto pf8 = ex<Power>(rational{1, 4}, rational{-1, 2});
      Power::flatten(pf8);
      REQUIRE(pf8 == ex<Constant>(rational{2}));

      // (9/16)^(3/2) = 27/64
      auto pf9 = ex<Power>(rational{9, 16}, rational{3, 2});
      Power::flatten(pf9);
      REQUIRE(pf9 == ex<Constant>(rational{27, 64}));

      // (-1)^(1/2) is imaginary; left as Power
      auto pf10 = ex<Power>(-1, rational{1, 2});
      Power::flatten(pf10);
      REQUIRE(pf10->is<Power>());

      // negative-real base with half-integer exponent: not folded
      auto pf11 = ex<Power>(-4, rational{1, 2});
      Power::flatten(pf11);
      REQUIRE(pf11->is<Power>());

      // complex base (imag != 0) with half-integer exponent: not folded
      auto pf12 = ex<Power>(Constant::scalar_type{1, 1}, rational{1, 2});
      Power::flatten(pf12);
      REQUIRE(pf12->is<Power>());

      // 8^(1/3) is not folded
      auto pf13 = ex<Power>(8, rational{1, 3});
      Power::flatten(pf13);
      REQUIRE(pf13->is<Power>());

      // (b^1)* = conj(b):
      auto pf14 = ex<Power>(L"y", rational{1});
      pf14->as<Power>().conjugate();
      Power::flatten(pf14);
      REQUIRE(pf14->is<Variable>());
      REQUIRE(pf14->as<Variable>().label() == L"y");
      REQUIRE(pf14->as<Variable>().conjugated());
    }
  }

  SECTION("hashing") {
    const auto ex5_init = std::vector<std::shared_ptr<Constant>>{
        std::make_shared<Constant>(1), std::make_shared<Constant>(2),
        std::make_shared<Constant>(3)};
    REQUIRE_NOTHROW(hash_value(ex5_init));
    REQUIRE(hash_value(ex5_init) != hash_value(ex<Constant>(1)));

    REQUIRE(hash_value(ex<Constant>(1)) == hash_value(ex<Constant>(1)));

    auto hasher = [](const std::shared_ptr<const Expr> &) -> unsigned int {
      return 0;
    };
    REQUIRE_NOTHROW(ex<Constant>(1)->hash_value(hasher) == 0);

    {  // Power
      const auto c2 = ex<Constant>(rational{1, 2});
      const auto vx = ex<Variable>(L"x");

      Power p1(c2, rational{1, 2});
      Power p3(c2, rational{1, 3});  // different exponent
      Power p4(vx, rational{1, 2});  // different base
      REQUIRE(sequant::hash::value(p1) != sequant::hash::value(p3));
      REQUIRE(sequant::hash::value(p1) != sequant::hash::value(p4));

      // Power(b, 1) has the same hash as b
      const auto v = ex<Variable>(L"u");
      Power p_one(v, rational{1});
      REQUIRE(sequant::hash::value(p_one) == sequant::hash::value(*v));

      // mutating the ExprPtr does not affect the Power's cached hash.
      auto shared = ex<Variable>(L"s");
      Power ps(shared, rational{2});
      const auto ps_hash_before = sequant::hash::value(ps);
      shared->as<Variable>().conjugate();
      REQUIRE(sequant::hash::value(ps) == ps_hash_before);
    }
  }

  SECTION("commutativity") {
    const auto ex1 = std::make_shared<VecExpr<std::shared_ptr<Constant>>>(
        std::initializer_list<std::shared_ptr<Constant>>{
            std::make_shared<Constant>(1), std::make_shared<Constant>(2),
            std::make_shared<Constant>(3)});
    const auto ex2 =
        ex<Constant>(1) + (ex<Constant>(2) + ex<Constant>(3)) * ex<Constant>(4);

    REQUIRE(ex1->is_cnumber());
    REQUIRE(ex2->is_cnumber());
    REQUIRE(ex1->commutes_with(*ex1));
    REQUIRE(ex1->commutes_with(*ex2));
    REQUIRE(ex2->commutes_with(*ex1));
    REQUIRE(ex2->commutes_with(*ex2));
  }

  SECTION("expr_ptr") {
    using namespace sequant;

    SECTION("constructors") {
      REQUIRE_NOTHROW(ExprPtr{});
      REQUIRE_NOTHROW(ExprPtr{std::make_shared<Constant>(2)});
      const auto ex_two = std::make_shared<Constant>(2);
      REQUIRE_NOTHROW(ExprPtr{ex_two});
      ExprPtr ex;
      REQUIRE_NOTHROW(ex = std::make_shared<Constant>(2));
      REQUIRE_NOTHROW(ex = ex_two);
      ExprPtr ex2 = std::make_shared<Constant>(1);
      REQUIRE_NOTHROW(ExprPtr{ex2});
      REQUIRE_NOTHROW(ex = ex2);
    }

    SECTION("basic use") {
      ExprPtr ex1 = ex<Constant>(1);
      REQUIRE_NOTHROW(to_latex(ex1));
      REQUIRE_NOTHROW(ex1->to_latex());
    }

    SECTION("clone") {
      ExprPtr ex1 = ex<Constant>(1);
      ExprPtr ex2;
      REQUIRE_NOTHROW(ex2 = ex1.clone());
      CHECK(ex1);
      CHECK(ex1->as<Constant>().value() == ex2->as<Constant>().value());
      auto ex1_ptr = ex1.get();
      REQUIRE_NOTHROW(ex2 = std::move(ex1).clone());
      CHECK(!ex1);
      CHECK(ex2);
      CHECK(ex2.get() == ex1_ptr);
    }

    SECTION("iteration") {
      const auto ex1 = ex<Dummy>();
      REQUIRE(begin(*ex1) == end(*ex1));
      REQUIRE(size(*ex1) == 0);
      REQUIRE(begin(ex1) == end(ex1));
      REQUIRE(cbegin(ex1) == cend(ex1));
      REQUIRE(size(ex1) == 0);

      const auto ex2 = ex<Constant>(2);
      REQUIRE(begin(*ex2) == end(*ex2));
      REQUIRE(size(*ex2) == 0);
      REQUIRE(begin(ex2) == end(ex2));
      REQUIRE(size(ex2) == 0);

      const auto ex3 = ex<VecExpr<double>>();
      REQUIRE(begin(*ex3) == end(*ex3));
      REQUIRE(size(*ex3) == 0);
      REQUIRE(begin(ex3) == end(ex3));
      REQUIRE(size(ex3) == 0);
      REQUIRE(begin(ex3->expr()) == end(ex3->expr()));
      REQUIRE(size(ex3->expr()) == 0);

      const auto ex4 =
          ex<VecExpr<double>>(std::initializer_list<double>{1.0, 2.0, 3.0});
      CHECK(begin(*ex4) == end(*ex4));
      CHECK(size(*ex4) == 0);
      CHECK(begin(ex4) == end(ex4));
      CHECK(size(ex4) == 0);
      CHECK(begin(ex4->expr()) == end(ex4->expr()));
      CHECK(size(ex4->expr()) == 0);
      CHECK(begin(ex4.as<VecExpr<double>>()) != end(ex4.as<VecExpr<double>>()));
      CHECK(size(ex4.as<VecExpr<double>>()) == 3);
    }

    SECTION("operators") {
      ExprPtr ex1 = ex<Constant>(1);
      ExprPtr ex2 = ex<Constant>(2);
      REQUIRE_NOTHROW(ex1 + ex2);
      REQUIRE_NOTHROW(ex1 - ex2);
      REQUIRE_NOTHROW(ex1 * ex2);
      REQUIRE_NOTHROW(ex1 += ex2);
      REQUIRE_NOTHROW(ex1 -= ex2);
      REQUIRE_NOTHROW(ex1 *= ex2);
      ex1 = ex<Constant>(1);
      ex1 += ex2;
      CHECK(ex1 == ex<Constant>(3));
      ex1 -= ex2;
      CHECK(ex1 == ex<Constant>(1));
      ex1 *= ex2;
      CHECK(ex1 == ex<Constant>(2));

      ExprPtr ex3;
      REQUIRE_NOTHROW(ex3 += ex2);
      CHECK(ex3 == ex2);
      ex3.reset();
      REQUIRE_NOTHROW(ex3 -= ex2);
      CHECK(ex3 == ex<Constant>(-1) * ex2);
      ex3.reset();
      REQUIRE_NOTHROW(ex3 *= ex2);
      CHECK(ex3 == ex2);

      SECTION("Overloads with basic numeric types") {
        ex1 = ex<Constant>(1);

        ExprPtr res = ex1 + 1;
        simplify(res);
        REQUIRE(res == ex<Constant>(2));

        res = 1 + ex1;
        simplify(res);
        REQUIRE(res == ex<Constant>(2));

        res = ex1 - 1;
        simplify(res);
        REQUIRE(res == ex<Constant>(0));

        res = 1 - ex1;
        simplify(res);
        REQUIRE(res == ex<Constant>(0));

        res = ex1 * 5.0;
        simplify(res);
        REQUIRE(res == ex<Constant>(5));

        res = 5.0 * ex1;
        simplify(res);
        REQUIRE(res == ex<Constant>(5));

        // This will be rewritten as (1/5.0) * ex1
        res = ex1 / 5.0;
        simplify(res);
        REQUIRE(res == ex<Constant>(rational(1, 5)));
      }

      SECTION("Overloads with Variables") {
        auto One = ex<Constant>(1);
        auto Two = ex<Constant>(2);

        ExprPtr res1 = One + L"x";
        simplify(res1);
        REQUIRE(res1 == simplify(One + ex<Variable>(L"x")));
        REQUIRE(simplify(L"x" + One) == res1);

        ExprPtr res2 = res1 - "x";
        simplify(res2);
        REQUIRE(res2 == One);

        ExprPtr res3 = L"x" - One;
        simplify(res3);
        REQUIRE(res3 == ex<Variable>(L"x") - One);

        ExprPtr res4 = Two * "y";
        simplify(res4);
        REQUIRE(res4 == simplify(Two * ex<Variable>("y")));
        REQUIRE(simplify("y" * Two) == res4);
      }

      SECTION("Divide by Constant") {
        ex1 = ex<Constant>(5);

        ExprPtr res = ex1 / Constant(3);
        simplify(res);

        REQUIRE(res == ex<Constant>(rational(5, 3)));
      }
    }
  }

  SECTION("ResultExpr") {
    SECTION("accessors") {
      SECTION("as_variable") {
        REQUIRE_THAT(deserialize<ResultExpr>(L"R = Var").result_as_variable(),
                     EquivalentTo(L"R"));

        REQUIRE_THAT(deserialize<ResultExpr>(L"R = Var").result_as_tensor(),
                     EquivalentTo(L"R{}:N-N-N"));
      }
      SECTION("as_tensor") {
        REQUIRE_THAT(deserialize<ResultExpr>(L"R{a1;i2;p3} = T{a1;i2;p3}")
                         .result_as_tensor(),
                     EquivalentTo(L"R{a1;i2;p3}"));
      }
    }
    SECTION("particle pairings") {
      std::vector<std::pair<Index, Index>> expected = {{L"a_1", L"i_1"}};
      auto pairings = deserialize<ResultExpr>(L"R{a1;i1} = t{a1;i1}")
                          .index_particle_grouping<std::pair<Index, Index>>();
      REQUIRE_THAT(pairings, ::Catch::Matchers::UnorderedRangeEquals(expected));

      expected = {{L"a_1", L"i_1"}, {L"a_2", L"i_2"}};
      pairings = deserialize<ResultExpr>(L"R{a1,a2;i1,i2} = t{a1,a2;i1,i2}")
                     .index_particle_grouping<std::pair<Index, Index>>();
      REQUIRE_THAT(pairings, ::Catch::Matchers::UnorderedRangeEquals(expected));

      // aux indices are ignored
      expected = {{L"a_1", L"i_1"}, {L"a_2", L"i_2"}};
      pairings =
          deserialize<ResultExpr>(L"R{a1,a2;i1,i2;p1} = t{a1,a2;i1,i2;p1}")
              .index_particle_grouping<std::pair<Index, Index>>();
      REQUIRE_THAT(pairings, ::Catch::Matchers::UnorderedRangeEquals(expected));
    }
  }
}
