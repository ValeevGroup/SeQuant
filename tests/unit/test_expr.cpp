//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"
#include "SeQuant/domain/mbpt/convention.hpp"

#include <iostream>
#include "SeQuant/core/hash.hpp"
#include "SeQuant/core/wick.hpp"

struct Dummy : public sequant::Expr {
  virtual ~Dummy() = default;
  std::wstring to_latex() const override { return L"{\\text{Dummy}}"; }
  std::wstring to_wolfram() const override { return L"Dummy[]"; }
  type_id_type type_id() const override { return get_type_id<Dummy>(); };
  sequant::ExprPtr clone() const override { return sequant::ex<Dummy>(); }
  bool static_equal(const sequant::Expr &that) const override { return true; }
};

template <typename T>
struct VecExpr : public std::vector<T>, public sequant::Expr {
  using base_type = std::vector<T>;
  using base_type::begin;
  using base_type::end;

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
  std::wstring to_wolfram() const override {
    std::wstring result = L"VecExpr[";
    size_t count = 1;
    for (const auto &e : *this) {
      const auto last_it = count == this->std::vector<T>::size();
      if constexpr (sequant::Expr::is_shared_ptr_of_expr_or_derived<T>::value) {
        result += e->to_wolfram() + (last_it ? L"" : L",");
      } else {
        result += std::to_wstring(e) + (last_it ? L"" : L",");
      }
      ++count;
    }
    result += L"]";
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
  std::wstring to_wolfram() const override {
    return L"Adjointable[" + std::to_wstring(v) + L"]";
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

TEST_CASE("Expr", "[elements]") {
  using namespace sequant;
  mbpt::set_default_convention();
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
    using ranges::size;

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
    REQUIRE(begin(*ex4) != end(*ex4));
    REQUIRE(size(*ex4) == 3);
    REQUIRE(begin(ex4->expr()) == end(ex4->expr()));
    REQUIRE(size(ex4->expr()) == 0);

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
    REQUIRE_THROWS_AS(ex->value<Dummy>(), std::invalid_argument);
    // sequant::rational is convertible to bool
    REQUIRE_NOTHROW(ex->value<bool>());
    REQUIRE_THROWS_AS(std::make_shared<Constant>(-2)->value<unsigned int>(),
                      std::range_error);
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
      REQUIRE_THROWS_AS(e->adjoint(), std::logic_error);
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
      REQUIRE(to_latex(e) ==
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

  SECTION("wolfram") {
    Product sp0{};
    sp0.append(2, std::make_shared<Dummy>());
    REQUIRE(to_wolfram(sp0) == L"Times[2,Dummy[]]");

    // VecExpr<ExprPtr>
    {
      const auto ex5_init = std::vector<std::shared_ptr<Constant>>{
          std::make_shared<Constant>(1), std::make_shared<Constant>(2),
          std::make_shared<Constant>(3)};
      auto ex6 =
          std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));
      REQUIRE(ex6->to_wolfram() == L"VecExpr[1,2,3]");
    }
  }

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
      //      std::wcout << "x = " << to_latex(x) << std::endl;
      REQUIRE(to_latex(x) ==
              L"{ \\bigl({{{3}}} + {{\\text{Dummy}}} + {{{3}}"
              L"{\\text{Dummy}}} + {{\\text{Dummy}}{\\text{Dummy}}}\\bigr) }");
      rapid_simplify(x);
      //      std::wcout << "x = " << to_latex(x) << std::endl;
    }
    {
      auto x =
          (ex<Constant>(1) +
           ex<Constant>(2) * (ex<Constant>(3) - ex<Dummy>())) *
          (ex<Constant>(5) * (ex<Constant>(6) + ex<Dummy>()) + ex<Dummy>());
      // std::wcout << "x = " << to_latex(x) << std::endl;
      REQUIRE(to_latex(x) ==
              L"{{ \\bigl({{{1}}} + {{{2}}{ \\bigl({{{3}}} - {"
              L"{\\text{Dummy}}}\\bigr) }}\\bigr) }{ \\bigl({{{5}}"
              L"{ \\bigl({{{6}}} + {\\text{Dummy}}\\bigr) }} + "
              L"{\\text{Dummy}}\\bigr) }}");
      expand(x);
      // std::wcout << "ex = " << to_latex(x) << std::endl;
      REQUIRE(to_latex(x) ==
              L"{ \\bigl({{{30}}} + {{{5}}{\\text{Dummy}}} + "
              L"{{\\text{Dummy}}} + {{{180}}} + {{{30}}"
              L"{\\text{Dummy}}} + {{{6}}{\\text{Dummy}}} - {{{60}}"
              L"{\\text{Dummy}}} - {{{10}}"
              L"{\\text{Dummy}}{\\text{Dummy}}} - {{{2}}"
              L"{\\text{Dummy}}{\\text{Dummy}}}\\bigr) }");
    }
  }

  SECTION("hashing") {
    const auto ex5_init = std::vector<std::shared_ptr<Constant>>{
        std::make_shared<Constant>(1), std::make_shared<Constant>(2),
        std::make_shared<Constant>(3)};
    REQUIRE_NOTHROW(hash_value(ex5_init));
    REQUIRE(hash_value(ex5_init) != hash_value(ex<Constant>(1)));

    auto hasher = [](const std::shared_ptr<const Expr> &) -> unsigned int {
      return 0;
    };
    REQUIRE_NOTHROW(ex<Constant>(1)->hash_value(hasher) == 0);
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

}  // TEST_CASE("Expr")

TEST_CASE("ExprPtr", "[elements]") {
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
    REQUIRE(ex1 == ex<Constant>(3));
    ex1 -= ex2;
    REQUIRE(ex1 == ex<Constant>(1));
    ex1 *= ex2;
    REQUIRE(ex1 == ex<Constant>(2));
  }
}  // TEST_CASE("ExprPtr")
