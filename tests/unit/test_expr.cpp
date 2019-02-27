//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <iostream>
#include "../../src/SeQuant/wick.hpp"
#include "../../src/SeQuant/hash.hpp"

struct Dummy : public sequant::Expr {
  virtual ~Dummy() = default;
  std::wstring to_latex() const override {
    return L"{\\text{Dummy}}";
  }
  std::wstring to_wolfram() const override {
    return L"Dummy[]";
  }
  type_id_type type_id() const override { return get_type_id<Dummy>(); };
  sequant::ExprPtr clone() const override { return sequant::ex<Dummy>(); }
  bool static_equal(const sequant::Expr &that) const override { return true; }
};

template<typename T>
struct VecExpr : public std::vector<T>, public sequant::Expr {
  using base_type = std::vector<T>;
  using base_type::begin;
  using base_type::end;

  VecExpr() = default;
  template<typename U>
  VecExpr(std::initializer_list<U> elements) : std::vector<T>(elements) {}
  template<typename Iter>
  VecExpr(Iter begin, Iter end) : std::vector<T>(begin, end) {}
  virtual ~VecExpr() = default;
  std::wstring to_latex() const override {
    std::wstring result = L"{\\text{VecExpr}\\{";
    for (const auto &e: *this) {
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
    for (const auto &e: *this) {
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

  type_id_type type_id() const override{
    return get_type_id<VecExpr<T>>();
  };

 private:
  cursor begin_cursor() const override {
    if constexpr (sequant::Expr::is_shared_ptr_of_expr<T>::value) {
      return base_type::empty() ? Expr::begin_cursor() : cursor{&base_type::at(0)};
    } else {
      return Expr::begin_cursor();
    }
  };
  cursor end_cursor() const override {
    if constexpr (sequant::Expr::is_shared_ptr_of_expr<T>::value) {
      return base_type::empty() ? Expr::end_cursor() : cursor{&base_type::at(0) + base_type::size()};
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
    return static_cast<const base_type&>(*this) == static_cast<const base_type&>(static_cast<const VecExpr&>(that));
  }

};

struct latex_visitor {
  void operator()(const std::shared_ptr<sequant::Expr>& expr) {
    result += expr->to_latex();
  }
  std::wstring result {};
};

TEST_CASE("Expr", "[elements]") {

  using namespace sequant;

  SECTION("constructors") {
    REQUIRE_NOTHROW(std::make_shared<Constant>(2));
    const auto ex2 = std::make_shared<Constant>(2);
    REQUIRE_NOTHROW(std::make_shared<VecExpr<double>>());
    const auto ex3 = std::make_shared<VecExpr<double>>();
    REQUIRE_NOTHROW(std::make_shared<VecExpr<double>>(std::initializer_list<double>{1.0, 2.0, 3.0}));
    const auto ex4 = std::make_shared<VecExpr<double>>(std::initializer_list<double>{1.0, 2.0, 3.0});
    REQUIRE_NOTHROW(std::make_shared<VecExpr<std::shared_ptr<Constant>>>(std::initializer_list<std::shared_ptr<Constant>>{
        std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0), std::make_shared<Constant>(3.0)}));
    const auto ex5 =
        std::make_shared<VecExpr<std::shared_ptr<Constant>>>(std::initializer_list<std::shared_ptr<Constant>>{
            std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0), std::make_shared<Constant>(3.0)});
    REQUIRE_NOTHROW(std::make_shared<Dummy>());
    const auto ex1 = std::make_shared<Dummy>();
  }

  SECTION("accessors") {
    {
      const auto ex = std::make_shared<Constant>(2);
      REQUIRE(ex->is_atom());
    }
    {
      const auto ex = std::make_shared<VecExpr<double>>(std::initializer_list<double>{1.0, 2.0, 3.0});
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
      const auto ex5 =
          std::make_shared<VecExpr<ExprPtr>>(ExprPtrList{
              std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0), std::make_shared<Constant>(3.0)});
      const auto ex0 = std::make_shared<Dummy>();

      // type ids get assigned in the order of use, which is program dependent, only check basic relations here
      REQUIRE(ex0->type_id() == Expr::get_type_id<Dummy>());
      REQUIRE(ex1->type_id() == Expr::get_type_id<Constant>());
      REQUIRE(ex4->type_id() == Expr::get_type_id<VecExpr<double>>());
      REQUIRE(ex4->type_id() < Expr::get_type_id<VecExpr<float>>());  // VecExpr<float> had not been used yet

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

    const auto ex4 = std::make_shared<VecExpr<double>>(std::initializer_list<double>{1.0, 2.0, 3.0});
    REQUIRE(begin(*ex4) != end(*ex4));
    REQUIRE(size(*ex4) == 3);
    REQUIRE(begin(ex4->expr()) == end(ex4->expr()));
    REQUIRE(size(ex4->expr()) == 0);

    const auto ex5_init =
        std::vector<std::shared_ptr<Constant>>{std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0),
                                               std::make_shared<Constant>(3.0)};
    const auto ex5 = std::make_shared<VecExpr<std::shared_ptr<Constant>>>(begin(ex5_init), end(ex5_init));
    REQUIRE(begin(*ex5) != end(*ex5));
    REQUIRE(size(*ex5) == 3);
    REQUIRE(begin(ex5->expr()) == end(ex5->expr()));
    REQUIRE(size(ex5->expr()) == 0);

    {
      auto ex6 = std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));
      REQUIRE(begin(*ex6) != end(*ex6));
      REQUIRE(size(*ex6) == 3);
      REQUIRE(begin(ex6->expr()) != end(ex6->expr()));
      REQUIRE(size(ex6->expr()) == 3);
      const auto& front_ptr = *begin(ex6->expr());
      auto front_ptr_cast = std::dynamic_pointer_cast<Constant>(front_ptr);
      REQUIRE(front_ptr_cast);
      REQUIRE(front_ptr_cast->value() == 1.0);
      {
        *begin(ex6->expr()) = std::make_shared<Constant>(4.0);
        auto &front_ptr = *begin(ex6->expr());
        auto front_ptr_cast = std::dynamic_pointer_cast<Constant>(front_ptr);
        REQUIRE(front_ptr_cast);
        REQUIRE(front_ptr_cast->value() == 4.0);
      }
    }
  }

  SECTION("scaled_product") {
    REQUIRE_NOTHROW(Product{});
    Product sp0{};
    REQUIRE(sp0.scalar() == 1.0);
    REQUIRE(sp0.factors().empty());

    REQUIRE_NOTHROW(sp0.append(2.0, std::make_shared<Dummy>()));
    REQUIRE(sp0.scalar() == 2.0);
    REQUIRE(sp0.factors().size() == 1);
    REQUIRE(*(sp0.factors()[0]) == Dummy{});

    REQUIRE(begin(sp0.expr()) != end(sp0.expr()));
    REQUIRE(ranges::size(sp0.expr()) == 1);
  }

  SECTION("latex") {
    Product sp0{};
    sp0.append(2.0, std::make_shared<Dummy>());
    REQUIRE(to_latex(sp0) == L"{{{2}} \\times {\\text{Dummy}}}");

    // VecExpr<shared_ptr<Expr>>
    {
      const auto ex5_init =
          std::vector<std::shared_ptr<Constant>>{std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0),
                                                 std::make_shared<Constant>(3.0)};
      auto ex6 =
          std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));
      REQUIRE(ex6->to_latex() ==
              L"{\\text{VecExpr}\\{{{{1}}} {{{2}}} {{{3}}} \\}}");
    }
  }

  SECTION("wolfram") {
    Product sp0{};
    sp0.append(2.0, std::make_shared<Dummy>());
    REQUIRE(to_wolfram(sp0) == L"Times[2,Dummy[]]");

    // VecExpr<shared_ptr<Expr>>
    {
      const auto ex5_init =
          std::vector<std::shared_ptr<Constant>>{std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0),
                                                 std::make_shared<Constant>(3.0)};
      auto ex6 = std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));
      REQUIRE(ex6->to_wolfram() == L"VecExpr[1,2,3]");
    }
  }

  SECTION("visitor") {
    {
      const auto ex5_init =
          std::vector<std::shared_ptr<Constant>>{std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0),
                                                 std::make_shared<Constant>(3.0)};
      ExprPtr ex6 = std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));

      auto ex = ex6 + ex6;

      latex_visitor v1{};
      ex->visit(v1);

//      std::wcout << "v1.result = " << v1.result << std::endl;
      REQUIRE(
          v1.result ==
          L"{{{1}}}{{{2}}}{{{3}}}{\\text{VecExpr}\\{{{{1}}} {{{2}}} {{{3}}} "
          L"\\}}{{{1}}}{{{2}}}{{{3}}}{\\text{VecExpr}\\{{{{1}}} {{{2}}} "
          L"{{{3}}} \\}}{ \\left({\\text{VecExpr}\\{{{{1}}} {{{2}}} {{{3}}} "
          L"\\}} + {\\text{VecExpr}\\{{{{1}}} {{{2}}} {{{3}}} \\}}\\right) }");

      latex_visitor v2{};
      ex->visit(v2, /* atoms_only = */ true);
      //      std::wcout << "v2.result = " << v2.result << std::endl;
      REQUIRE(v2.result ==
              L"{{{1}}}{{{2}}}{{{3}}}{{{1}}}{{{"
              L"2}}}{{{3}}}");
    }
  }

  SECTION("range") {

    {
      REQUIRE_NOTHROW(expr_range{});
      expr_range exrng{};
      REQUIRE(ranges::begin(exrng) == ranges::begin(exrng));
      REQUIRE(ranges::begin(exrng) == ranges::end(exrng));
    }

    // compares indices in address provided by cursor::address() to a list of indices
    auto compare = [](const container::svector<std::pair<ExprPtr*,int64_t>>& address1,
        std::initializer_list<int> address2) {
      return address1.size() == address2.size() &&
          std::equal(begin(address1), end(address1), begin(address2), [](const auto& parent_and_index1, const auto& index2) {
            return parent_and_index1.second == index2;
          });
    };

    {
      auto x = (ex<Constant>(1.0) + ex<Constant>(2.0)) * (ex<Constant>(3.0) + ex<Constant>(4.0));
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
      auto x = (ex<Constant>(1.0) +
                ex<Constant>(2.0) * (ex<Constant>(3.0) - ex<Dummy>())) *
               (ex<Constant>(5.0) +
                (ex<Constant>(6.0) + ex<Dummy>()) * ex<Constant>(8.0));
      REQUIRE_NOTHROW(expr_range{x});
      expr_range exrng{x};
      REQUIRE(std::distance(ranges::begin(exrng), ranges::end(exrng)) == 6);

      auto i = 0;
      for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
        switch (i) {
          case 0:
            REQUIRE(to_latex(*it) == L"{{{1}}}");
            REQUIRE( compare(ranges::get_cursor(it).address(), {0,0}) );
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
            REQUIRE( ranges::get_cursor(it).ordinal() == 5 );
            break;
        }
        ++i;
      }
    }
  }

  SECTION("expand") {
    {
      auto x =
          (ex<Constant>(1.0) + ex<Dummy>()) * (ex<Constant>(3.0) + ex<Dummy>());
      REQUIRE(to_latex(x) ==
              L"{{ \\left({{{1}}} + {\\text{Dummy}}\\right) }{ \\left({{{3}}} "
              L"+ {\\text{Dummy}}\\right) }}");
      expand(x);
      //      std::wcout << "x = " << to_latex(x) << std::endl;
      REQUIRE(to_latex(x) ==
              L"{ \\left({{{3}} \\times } + {{\\text{Dummy}}} + {{{3}} \\times "
              L"{\\text{Dummy}}} + {{\\text{Dummy}}{\\text{Dummy}}}\\right) }");
      simplify(x);
      //      std::wcout << "x = " << to_latex(x) << std::endl;
    }
    {
      auto x =
          (ex<Constant>(1.0) +
           ex<Constant>(2.0) * (ex<Constant>(3.0) - ex<Dummy>())) *
          (ex<Constant>(5.0) * (ex<Constant>(6.0) + ex<Dummy>()) + ex<Dummy>());
      //      std::wcout << "x = " << to_latex(x) << std::endl;
      REQUIRE(to_latex(x) ==
              L"{{ \\left({{{1}}} + {{{2}} \\times { \\left({{{3}}} + {{{-1}} "
              L"\\times {\\text{Dummy}}}\\right) }}\\right) }{ \\left({{{5}} "
              L"\\times { \\left({{{6}}} + {\\text{Dummy}}\\right) }} + "
              L"{\\text{Dummy}}\\right) }}");
      expand(x);
      //      std::wcout << "ex = " << to_latex(x) << std::endl;
      REQUIRE(to_latex(x) ==
              L"{ \\left({{{30}} \\times } + {{{5}} \\times {\\text{Dummy}}} + "
              L"{{\\text{Dummy}}} + {{{180}} \\times } + {{{30}} \\times "
              L"{\\text{Dummy}}} + {{{6}} \\times {\\text{Dummy}}} + {{{-60}} "
              L"\\times {\\text{Dummy}}} + {{{-10}} \\times "
              L"{\\text{Dummy}}{\\text{Dummy}}} + {{{-2}} \\times "
              L"{\\text{Dummy}}{\\text{Dummy}}}\\right) }");
    }
  }

  SECTION("hashing") {
    const auto ex5_init =
        std::vector<std::shared_ptr<Constant>>{std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0),
                                               std::make_shared<Constant>(3.0)};
    using boost::hash_value;
    REQUIRE_NOTHROW(hash_value(ex5_init));
    REQUIRE(hash_value(ex5_init) != hash_value(ex<Constant>(1)));
  }

  SECTION("commutativity") {
    const auto ex1 = std::make_shared<VecExpr<std::shared_ptr<Constant>>>(
        std::initializer_list<std::shared_ptr<Constant>>{
            std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0),
            std::make_shared<Constant>(3.0)});
    const auto ex2 =
        ex<Constant>(1.0) +
        (ex<Constant>(2.0) + ex<Constant>(3.0)) * ex<Constant>(4.0);

    REQUIRE(ex1->is_cnumber());
    REQUIRE(ex2->is_cnumber());
    REQUIRE(ex1->commutes_with(*ex1));
    REQUIRE(ex1->commutes_with(*ex2));
    REQUIRE(ex2->commutes_with(*ex1));
    REQUIRE(ex2->commutes_with(*ex2));
  }

}  // TEST_CASE("Expr"
