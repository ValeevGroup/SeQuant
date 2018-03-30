//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <iostream>
#include "../../src/SeQuant2/wick.hpp"

struct Dummy : public sequant2::Expr {
  virtual ~Dummy() = default;
  std::wstring to_latex() const override {
    return L"{\\text{Dummy}}";
  }
  type_id_type type_id() const override {
    return get_type_id<Dummy>();
  };
  bool static_compare(const sequant2::Expr& that) const override {
    return true;
  }
};

template<typename T>
struct VecExpr : public std::vector<T>, public sequant2::Expr {
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
      if constexpr (sequant2::Expr::is_shared_ptr_of_expr_or_derived<T>::value) {
        result += e->to_latex() + L" ";
      } else {
        result += std::to_wstring(e) + L" ";
      }
    }
    result += L"\\}}";
    return result;
  }

  type_id_type type_id() const override{
    return get_type_id<VecExpr<T>>();
  };

 private:
  cursor begin_cursor() const override {
    if constexpr (sequant2::Expr::is_shared_ptr_of_expr<T>::value) {
      return base_type::empty() ? Expr::begin_cursor() : cursor{&base_type::at(0)};
    } else {
      return Expr::begin_cursor();
    }
  };
  cursor end_cursor() const override {
    if constexpr (sequant2::Expr::is_shared_ptr_of_expr<T>::value) {
      return base_type::empty() ? Expr::end_cursor() : cursor{&base_type::at(0) + base_type::size()};
    } else {
      return Expr::end_cursor();
    }
  };

  bool static_compare(const sequant2::Expr& that) const override {
    return static_cast<const base_type&>(*this) == static_cast<const base_type&>(static_cast<const VecExpr&>(that));
  }

};

struct latex_visitor {
  void operator()(const std::shared_ptr<sequant2::Expr>& expr) {
    result += expr->to_latex();
  }
  std::wstring result {};
};

TEST_CASE("Expr", "[elements]") {

  using namespace sequant2;

  SECTION("expr") {
    REQUIRE_NOTHROW(std::make_shared<Dummy>());
    const auto ex1 = std::make_shared<Dummy>();
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
  }

  SECTION("comparison") {

    {
      const auto ex0 = std::make_shared<Dummy>();
      const auto ex1 = std::make_shared<Constant>(1);
      const auto ex2 = std::make_shared<Constant>(2);
      const auto ex3 = std::make_shared<Constant>(1);
      const auto ex4 = std::make_shared<VecExpr<double>>();
      const auto ex5 =
          std::make_shared<VecExpr<ExprPtr>>(ExprPtrList{
              std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0), std::make_shared<Constant>(3.0)});

      REQUIRE(ex0->type_id() == 1);  // Dummy does not have a type_id
      REQUIRE(ex1->type_id() == 2);  // Constant does
      REQUIRE(ex4->type_id() == 3);  // so does VecExpr<double>
      REQUIRE(ex5->type_id() == 4);  // so does VecExpr<ExprPtr>

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
    REQUIRE(to_latex(sp0) == L"{{2.000000} \\times {\\text{Dummy}}}");

    // VecExpr<shared_ptr<Expr>>
    {
      const auto ex5_init =
          std::vector<std::shared_ptr<Constant>>{std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0),
                                                 std::make_shared<Constant>(3.0)};
      auto ex6 = std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));
      REQUIRE(ex6->to_latex() == L"{\\text{VecExpr}\\{{1.000000} {2.000000} {3.000000} \\}}");
    }
  }

  SECTION("visitor") {
    {
      latex_visitor v{};

      const auto ex5_init =
          std::vector<std::shared_ptr<Constant>>{std::make_shared<Constant>(1.0), std::make_shared<Constant>(2.0),
                                                 std::make_shared<Constant>(3.0)};
      ExprPtr ex6 = std::make_shared<VecExpr<ExprPtr>>(begin(ex5_init), end(ex5_init));

      auto ex = ex6 + ex6;
      ex->visit(v);

      REQUIRE(v.result == L"{1.000000}{2.000000}{3.000000}{\\text{VecExpr}\\{{1.000000} {2.000000} {3.000000} \\}}{\\text{VecExpr}\\{{1.000000} {2.000000} {3.000000} \\}}{1.000000}{2.000000}{3.000000}{\\text{VecExpr}\\{{1.000000} {2.000000} {3.000000} \\}}{\\text{VecExpr}\\{{1.000000} {2.000000} {3.000000} \\}}{ \\left({\\text{VecExpr}\\{{1.000000} {2.000000} {3.000000} \\}} + {\\text{VecExpr}\\{{1.000000} {2.000000} {3.000000} \\}}\\right) }");
    }
  }

}  // TEST_CASE("Expr"