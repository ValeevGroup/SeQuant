//
// Created by Eduard Valeyev on 2019-02-06.
//

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_iterator.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>

namespace sequant {

ExprIterator Expr::begin() { return begin_subexpr(); }

ExprIterator Expr::end() { return end_subexpr(); }

ConstExprIterator Expr::begin() const { return begin_subexpr(); }

ConstExprIterator Expr::end() const { return end_subexpr(); }

ConstExprIterator Expr::cbegin() const { return begin_subexpr(); }

ConstExprIterator Expr::cend() const { return end_subexpr(); }

ExprIterator Expr::begin_subexpr() { return ExprIterator{}; }

ExprIterator Expr::end_subexpr() { return ExprIterator{}; }

ConstExprIterator Expr::begin_subexpr() const { return ConstExprIterator{}; }

ConstExprIterator Expr::end_subexpr() const { return ConstExprIterator{}; }

std::size_t Expr::size() const { return end() - begin(); }

bool Expr::empty() const { return size() == 0; }

ExprPtr &Expr::operator[](std::size_t idx) {
  SEQUANT_ASSERT(idx < size());
  return begin()[idx];
}

const ExprPtr &Expr::operator[](std::size_t idx) const {
  SEQUANT_ASSERT(idx < size());
  return begin()[idx];
}

ExprPtr &Expr::at(std::size_t idx) { return (*this)[idx]; }

const ExprPtr &Expr::at(std::size_t idx) const { return (*this)[idx]; }

ExprPtr &Expr::front() { return at(0); }

const ExprPtr &Expr::front() const { return at(0); }

ExprPtr &Expr::back() { return at(size() - 1); }

const ExprPtr &Expr::back() const { return at(size() - 1); }

std::wstring Expr::to_latex() const {
  throw Exception("to_latex not implemented for " + type_name());
}

ExprPtr Expr::clone() const { return unique_copy(); }

bool proportional_to::operator()(const ExprPtr &expr1,
                                 const ExprPtr &expr2) const {
  if (expr1->type_id() !=
      expr2->type_id()) {  // if expr1 is a Product with single factor == expr2,
                           // or vice versa
    if (expr1.is<Product>()) {
      return expr1.as<Product>().factors().size() == 1 &&
             expr1.as<Product>().factors().front() == expr2;
    } else if (expr2.is<Product>()) {
      return expr2.as<Product>().factors().size() == 1 &&
             expr2.as<Product>().factors().front() == expr1;
    } else
      return false;
  }

  // expr1 and expr2 are same type

  if (expr1.is<Constant>()) {
    return true;
  }
  if (expr1.is<Product>()) {
    return expr1->hash_value() == expr2->hash_value() &&
           expr1.as<Product>().factors() == expr2.as<Product>().factors();
  }
  return expr1 == expr2;
}

}  // namespace sequant
