//
// Created by Eduard Valeyev on 2019-02-06.
//

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_iterator.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/tree_index.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <boost/core/demangle.hpp>

#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <typeindex>
#include <utility>

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

ExprPtr &ExprPtr::operator[](const TreeIndex &idx) {
  return idx.select_from(*this);
}

const ExprPtr &ExprPtr::operator[](const TreeIndex &idx) const {
  return idx.select_from(*this);
}

void Expr::throw_out_of_range(std::size_t idx) const {
  std::ostringstream oss;
  oss << "Expr::at(" << idx << "): index out of range (size=" << size()
      << ", type_name=" << type_name() << ")";
  throw Exception(oss.str());
}

ExprPtr &Expr::at(std::size_t idx) {
  if (idx >= size()) throw_out_of_range(idx);
  return begin()[static_cast<std::ptrdiff_t>(idx)];
}

const ExprPtr &Expr::at(std::size_t idx) const {
  if (idx >= size()) throw_out_of_range(idx);
  return begin()[static_cast<std::ptrdiff_t>(idx)];
}

ExprPtr &Expr::front() { return at(0); }

const ExprPtr &Expr::front() const { return at(0); }

ExprPtr &Expr::back() { return at(size() - 1); }

const ExprPtr &Expr::back() const { return at(size() - 1); }

std::wstring Expr::to_latex() const {
  throw Exception("to_latex not implemented for " + type_name());
}

void Expr::register_type_id(type_id_type id, const std::string &name,
                            std::type_index type) {
  static std::mutex mutex;
  static std::map<type_id_type, std::pair<std::string, std::type_index>>
      registry;
  std::scoped_lock lock(mutex);
  const auto [it, inserted] = registry.try_emplace(id, name, type);
  if (!inserted && it->second.second != type)
    throw Exception("Expr types " +
                    boost::core::demangle(it->second.second.name()) +
                    " (name \"" + it->second.first + "\") and " +
                    boost::core::demangle(type.name()) + " (name \"" + name +
                    "\") have the same type id " + std::to_string(id) +
                    "; declare a distinct static_type_name() for one of them");
}

Expr &Expr::operator[](const TreeIndex &idx) { return idx.select_from(*this); }

const Expr &Expr::operator[](const TreeIndex &idx) const {
  return idx.select_from(*this);
}

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
