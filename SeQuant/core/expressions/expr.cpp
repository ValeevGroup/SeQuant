//
// Created by Eduard Valeyev on 2019-02-06.
//

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_iterator.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network/typedefs.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/transform.hpp>

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

void Sum::adjoint() {
  using namespace ranges;
  auto adj_summands = summands() | views::transform([](auto &&expr) {
                        return ::sequant::adjoint(expr);
                      });
  *this = Sum(ranges::begin(adj_summands), ranges::end(adj_summands));
}

ExprPtr Sum::canonicalize_impl(bool multipass, CanonicalizeOptions opts) {
  if (Logger::instance().canonicalize)
    std::wcout << "Sum::canonicalize_impl: input = "
               << to_latex_align(shared_from_this()) << std::endl;

  const auto npasses = multipass ? 2 : 1;
  for (auto pass = 0; pass != npasses; ++pass) {
    const auto rapid = (pass % 2 == 0);

    // canonicalizing TNs in a sum requires treating named indices as
    // meaningful/distinct
    auto opts_copy = opts;
    opts_copy.ignore_named_index_labels =
        CanonicalizeOptions::IgnoreNamedIndexLabel::No;
    if (rapid) {
      opts_copy.method = CanonicalizationMethod::Lexicographic;
    } else
      opts_copy.method = opts.method | CanonicalizationMethod::Topological;

    // recursively canonicalize summands ...
    // using for_each and direct access to summands
    sequant::for_each(summands_, [pass, &opts_copy, &rapid](ExprPtr &summand) {
      ExprPtr bp;
      if (rapid) {
        bp = summand->rapid_canonicalize(opts_copy);
      } else {
        bp = summand->canonicalize(opts_copy);
      }
      if (bp) {
        SEQUANT_ASSERT(bp->template is<Constant>());
        summand = ex<Product>(std::static_pointer_cast<Constant>(bp)->value(),
                              ExprPtrList{summand});
      }
    });
    if (Logger::instance().canonicalize)
      std::wcout << "Sum::canonicalize_impl (pass=" << pass
                 << "): after canonicalizing summands = "
                 << to_latex_align(shared_from_this()) << std::endl;

    HashingAccumulator acc;
    for (auto &summand : summands_) {
      acc.append(summand);
    }

    // last pass? sort by hash then by Expr::operator<
    // N.B. no point in differentiating between canonicalization methods here
    // since need to sort in both cases
    auto new_sum =
        (pass == npasses - 1) ? acc.make_canonicalized_sum() : acc.make_sum();
    this->swap(*new_sum);

    if (Logger::instance().canonicalize)
      std::wcout << "Sum::canonicalize_impl (pass=" << pass
                 << "): after reducing summands = "
                 << to_latex_align(shared_from_this()) << std::endl;
  }

  return {};  // side effects are absorbed into summands
}

HashingAccumulator &HashingAccumulator::append(ExprPtr summand, bool flatten) {
  // flatten, if needed
  if (flatten && summand.is<Sum>()) {
    for (auto &subsummand : summand.as<Sum>().summands()) {
      this->append(subsummand, flatten);
    }
    return *this;
  }

  // process summand as a whole
  auto it = summands_.find(summand);
  if (it == summands_.end()) {
    summands_.emplace(summand);
  } else {  // found existing term with the same hash
    auto existing_summand = *it;
    if (summand.template is<Product>()) {
      if (existing_summand.is<Product>()) {
        // both are products - add them
        existing_summand.as<Product>().add_identical(
            summand.template as<Product>());
      } else {
        // convert existing term to product and add
        auto product_copy = std::make_shared<Product>(summand->clone());
        product_copy->add_identical(existing_summand);
        summands_.erase(it);
        summands_.emplace(std::move(product_copy));
      }
    } else {
      if (existing_summand.is<Product>()) {
        existing_summand.as<Product>().add_identical(summand);
      } else {
        // neither is a product - create new product
        auto product_form = std::make_shared<Product>();
        product_form->append(2, summand.template as<Expr>());
        summands_.erase(it);
        summands_.emplace(std::move(product_form));
      }
    }
  }

  return *this;
}

SumPtr HashingAccumulator::make_sum_impl(bool canonicalize) {
  Sum::summands_type summands;
  summands.reserve(summands_.size());
  for (auto summand : summands_) {
    if (!summand->is_zero()) {
      summands.push_back(summand);
    }
  }

  if (canonicalize) {
    ranges::sort(summands, [](const auto &e1, const auto &e2) {
      if (e1->hash_value() == e2->hash_value()) {
        return e1 < e2;
      } else {
        return e1->hash_value() < e2->hash_value();
      }
    });
  }

  return std::make_shared<Sum>(std::move(summands), Sum::move_only_tag{});
}

SumPtr HashingAccumulator::make_sum() { return make_sum_impl(false); }

SumPtr HashingAccumulator::make_canonicalized_sum() {
  return make_sum_impl(true);
}

ExprPtr HashingAccumulator::make_expr(bool canonicalize) {
  if (summands_.size() == 0) {
    return ex<Constant>(0);
  } else if (summands_.size() == 1)
    return *(summands_.begin());
  else
    return make_sum_impl(canonicalize);
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
