//
// Created by Eduard Valeyev on 2019-02-06.
//

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>

#include <range/v3/all.hpp>

#include <thread>
#include <vector>

namespace sequant {

ExprPtr ExprPtr::clone() const & {
  if (!*this) return {};
  return ExprPtr(as_shared_ptr()->clone());
}

ExprPtr ExprPtr::clone() && noexcept { return std::move(*this); }

ExprPtr::base_type &ExprPtr::as_shared_ptr() & {
  return static_cast<base_type &>(*this);
}
const ExprPtr::base_type &ExprPtr::as_shared_ptr() const & {
  return static_cast<const base_type &>(*this);
}
ExprPtr::base_type &&ExprPtr::as_shared_ptr() && {
  return static_cast<base_type &&>(*this);
}

Expr &ExprPtr::operator*() & {
  assert(this->operator bool());
  return *(this->get());
}

const Expr &ExprPtr::operator*() const & {
  assert(this->operator bool());
  return *(this->get());
}

Expr &&ExprPtr::operator*() && {
  assert(this->operator bool());
  return std::move(*(this->get()));
}

ExprPtr &ExprPtr::operator+=(const ExprPtr &other) {
  if (!*this) {
    *this = other.clone();
  } else if (as_shared_ptr()->is<Sum>()) {
    as_shared_ptr()->operator+=(*other);
  } else if (as_shared_ptr()->is<Constant>() && other->is<Constant>()) {
    *this = ex<Constant>(this->as<Constant>().value() +
                         other->as<Constant>().value());
  } else {
    *this = ex<Sum>(ExprPtrList{*this, other});
  }
  return *this;
}

ExprPtr &ExprPtr::operator-=(const ExprPtr &other) {
  if (!*this) {
    *this = ex<Constant>(-1) * other.clone();
  } else if (as_shared_ptr()->is<Sum>()) {
    as_shared_ptr()->operator-=(*other);
  } else if (as_shared_ptr()->is<Constant>() && other->is<Constant>()) {
    *this = ex<Constant>(this->as<Constant>().value() -
                         other->as<Constant>().value());
  } else {
    *this = ex<Sum>(ExprPtrList{*this, ex<Product>(-1, ExprPtrList{other})});
  }
  return *this;
}

ExprPtr &ExprPtr::operator*=(const ExprPtr &other) {
  if (!*this) {
    *this = other.clone();
  } else if (as_shared_ptr()->is<Product>()) {
    as_shared_ptr()->operator*=(*other);
  } else if (as_shared_ptr()->is<Constant>() && other->is<Constant>()) {
    *this = ex<Constant>(this->as<Constant>().value() *
                         other->as<Constant>().value());
  } else {
    *this = ex<Product>(ExprPtrList{*this, other});
  }
  return *this;
}

std::size_t ExprPtr::size() const { return this->get()->size(); }

std::wstring ExprPtr::to_latex() const { return as_shared_ptr()->to_latex(); }

std::logic_error Expr::not_implemented(const char *fn) const {
  std::ostringstream oss;
  oss << "Expr::" << fn
      << " not implemented in this derived class (type_name=" << type_name()
      << ")";
  return std::logic_error(oss.str().c_str());
}

std::wstring Expr::to_latex() const { throw not_implemented("to_latex"); }

std::wstring Expr::to_wolfram() const { throw not_implemented("to_wolfram"); }

ExprPtr Expr::clone() const { throw not_implemented("clone"); }

void Expr::adjoint() { throw not_implemented("adjoint"); }

Expr &Expr::operator*=(const Expr &that) {
  throw not_implemented("operator*=");
}

Expr &Expr::operator^=(const Expr &that) {
  throw not_implemented("operator^=");
}

Expr &Expr::operator+=(const Expr &that) {
  throw not_implemented("operator+=");
}

Expr &Expr::operator-=(const Expr &that) {
  throw not_implemented("operator-=");
}

ExprPtr adjoint(const ExprPtr &expr) {
  auto result = expr->clone();
  result->adjoint();
  return result;
}

void Constant::adjoint() {
  value_ = conj(value_);
  reset_hash_value();
}

std::wstring_view Variable::label() const { return label_; }

void Variable::set_label(std::wstring label) {
  label_ = std::move(label);
  reset_hash_value();
}

void Variable::conjugate() { conjugated_ = !conjugated_; }

bool Variable::conjugated() const { return conjugated_; }

std::wstring Variable::to_latex() const {
  std::wstring result = L"{" + utf_to_latex(label_) + L"}";
  if (conjugated_) result = L"{" + result + L"^*" + L"}";
  return result;
}

ExprPtr Variable::clone() const { return ex<Variable>(*this); }

void Variable::adjoint() { conjugate(); }

bool Product::is_commutative() const {
  bool result = true;
  const auto nfactors = size();
  for (size_t f = 0; f != nfactors; ++f) {
    for (size_t s = 1; result && s != nfactors; ++s) {
      result &= factors_[f]->commutes_with(*factors_[s]);
    }
  }
  return result;
}

ExprPtr Product::canonicalize_impl(bool rapid) {
  // recursively canonicalize subfactors ...
  ranges::for_each(factors_, [this](auto &factor) {
    auto bp = factor->canonicalize();
    if (bp) {
      assert(bp->template is<Constant>());
      this->scalar_ *= std::static_pointer_cast<Constant>(bp)->value();
    }
  });

  if (Logger::instance().canonicalize) {
    std::wcout << "Product canonicalization(" << (rapid ? "fast" : "slow")
               << ") input: " << to_latex() << std::endl;
  }

  // pull out all Variables to the front
  auto variables = factors_ | ranges::views::filter([](const auto &factor) {
                     return factor.template is<Variable>();
                   }) |
                   ranges::to_vector;
  // sort variables
  ranges::sort(variables, [](const auto &first, const auto &second) {
    return first.template as<Variable>().label() <
           second.template as<Variable>().label();
  });

  factors_ = factors_ | ranges::views::filter([](const auto &factor) {
               return !factor.template is<Variable>();
             }) |
             ranges::to<decltype(factors_)>;

  // if there are no factors, insert variables back and return
  if (factors_.empty()) {
    factors_.insert(factors_.begin(), variables.begin(), variables.end());
    return {};
  }

  auto contains_nontensors = ranges::any_of(factors_, [](const auto &factor) {
    return std::dynamic_pointer_cast<AbstractTensor>(factor) == nullptr;
  });
  if (!contains_nontensors) {  // tensor network canonization is a special case
                               // that's done in
                               // TensorNetwork
    TensorNetworkV2 tn(factors_);
    auto canon_factor =
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), rapid);
    const auto &tensors = tn.tensors();
    using std::size;
    assert(size(tensors) == size(factors_));
    using std::begin;
    using std::end;
    std::transform(begin(tensors), end(tensors), begin(factors_),
                   [](const auto &tptr) {
                     auto exprptr = std::dynamic_pointer_cast<Expr>(tptr);
                     assert(exprptr);
                     return exprptr;
                   });
    if (canon_factor) scalar_ *= canon_factor->as<Constant>().value();
    this->reset_hash_value();
  } else {  // if contains non-tensors, do commutation-checking resort

    // comparer that respects cardinal tensor labels
    auto &cardinal_tensor_labels =
        TensorCanonicalizer::cardinal_tensor_labels();
    auto local_compare = [&cardinal_tensor_labels](const ExprPtr &first,
                                                   const ExprPtr &second) {
      if (first->is<Labeled>() && second->is<Labeled>()) {
        const auto first_label = first->as<Tensor>().label();
        const auto second_label = second->as<Tensor>().label();
        if (first_label == second_label) return *first < *second;
        const auto first_is_cardinal_it = ranges::find_if(
            cardinal_tensor_labels,
            [&first_label](const std::wstring &l) { return l == first_label; });
        const auto first_is_cardinal =
            first_is_cardinal_it != ranges::end(cardinal_tensor_labels);
        const auto second_is_cardinal_it = ranges::find_if(
            cardinal_tensor_labels, [&second_label](const std::wstring &l) {
              return l == second_label;
            });
        const auto second_is_cardinal =
            second_is_cardinal_it != ranges::end(cardinal_tensor_labels);
        if (first_is_cardinal && second_is_cardinal)
          return first_is_cardinal_it < second_is_cardinal_it;
        else if (first_is_cardinal && !second_is_cardinal)
          return true;
        else if (!first_is_cardinal && second_is_cardinal)
          return false;
        else {
          assert(!first_is_cardinal && !second_is_cardinal);
          return *first < *second;
        }
      } else
        return *first < *second;
    };

    // ... then resort, respecting commutativity
    using std::begin;
    using std::end;
    if (static_commutativity()) {
      if (is_commutative()) {
        std::stable_sort(begin(factors_), end(factors_), local_compare);
      }
    } else {
      // must do bubble sort if not commuting to avoid swapping elements across
      // a noncommuting element
      bubble_sort(
          begin(factors_), end(factors_),
          [&local_compare](const ExprPtr &first, const ExprPtr &second) {
            bool result = (first->commutes_with(*second))
                              ? local_compare(first, second)
                              : false;
            return result;
          });
    }
  }
  // reinsert Variables at the front
  factors_.insert(factors_.begin(), variables.begin(), variables.end());

  // TODO evaluate product of Tensors (turn this into Products of Products)

  if (Logger::instance().canonicalize)
    std::wcout << "Product canonicalization(" << (rapid ? "fast" : "slow")
               << ") result: " << to_latex() << std::endl;

  return {};  // side effects are absorbed into the scalar_
}

void Product::adjoint() {
  assert(static_commutativity() == false);  // assert no slicing
  auto adj_scalar = conj(scalar());
  using namespace ranges;
  auto adj_factors =
      factors() | views::reverse |
      views::transform([](auto &expr) { return ::sequant::adjoint(expr); });
  using std::swap;
  *this =
      Product(adj_scalar, ranges::begin(adj_factors), ranges::end(adj_factors));
}

ExprPtr Product::canonicalize() {
  return this->canonicalize_impl(/* rapid = */ false);
}

ExprPtr Product::rapid_canonicalize() {
  return this->canonicalize_impl(/* rapid = */ true);
}

void CProduct::adjoint() {
  auto adj_scalar = conj(scalar());
  using namespace ranges;
  // no need to reverse for commutative product
  auto adj_factors = factors() | views::transform([](auto &&expr) {
                       return ::sequant::adjoint(expr);
                     });
  *this = CProduct(adj_scalar, ranges::begin(adj_factors),
                   ranges::end(adj_factors));
}

void NCProduct::adjoint() {
  auto adj_scalar = conj(scalar());
  using namespace ranges;
  // no need to reverse for commutative product
  auto adj_factors =
      factors() | views::reverse |
      views::transform([](auto &&expr) { return ::sequant::adjoint(expr); });
  *this = NCProduct(adj_scalar, ranges::begin(adj_factors),
                    ranges::end(adj_factors));
}

void Sum::adjoint() {
  using namespace ranges;
  auto adj_summands = summands() | views::transform([](auto &&expr) {
                        return ::sequant::adjoint(expr);
                      });
  *this = Sum(ranges::begin(adj_summands), ranges::end(adj_summands));
}

ExprPtr Sum::canonicalize_impl(bool multipass) {
  if (Logger::instance().canonicalize)
    std::wcout << "Sum::canonicalize_impl: input = "
               << to_latex_align(shared_from_this()) << std::endl;

  const auto npasses = multipass ? 3 : 1;
  for (auto pass = 0; pass != npasses; ++pass) {
    // recursively canonicalize summands ...
    const auto nsubexpr = ranges::size(*this);
    for (std::size_t i = 0; i != nsubexpr; ++i) {
      auto bp = (pass % 2 == 0) ? summands_[i]->rapid_canonicalize()
                                : summands_[i]->canonicalize();
      if (bp) {
        assert(bp->template is<Constant>());
        summands_[i] =
            ex<Product>(std::static_pointer_cast<Constant>(bp)->value(),
                        ExprPtrList{summands_[i]});
      }
    };

    if (Logger::instance().canonicalize)
      std::wcout << "Sum::canonicalize_impl (pass=" << pass
                 << "): after canonicalizing summands = "
                 << to_latex_align(shared_from_this()) << std::endl;

    // ... then resort according to size, then hash values
    using std::begin;
    using std::end;
    std::stable_sort(begin(summands_), end(summands_),
                     [](const auto &first, const auto &second) {
                       const auto first_size = sequant::size(first);
                       const auto second_size = sequant::size(second);

                       return (first_size == second_size)
                                  ? *first < *second
                                  : first_size < second_size;
                     });

    if (Logger::instance().canonicalize)
      std::wcout << "Sum::canonicalize_impl (pass=" << pass
                 << "): after hash-sorting summands = "
                 << to_latex_align(shared_from_this()) << std::endl;

    // ... then reduce terms whose hash values are identical
    auto first_it = begin(summands_);
    auto hash_comparer = [](const auto &first, const auto &second) {
      return first->hash_value() == second->hash_value();
    };
    while ((first_it = std::adjacent_find(first_it, end(summands_),
                                          hash_comparer)) != end(summands_)) {
      assert((*first_it)->hash_value() == (*(first_it + 1))->hash_value());
      // find first element whose hash is not equal to (*first_it)->hash_value()
      auto plast_it = std::find_if_not(
          first_it + 1, end(summands_), [first_it](const auto &elem) {
            return (*first_it)->hash_value() == elem->hash_value();
          });
      const auto nidentical = plast_it - first_it;
      assert(nidentical > 1);
      // combine all identical summands into Product
      auto reduce_range = [first_it, this, nidentical](auto &begin, auto &end) {
        if ((*first_it)->is<Product>()) {  // handle group of Products
          auto &prod = (*first_it)->as<Product>();
          for (auto it = begin + 1; it != end; ++it) {
            // convert to Product if not already
            if (!(*it)->template is<Product>()) {
              *it = std::make_shared<Product>(1, ExprPtrList{*it});
            }
            prod.add_identical((*it)->template as<Product>());
          }
          this->summands_.erase(prod.is_zero() ? first_it : first_it + 1, end);
        } else {  // handle all other types
          auto product_form = std::make_shared<Product>();
          product_form->append(nidentical, (*first_it)->as<Expr>());
          *first_it = product_form;
          this->summands_.erase(first_it + 1, end);
        }
      };
      reduce_range(first_it, plast_it);
    }

    if (Logger::instance().canonicalize)
      std::wcout << "Sum::canonicalize_impl (pass=" << pass
                 << "): after reducing summands = "
                 << to_latex_align(shared_from_this()) << std::endl;
  }

  return {};  // side effects are absorbed into summands
}

}  // namespace sequant
