//
// Created by Eduard Valeyev on 2019-02-06.
//

#include "expr.hpp"
#include "tensor_network.hpp"
#include "tensor.hpp"
#include "utility.hpp"

namespace sequant {

std::logic_error
Expr::not_implemented(const char* fn) const {
  std::ostringstream oss;
  oss << "Expr::" << fn << " not implemented in this derived class (type_name=" << type_name() << ")";
  return std::logic_error(oss.str().c_str());
}

std::wstring Expr::to_latex() const {
  throw not_implemented("to_latex");
}

std::wstring Expr::to_wolfram() const{
  throw not_implemented("to_wolfram");
}

ExprPtr Expr::clone() const {
  throw not_implemented("clone");
}

void Expr::adjoint() {
  throw not_implemented("adjoint");
}

Expr& Expr::operator*=(const Expr &that) {
  throw not_implemented("operator*=");
}

Expr& Expr::operator^=(const Expr &that) {
  throw not_implemented("operator^=");
}

Expr& Expr::operator+=(const Expr &that) {
  throw not_implemented("operator+=");
}

Expr& Expr::operator-=(const Expr &that) {
  throw not_implemented("operator-=");
}

ExprPtr adjoint(const ExprPtr& expr) {
  auto result = expr->clone();
  result->adjoint();
  return result;
}

void Constant::adjoint() {
  value_ = conj(value_);
  reset_hash_value();
}

bool Product::is_commutative() const {
  bool result = true;
  const auto nfactors = size();
  for(size_t f=0; f!= nfactors; ++f) {
    for(size_t s=1; result && s != nfactors; ++s) {
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

  if (Logger::get_instance().canonicalize) {
    std::wcout << "Product canonicalization(" << (rapid ? "fast" : "slow")
               << ") input: " << to_latex() << std::endl;
  }

  try {  // tensor network canonization is a special case that's done in TensorNetwork
    TensorNetwork tn(factors_);
    auto canon_factor =
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), rapid);
    const auto &tensors = tn.tensors();
    using std::size;
    assert(size(tensors) == size(factors_));
    using std::begin;
    using std::end;
    std::transform(begin(tensors), end(tensors), begin(factors_), [](const auto& tptr) {
      auto exprptr = std::dynamic_pointer_cast<Expr>(tptr);
      assert(exprptr);
      return exprptr;
    });
    if (canon_factor) scalar_ *= canon_factor->as<Constant>().value();
    this->reset_hash_value();
  } catch (std::logic_error
               &) {  // if contains non-tensors, do commutation-checking resort

    // comparer that respects cardinal tensor labels
    auto &cardinal_tensor_labels =
        TensorCanonicalizer::cardinal_tensor_labels();
    auto local_compare = [&cardinal_tensor_labels](const ExprPtr &first,
                                                   const ExprPtr &second) {
      if (first->is<Tensor>() && second->is<Tensor>()) {
        const auto first_label = first->as<Tensor>().label();
        const auto second_label = second->as<Tensor>().label();
        const bool first_is_cardinal = ranges::any_of(
            cardinal_tensor_labels,
            [&first_label](const std::wstring &l) { return l == first_label; });
        const bool second_is_cardinal = ranges::any_of(
            cardinal_tensor_labels, [&second_label](const std::wstring &l) {
              return l == second_label;
            });
        if (!(first_is_cardinal ^ second_is_cardinal))
          return *first < *second;
        else if (first_is_cardinal)
          return true;
        else  // if (second_is_cardinal)
          return false;
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
                              : true;
            return result;
          });
    }
  }

  // TODO evaluate product of Tensors (turn this into Products of Products)

  if (Logger::get_instance().canonicalize)
    std::wcout << "Product canonicalization(" << (rapid ? "fast" : "slow")
               << ") result: " << to_latex() << std::endl;

  return {};  // side effects are absorbed into the scalar_
}

void Product::adjoint() {
  assert(static_commutativity() == false);  // assert no slicing
  auto adj_scalar = conj(scalar());
  using namespace ranges;
  auto adj_factors = factors() | views::reverse | views::transform([](auto& expr) { return ::sequant::adjoint(expr); });
  using std::swap;
  *this = Product(adj_scalar, ranges::begin(adj_factors), ranges::end(adj_factors));
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
  auto adj_factors = factors() | views::transform([](auto&& expr) { return ::sequant::adjoint(expr); });
  *this = CProduct(adj_scalar, ranges::begin(adj_factors), ranges::end(adj_factors));
}

void NCProduct::adjoint() {
  auto adj_scalar = conj(scalar());
  using namespace ranges;
  // no need to reverse for commutative product
  auto adj_factors = factors() | views::reverse | views::transform([](auto&& expr) { return ::sequant::adjoint(expr); });
  *this = NCProduct(adj_scalar, ranges::begin(adj_factors), ranges::end(adj_factors));
}

void Sum::adjoint() {
  using namespace ranges;
  auto adj_summands = summands() | views::transform([](auto&& expr) { return ::sequant::adjoint(expr); });
  *this = Sum(ranges::begin(adj_summands), ranges::end(adj_summands));
}

ExprPtr Sum::canonicalize_impl(bool multipass) {

  if (Logger::get_instance().canonicalize) std::wcout << "Sum::canonicalize_impl: input = " << to_latex_align(shared_from_this()) << std::endl;

  const auto npasses = multipass ? 3 : 1;
  for (auto pass = 0; pass != npasses; ++pass) {
    // recursively canonicalize summands ...
    const auto nsubexpr = ranges::size(*this);
    for (std::size_t i = 0; i != nsubexpr; ++i) {
      auto bp = (pass % 2 == 0) ? summands_[i]->rapid_canonicalize() : summands_[i]->canonicalize();
      if (bp) {
        assert(bp->template is<Constant>());
        summands_[i] = ex<Product>(std::static_pointer_cast<Constant>(bp)->value(), ExprPtrList{summands_[i]});
      }
    };

    if (Logger::get_instance().canonicalize) std::wcout << "Sum::canonicalize_impl (pass=" << pass << "): after canonicalizing summands = " << to_latex_align(shared_from_this()) << std::endl;

    // ... then resort according to size, then hash values
    using std::begin;
    using std::end;
    std::stable_sort(begin(summands_), end(summands_), [](const auto &first, const auto &second) {
      const auto first_size = ranges::size(*first);
      const auto second_size = ranges::size(*second);

      return (first_size == second_size) ? *first < *second : first_size < second_size;
    });

    if (Logger::get_instance().canonicalize) std::wcout << "Sum::canonicalize_impl (pass=" << pass << "): after hash-sorting summands = " << to_latex_align(shared_from_this()) << std::endl;

    // ... then reduce terms whose hash values are identical
    auto first_it = begin(summands_);
    auto hash_comparer = [](const auto &first, const auto &second) {
      return first->hash_value() == second->hash_value();
    };
    while ((first_it = std::adjacent_find(first_it, end(summands_), hash_comparer)) != end(summands_)) {
      assert((*first_it)->hash_value() == (*(first_it + 1))->hash_value());
      // find first element whose hash is not equal to (*first_it)->hash_value()
      auto plast_it = std::find_if_not(first_it + 1, end(summands_), [first_it](const auto &elem) {
        return (*first_it)->hash_value() == elem->hash_value();
      });
      assert(plast_it - first_it > 1);
      auto reduce_range = [first_it, this](auto &begin, auto &end) {

        // Converting tensor into a product; TODO: Use separate add_identical function for Tensor
        if ((*first_it)->template is<Tensor>()){
          assert((*first_it)->template is<Tensor>());
          Product tensor_as_Product{};
          tensor_as_Product.append(1.0,(*first_it)->as<Tensor>());
          (*first_it) = std::make_shared<Product>(tensor_as_Product);
        }

        assert((*first_it)->template is<Product>());
        for (auto it = begin; it != end; ++it) {
          if (it != first_it) {
            if ((*it)->template is<Tensor>()){
              assert((*it)->template is<Tensor>());
              Product tensor_as_Product{};
              tensor_as_Product.append(1.0,(*it)->template as<Tensor>());
              (*it) = std::make_shared<Product>(tensor_as_Product);
            }
            assert((*it)->template is<Product>());
            std::static_pointer_cast<Product>(*first_it)->add_identical(std::static_pointer_cast<Product>(*it));
          }
        }
        this->summands_.erase(first_it + 1, end);
      };
      reduce_range(first_it, plast_it);
    }

    if (Logger::get_instance().canonicalize) std::wcout << "Sum::canonicalize_impl (pass=" << pass << "): after reducing summands = " << to_latex_align(shared_from_this()) << std::endl;

  }

  return {};  // side effects are absorbed into summands
}


}  // namespace sequant
