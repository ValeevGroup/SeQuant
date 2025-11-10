//
// Created by Eduard Valeyev on 2019-02-06.
//

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network/typedefs.hpp>
#include <SeQuant/core/utility/macros.hpp>

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
  SEQUANT_ASSERT(this->operator bool());
  return *(this->get());
}

const Expr &ExprPtr::operator*() const & {
  SEQUANT_ASSERT(this->operator bool());
  return *(this->get());
}

Expr &&ExprPtr::operator*() && {
  SEQUANT_ASSERT(this->operator bool());
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

Expr &Expr::operator*=(const Expr &) { throw not_implemented("operator*="); }

Expr &Expr::operator^=(const Expr &) { throw not_implemented("operator^="); }

Expr &Expr::operator+=(const Expr &) { throw not_implemented("operator+="); }

Expr &Expr::operator-=(const Expr &) { throw not_implemented("operator-="); }

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

ExprPtr Product::canonicalize_impl(CanonicalizeOptions opts) {
  // recursively canonicalize non-tensor subfactors (tensors will be
  // canonicalized as part of the TN built of all tensor factors of this) ...
  ranges::for_each(factors_, [this, opts](auto &factor) {
    if (factor.template is<AbstractTensor>()) {
      return;
    }
    auto bp = factor->canonicalize(opts);
    if (bp) {
      SEQUANT_ASSERT(bp->template is<Constant>());
      this->scalar_ *= std::static_pointer_cast<Constant>(bp)->value();
    }
  });

  if (Logger::instance().canonicalize) {
    std::wcout << "Product canonicalization(" << to_wstring(opts.method)
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
    auto make_canonical_tn = [this, &opts](auto *tn_null_ptr) {
      using TN = std::decay_t<std::remove_pointer_t<decltype(tn_null_ptr)>>;
      ExprPtr canon_factor;
      TN tn(this->factors_);
      if constexpr (TN::version() == 3) {
        canon_factor = tn.canonicalize(
            TensorCanonicalizer::cardinal_tensor_labels(), opts);
      } else {
        using NamedIndexSet = tensor_network::NamedIndexSet;
        std::shared_ptr<NamedIndexSet> named_indices =
            !opts.named_indices
                ? nullptr
                : std::make_shared<NamedIndexSet>(opts.named_indices->begin(),
                                                  opts.named_indices->end());
        canon_factor = tn.canonicalize(
            TensorCanonicalizer::cardinal_tensor_labels(),
            opts.method == CanonicalizationMethod::Rapid, named_indices.get());
      }
      return std::pair{std::move(tn), canon_factor};
    };
    using TN = TensorNetwork;
    auto [tn, canon_factor] = make_canonical_tn(static_cast<TN *>(nullptr));

    const auto &tensors = tn.tensors();
    using std::size;
    SEQUANT_ASSERT(size(tensors) == size(factors_));
    using std::begin;
    using std::end;
    std::transform(begin(tensors), end(tensors), begin(factors_),
                   [](const auto &tptr) {
                     auto exprptr = std::dynamic_pointer_cast<Expr>(tptr);
                     SEQUANT_ASSERT(exprptr);
                     return exprptr;
                   });
    if (canon_factor) scalar_ *= canon_factor->template as<Constant>().value();
    this->reset_hash_value();
  } else {  // if contains non-tensors, do commutation-checking resort

    // comparer that respects cardinal tensor labels
    auto &cardinal_tensor_labels =
        TensorCanonicalizer::cardinal_tensor_labels();
    auto local_compare = [&cardinal_tensor_labels](const ExprPtr &first,
                                                   const ExprPtr &second) {
      if (first->is<Labeled>() && second->is<Labeled>()) {
        const auto first_label = first->as<Labeled>().label();
        const auto second_label = second->as<Labeled>().label();
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
          SEQUANT_ASSERT(!first_is_cardinal && !second_is_cardinal);
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
    std::wcout << "Product canonicalization(" << to_wstring(opts.method)
               << ") result: " << to_latex() << std::endl;

  return {};  // side effects are absorbed into the scalar_
}

void Product::adjoint() {
  SEQUANT_ASSERT(static_commutativity() == false);  // assert no slicing
  auto adj_scalar = conj(scalar());
  using namespace ranges;
  auto adj_factors =
      factors() | views::reverse |
      views::transform([](auto &expr) { return ::sequant::adjoint(expr); });
  using std::swap;
  *this =
      Product(adj_scalar, ranges::begin(adj_factors), ranges::end(adj_factors));
}

ExprPtr Product::canonicalize(CanonicalizeOptions opt) {
  return this->canonicalize_impl(opt);
}

ExprPtr Product::rapid_canonicalize(CanonicalizeOptions opt) {
  SEQUANT_ASSERT(opt.method == CanonicalizationMethod::Rapid);
  return this->canonicalize_impl(opt);
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

    // used to combine terms proportional to each other
    container::unordered_set<ExprPtr, sequant::hash::_<ExprPtr>,
                             proportional_to>
        summands;

    // process each summand
    for (auto &summand : summands_) {
      auto it = summands.find(summand);
      if (it == summands.end()) {
        summands.emplace(summand);
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
            summands.erase(it);
            summands.emplace(std::move(product_copy));
          }
        } else {
          if (existing_summand.is<Product>()) {
            existing_summand.as<Product>().add_identical(summand);
          } else {
            // neither is a product - create new product
            auto product_form = std::make_shared<Product>();
            product_form->append(2, summand.template as<Expr>());
            summands.erase(it);
            summands.emplace(std::move(product_form));
          }
        }
      }
    }
    decltype(summands_) new_summands;
    new_summands.reserve(summands.size());
    for (auto term : summands) {
      if (!term->template is<Product>() ||
          !term->template as<Product>().is_zero()) {
        new_summands.push_back(term);
      }
    }

    // last pass? sort by hash then by Expr::operator<
    // N.B. no point in differentiating between canonicalization methods here
    // since need to sort in both cases
    if (pass == npasses - 1) {
      ranges::sort(new_summands, [](const auto &e1, const auto &e2) {
        if (e1->hash_value() == e2->hash_value()) {
          return e1 < e2;
        } else {
          return e1->hash_value() < e2->hash_value();
        }
      });
    }
    summands_.swap(new_summands);

    if (Logger::instance().canonicalize)
      std::wcout << "Sum::canonicalize_impl (pass=" << pass
                 << "): after reducing summands = "
                 << to_latex_align(shared_from_this()) << std::endl;
  }

  return {};  // side effects are absorbed into summands
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

  if (expr1.is<Variable>()) {
    return true;
  }
  if (expr1.is<Product>()) {
    return expr1->hash_value() == expr2->hash_value() &&
           expr1.as<Product>().factors() == expr2.as<Product>().factors();
  }
  return expr1 == expr2;
}

}  // namespace sequant
