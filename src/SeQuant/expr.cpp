//
// Created by Eduard Valeyev on 2019-02-06.
//

#include "./expr.hpp"
#include "./tensor_network.hpp"
#include "./utility.hpp"

namespace sequant {

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

std::shared_ptr<Expr> Product::canonicalize_impl(bool rapid) {
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

  try {  // tensor tensor canonization is a special case that's done in TensorNetwork
    TensorNetwork tn(factors_);
    auto canon_factor =
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), rapid);
    const auto &tensors = tn.tensors();
    using std::size;
    assert(size(tensors) == size(factors_));
    using std::begin;
    using std::end;
    std::copy(begin(tensors), end(tensors), begin(factors_));
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

  // TODO factorize product of Tensors (turn this into Products of Products

  if (Logger::get_instance().canonicalize)
    std::wcout << "Product canonicalization(" << (rapid ? "fast" : "slow")
               << ") result: " << to_latex() << std::endl;

  return {};  // side effects are absorbed into the scalar_
}

std::shared_ptr<Expr> Product::canonicalize() {
  return this->canonicalize_impl(/* rapid = */ false);
}

std::shared_ptr<Expr> Product::rapid_canonicalize() {
  return this->canonicalize_impl(/* rapid = */ true);
}

// std::shared_ptr<Expr> Product::rapid_canonicalize() {
//  // recursively canonicalize subfactors ...
//  ranges::for_each(factors_, [this](auto &factor) {
//    auto bp = factor->canonicalize();
//    if (bp) {
//      assert(bp->template is<Constant>());
//      this->scalar_ *= std::static_pointer_cast<Constant>(bp)->value();
//    }
//  });
//
//  // ... then resort
//  using std::begin;
//  using std::end;
//  // default sorts by type, then by hash
//  // TODO for same types see if that type's operator< is defined, otherwise
//  use hashes std::stable_sort(begin(factors_), end(factors_), [](const auto
//  &first, const auto &second) {
//    const auto first_id = first->type_id();
//    const auto second_id = second->type_id();
//    if (first_id == second_id) {
//      return first->hash_value() < second->hash_value();
//    } else // first_id != second_id
//      return first_id < second_id;
//  });
//
//  // TODO factorize product of Tensors (turn this into Products of Products
//
//  return {};  // side effects are absorbed into the scalar_
//}

ExprPtr Sum::canonicalize_impl(bool multipass) {

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

    // ... then resort according to hash values
    using std::begin;
    using std::end;
    std::stable_sort(begin(summands_), end(summands_), [](const auto &first, const auto &second) {
      return *first < *second;
    });

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
        assert((*first_it)->template is<Product>());
        for (auto it = begin; it != end; ++it) {
          if (it != first_it) {
            assert((*it)->template is<Product>());
            std::static_pointer_cast<Product>(*first_it)->add_identical(std::static_pointer_cast<Product>(*it));
          }
        }
        this->summands_.erase(first_it + 1, end);
      };
      reduce_range(first_it, plast_it);
    }
  }

  return {};  // side effects are absorbed into summands
}


}  // namespace sequant
