//
// Created by Eduard Valeyev on 2019-02-06.
//

#include "./expr.hpp"
#include "./tensor_network.hpp"

namespace sequant2 {

bool debug_canonicalize = false;

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

  // ... then resort, respecting commutativity
  using std::begin;
  using std::end;
  if (static_commutativity()) {
    if (is_commutative()) {
      std::stable_sort(begin(factors_), end(factors_));
    }
  }
  else {
    // must do bubble sort if not commuting to avoid swapping elements across a noncommuting element
    bubble_sort(begin(factors_), end(factors_), [](const ExprPtr& first, const ExprPtr& second) {
      bool result = (first->commutes_with(*second)) ? (*first < *second) : true;
      return result;
    });
  }

  try {
    if (debug_canonicalize)
      std::wcout << "Product canonicalization input: " << to_latex() << std::endl;
    TensorNetwork tn(factors_);
    auto canon_factor =
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), rapid);
    const auto &tensors = tn.tensors();
    using std::size;
    assert(size(tensors) == size(factors_));
    using std::begin;
    using std::end;
    std::copy(begin(tensors), end(tensors), begin(factors_));
    if (canon_factor)
      scalar_ *= canon_factor->as<Constant>().value();
    this->reset_hash_value();
    if (debug_canonicalize)
      std::wcout << "Product canonicalization result: " << to_latex() << std::endl;
  } catch (std::logic_error &) {
  }  // do nothing if contains non-tensors

  // TODO factorize product of Tensors (turn this into Products of Products

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

}  // namespace sequant2
