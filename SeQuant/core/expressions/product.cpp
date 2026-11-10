#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network/typedefs.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/equal.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>

namespace sequant {

Product::Product(ExprPtrList factors, Flatten flatten_tag) {
  using std::begin;
  using std::end;
  for (auto it = begin(factors); it != end(factors); ++it)
    append(1, *it, flatten_tag);
}

Product &Product::append(ExprPtr factor, Flatten flatten_tag) {
  return this->append(1, factor, flatten_tag);
}

const Product::scalar_type &Product::scalar() const { return scalar_; }

bool Product::is_zero() const { return Constant::is_zero(this->scalar()); }

const Product::factors_type &Product::factors() const { return factors_; }
Product::factors_type &Product::factors() {
  // N.B. handing out a mutable reference to factors_ invalidates the memoized
  // hash.  Unlike begin_subexpr()/end_subexpr(), this must not be guarded by
  // `!factors_.empty()`: those hand out iterators, and an empty range yields
  // nothing to dereference, whereas the caller can grow an empty factors_
  // through this reference (e.g. push_back()).
  // N.B. the canonical mark need not be reset: changes of the factors are
  // detected by Expr::is_canonical()
  hash_value_.reset();
  return factors_;
}

const ExprPtr &Product::factor(size_t i) const { return factors_.at(i); }

bool Product::empty() const { return factors_.empty(); }

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

ExprPtr Product::canonicalize_impl(const NamedIndexSpec &named_indices,
                                   CanonicalizeOptions opts) {
  reset_canonical_mark();

  // a tree of products of tensors, e.g. a product parenthesized for
  // evaluation, is one tensor network: its tensors are canonicalized together,
  // so that the dummy indices of the nested products are named consistently,
  // and its parenthesization is preserved
  if (container::svector<ExprPtr *> leaves;
      collect_nested_tensor_leaves(leaves)) {
    canonicalize_nested(leaves, named_indices, opts);
    return {};
  }

  // recursively canonicalize non-tensor subfactors (tensors will be
  // canonicalized as part of the TN built of all tensor factors of this) ...
  // an index that a subfactor shares with another factor is a dummy of this
  // product, but external to the subfactor, so it is named in the subfactor;
  // the other indices of the other factors are named, too: they do not occur
  // in the subfactor, but naming them keeps its dummies from being renamed to
  // them (which would contract them with the other factors)
  // N.B. the indices of each factor are collected once
  container::svector<container::set<Index>> factor_indices;
  auto subfactor_named_indices = [this, &named_indices,
                                  &factor_indices](std::size_t f) {
    if (factor_indices.empty())
      for (const auto &factor : factors_)
        factor_indices.push_back(get_used_indices(factor));
    container::set<Index> others;
    for (std::size_t g = 0; g != factors_.size(); ++g) {
      if (g == f) continue;
      others.insert(factor_indices[g].begin(), factor_indices[g].end());
    }
    container::set<Index> result;
    if (named_indices.deduce()) {
      // only sums and products have dummies that could clash; if the other
      // factors only have indices of the subfactor, they cannot clash with
      // them: keep deducing
      if (!factors_[f].is<Sum>() && !factors_[f].is<Product>())
        return named_indices;
      if (std::ranges::all_of(others, [&](const Index &idx) {
            return factor_indices[f].contains(idx);
          }))
        return named_indices;
      // else name the indices that occur once in the subfactor, i.e. those
      // that the deduction would name
      const auto unique = get_unique_indices(factors_[f]);
      result.insert(unique.bra.begin(), unique.bra.end());
      result.insert(unique.ket.begin(), unique.ket.end());
      result.insert(unique.aux.begin(), unique.aux.end());
    } else {
      result = named_indices.indices();
    }
    result.insert(others.begin(), others.end());
    return detail::explicit_named_indices(std::move(result));
  };
  for (std::size_t f = 0; f != factors_.size(); ++f) {
    auto &factor = factors_[f];
    if (factor.template is<AbstractTensor>()) {
      continue;
    }
    auto bp = factor->canonicalize(subfactor_named_indices(f), opts);
    // the canonicalization may have renamed the dummies of the factor
    if (!factor_indices.empty()) factor_indices[f] = get_used_indices(factor);
    if (bp) {
      SEQUANT_ASSERT(bp->template is<Constant>());
      this->scalar_ *= std::static_pointer_cast<Constant>(bp)->value();
    }
  }

  if (Logger::instance().canonicalize) {
    std::wcout << "Product canonicalization(" << to_wstring(opts.method)
               << ") input: " << to_latex() << std::endl;
  }

  // pull out all scalar factors to the front
  auto is_scalar = [](const auto &factor) { return factor->is_scalar(); };
  auto scalars =
      factors_ | ranges::views::filter(is_scalar) | ranges::to_vector;
  // scalars commute, so we can reorder them freely
  ranges::sort(scalars, [](const auto &first, const auto &second) {
    return *first < *second;
  });

  factors_ = factors_ | ranges::views::filter([&is_scalar](const auto &factor) {
               return !is_scalar(factor);
             }) |
             ranges::to<decltype(factors_)>;

  // if there are no factors, insert scalars back and return
  if (factors_.empty()) {
    factors_.insert(factors_.begin(), scalars.begin(), scalars.end());
    return {};
  }

  auto contains_nontensors = ranges::any_of(factors_, [](const auto &factor) {
    return std::dynamic_pointer_cast<AbstractTensor>(factor) == nullptr;
  });
  if (!contains_nontensors) {  // tensor network canonization is a special case
                               // that's done in
                               // TensorNetwork
    TensorNetwork tn(this->factors_);
    const ExprPtr canon_factor = tn.canonicalize(
        detail::cardinal_tensor_labels(opts), named_indices, opts);

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
    if (canon_factor) {
      const auto &factor = canon_factor->template as<Constant>();
      scalar_ *= factor.value();
      // a network that vanishes by symmetry is not canonicalized; its tensors
      // would make the zero product print, hash and compare as nonzero ones
      if (factor.is_zero()) factors_.clear();
    }
    this->reset_hash_value();
  } else {  // if contains non-tensors, do commutation-checking resort
    // the scalars were pulled out above
    sort_nonscalar_factors(detail::cardinal_tensor_labels(opts));
  }
  // reinsert scalar factors at the front
  factors_.insert(factors_.begin(), scalars.begin(), scalars.end());

  // TODO evaluate product of Tensors (turn this into Products of Products)

  if (Logger::instance().canonicalize)
    std::wcout << "Product canonicalization(" << to_wstring(opts.method)
               << ") result: " << to_latex() << std::endl;

  return {};  // side effects are absorbed into the scalar_
}

namespace {

/// @return a key of @p tensor that does not depend on its indices
std::size_t nested_leaf_key(const AbstractTensor &tensor) {
  auto key = hash::value(tensor._label());
  hash::combine(key, tensor._bra_rank());
  hash::combine(key, tensor._ket_rank());
  hash::combine(key, tensor._aux_rank());
  hash::combine(key, static_cast<int>(tensor._symmetry()));
  hash::combine(key, static_cast<int>(tensor._braket_symmetry()));
  hash::combine(key, static_cast<int>(tensor._column_symmetry()));
  return key;
}

/// @return a key of the nested product @p product (see
/// Product::collect_nested_tensor_leaves()) that depends neither on the
/// indices nor on the order of its factors
/// @note the scalar and the scalar factors (e.g. Variable objects) are part of
/// the key: nested products that differ only by them are not interchangeable
std::size_t nested_product_key(const Product &product) {
  container::svector<std::size_t> keys;
  for (const auto &factor : product.factors()) {
    if (factor.is<AbstractTensor>())
      keys.push_back(nested_leaf_key(factor.as<AbstractTensor>()));
    else if (factor.is<Product>())
      keys.push_back(nested_product_key(factor.as<Product>()));
    else
      keys.push_back(factor->hash_value());
  }
  std::ranges::sort(keys);
  auto key = hash::range(keys.begin(), keys.end());
  hash::combine(key, Constant(product.scalar()).hash_value());
  return key;
}

/// collects the groups of the tensors of @p product (see
/// TensorNetworkV3::set_tensor_groups()): each nested product is a group, keyed
/// by its structure (see nested_product_key()), within the group of the
/// product that contains it; the tensors of the outermost product are in no
/// group
/// @param group the group of the tensor factors of @p product
/// @param groups the groups, to which those of the nested products are
/// appended
/// @param tensor_groups the group of each tensor, appended in the order of
/// Product::collect_nested_tensor_leaves()
void collect_nested_tensor_groups(
    const Product &product, const std::optional<std::size_t> &group,
    container::svector<TensorNetwork::TensorGroup> &groups,
    container::svector<std::optional<std::size_t>> &tensor_groups) {
  for (const auto &factor : product.factors()) {
    if (factor.is<AbstractTensor>()) {
      tensor_groups.push_back(group);
    } else if (factor.is<Product>()) {
      groups.push_back(
          {.key = nested_product_key(factor.as<Product>()), .parent = group});
      collect_nested_tensor_groups(factor.as<Product>(), groups.size() - 1,
                                   groups, tensor_groups);
    }
  }
}

}  // namespace

bool Product::collect_nested_tensor_leaves(
    container::svector<ExprPtr *> &leaves) {
  // N.B. mutable iteration resets the memoized hash of each product visited
  auto collect = [&leaves](auto &self, Product &product, bool &nested) -> bool {
    for (auto &factor : product) {
      if (factor.is<AbstractTensor>()) {
        leaves.push_back(&factor);
      } else if (factor.is<Product>()) {
        nested = true;
        if (!self(self, factor.as<Product>(), nested)) return false;
      } else if (!factor->is_scalar()) {
        return false;
      }
    }
    return true;
  };
  bool nested = false;
  return collect(collect, *this, nested) && nested;
}

void Product::canonicalize_nested(const container::svector<ExprPtr *> &leaves,
                                  const NamedIndexSpec &named_indices,
                                  const CanonicalizeOptions &opts) {
  if (Logger::instance().canonicalize) {
    std::wcout << "Product canonicalization(" << to_wstring(opts.method)
               << ") of nested products, input: " << to_latex() << std::endl;
  }

  const auto &cardinal_tensor_labels = detail::cardinal_tensor_labels(opts);
  TensorNetwork tn(leaves | ranges::views::transform([](const ExprPtr *leaf) {
                     return *leaf;
                   }) |
                   ranges::to<container::svector<ExprPtr>>);
  // the parenthesization breaks symmetries of the network: only the
  // relabelings that keep the tensors of each nested product together are
  // considered, so that the result does not depend on the input
  {
    container::svector<TensorNetwork::TensorGroup> groups;
    container::svector<std::optional<std::size_t>> tensor_groups;
    collect_nested_tensor_groups(*this, std::nullopt, groups, tensor_groups);
    tn.set_tensor_groups(std::move(groups), std::move(tensor_groups));
  }
  const ExprPtr byproduct =
      tn.canonicalize(cardinal_tensor_labels, named_indices, opts);
  if (byproduct) {
    const auto &factor = byproduct->as<Constant>();
    scalar_ *= factor.value();
    if (factor.is_zero()) {
      // see the tensor network branch of canonicalize_impl()
      factors_.clear();
      reset_hash_value();
      return;
    }
  }

  // put each canonicalized tensor where its input was: the tensor network
  // reorders its tensors, but the parenthesization must be kept
  const auto &tensors = tn.tensors();
  const auto &input_ordinals = tn.tensor_input_ordinals();
  SEQUANT_ASSERT(tensors.size() == leaves.size());
  for (std::size_t t = 0; t != tensors.size(); ++t) {
    auto exprptr = std::dynamic_pointer_cast<Expr>(tensors[t]);
    SEQUANT_ASSERT(exprptr);
    *leaves.at(input_ordinals.at(t)) = std::move(exprptr);
  }

  sort_factors(cardinal_tensor_labels, /* recursive = */ true);

  if (Logger::instance().canonicalize)
    std::wcout << "Product canonicalization(" << to_wstring(opts.method)
               << ") of nested products, result: " << to_latex() << std::endl;
}

void Product::sort_factors(
    const container::vector<std::wstring> &cardinal_tensor_labels,
    bool recursive) {
  if (recursive) {
    for (auto &factor : factors_) {
      if (factor.is<Product>())
        factor.as<Product>().sort_factors(cardinal_tensor_labels, true);
    }
  }

  // scalars commute, so they go to the front, in order
  auto is_scalar = [](const auto &factor) { return factor->is_scalar(); };
  auto scalars =
      factors_ | ranges::views::filter(is_scalar) | ranges::to_vector;
  ranges::sort(scalars, [](const auto &first, const auto &second) {
    return *first < *second;
  });
  factors_ = factors_ | ranges::views::filter([&is_scalar](const auto &factor) {
               return !is_scalar(factor);
             }) |
             ranges::to<decltype(factors_)>;

  sort_nonscalar_factors(cardinal_tensor_labels);

  factors_.insert(factors_.begin(), scalars.begin(), scalars.end());
  reset_hash_value();
}

void Product::sort_nonscalar_factors(
    const container::vector<std::wstring> &cardinal_tensor_labels) {
  // comparer that respects cardinal tensor labels
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
          cardinal_tensor_labels,
          [&second_label](const std::wstring &l) { return l == second_label; });
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
    bubble_sort(begin(factors_), end(factors_),
                [&local_compare](const ExprPtr &first, const ExprPtr &second) {
                  bool result = (first->commutes_with(*second))
                                    ? local_compare(first, second)
                                    : false;
                  return result;
                });
  }
  reset_hash_value();
}

void Product::adjoint() {
  SEQUANT_ASSERT(static_commutativity() == false);  // assert no slicing
  auto adj_scalar = conj(scalar());
  using namespace ranges;
  auto adj_factors =
      factors() | views::reverse |
      views::transform([](auto &expr) { return ::sequant::adjoint(expr); });
  *this =
      Product(adj_scalar, ranges::begin(adj_factors), ranges::end(adj_factors));
}

ExprPtr Product::canonicalize(const NamedIndexSpec &named_indices,
                              CanonicalizeOptions opt) {
  if (this->is_canonical(named_indices, opt)) return {};
  const auto contexts_version = current_contexts_version();
  auto byproduct = this->canonicalize_impl(named_indices, opt);
  this->mark_canonical(named_indices, opt, contexts_version);
  return byproduct;
}

ExprPtr Product::rapid_canonicalize(const NamedIndexSpec &named_indices,
                                    CanonicalizeOptions opt) {
  SEQUANT_ASSERT(opt.method == CanonicalizationMethod::Rapid);
  return this->canonicalize_impl(named_indices, opt);
}

bool Product::static_commutativity() const { return false; }

std::wstring Product::to_latex() const { return to_latex(false); }

std::wstring Product::to_latex(bool negate) const {
  std::wstring result;
  result = L"{";
  if (!scalar().is_zero()) {
    const auto scal = negate ? -scalar() : scalar();
    if (!scal.is_identity()) {
      // replace -1 prefactor by -
      if (!(negate ? scalar() : -scalar()).is_identity()) {
        result += io::latex::to_string(scal);
      } else {
        result += L"{-}";
      }
    }
    for (const auto &i : factors()) {
      if (i->is<Product>())
        result += L"\\bigl(" + i->to_latex() + L"\\bigr)";
      else
        result += i->to_latex();
    }
  }
  result += L"}";
  return result;
}

Product::type_id_type Product::type_id() const {
  return get_type_id<Product>();
}

/// @return an identical clone of this Product (a deep copy allocated on the
///         heap)
/// @note this does not flatten the product
ExprPtr Product::clone() const { return ex<Product>(this->deep_copy()); }

Product Product::deep_copy() const {
  Product result(this->scalar(), ExprPtrList{});
  // append() clones each factor; Flatten::No preserves nested products.
  for (const auto &factor : factors()) {
    result.append(1, factor, Flatten::No);
  }
  result.copy_canonical_mark(*this);
  return result;
}

Product &Product::operator*=(const Expr &that) {
  if (!that.is<Constant>()) {
    this->append(1, const_cast<Expr &>(that).shared_from_this());
  } else {
    scalar_ *= that.as<Constant>().value();
    reset_canonical_mark();
  }
  return *this;
}

void Product::add_identical(const Product &other) {
  SEQUANT_ASSERT(ranges::equal(this->factors(), other.factors()));
  scalar_ += other.scalar_;
  reset_canonical_mark();
}

void Product::add_identical(const std::shared_ptr<Product> &other) {
  SEQUANT_ASSERT(ranges::equal(this->factors(), other->factors()));
  scalar_ += other->scalar_;
  reset_canonical_mark();
}

void Product::add_identical(const ExprPtr &other) {
  if (other.is<Product>()) return this->add_identical(other.as<Product>());

  // only makes sense if this has a single factor
  SEQUANT_ASSERT(this->factors_.size() == 1 && this->factors_[0] == other);
  scalar_ += 1;
  reset_canonical_mark();
}

ExprIterator Product::begin_subexpr() {
  // N.B. the canonical mark need not be reset: changes of the factors are
  // detected by Expr::is_canonical()
  if (!factors_.empty()) {
    hash_value_.reset();
  }

  return ExprIterator{factors_.data()};
}

ExprIterator Product::end_subexpr() {
  // N.B. handing out a mutable iterator into factors_ invalidates the
  // memoized hash, regardless of which end of the range it points at
  // (`*(--end())` mutates just as `*begin()` does)
  if (!factors_.empty()) {
    hash_value_.reset();
  }

  return ExprIterator{factors_.data() + factors_.size()};
}

ConstExprIterator Product::begin_subexpr() const {
  return ConstExprIterator{factors_.data()};
}

ConstExprIterator Product::end_subexpr() const {
  return ConstExprIterator{factors_.data() + factors_.size()};
}

/// @return the hash of this object, by hashing only the factors,
/// not the scalar to make possible rapid finding of Products that only
/// differ by a factor
/// @note this ensures that hash of a Product involving a single factor is
/// identical to the hash of the factor itself.
Expr::hash_type Product::memoizing_hash() const {
  auto compute_hash = [this]() {
    if (factors().size() == 1)
      return factors_[0]->hash_value();
    else {
      auto deref_factors =
          factors() |
          ranges::views::transform(
              [](const ExprPtr &ptr) -> const Expr & { return *ptr; });
      auto value =
          hash::range(ranges::begin(deref_factors), ranges::end(deref_factors));
      return value;
    }
  };

  if (!hash_value_) {
    hash_value_ = compute_hash();
  } else {
    SEQUANT_ASSERT(*hash_value_ == compute_hash());
  }

  return *hash_value_;
}

bool Product::static_equal(const Expr &that) const {
  const auto &that_cast = static_cast<const Product &>(that);
  if (scalar() == that_cast.scalar() &&
      factors().size() == that_cast.factors().size()) {
    if (this->empty()) return true;
    // compare hash values first
    if (this->hash_value() ==
        that.hash_value())  // hash values agree -> do full comparison
      return std::equal(begin_subexpr(), end_subexpr(), that.begin_subexpr(),
                        expr_ptr_comparer);
    else
      return false;
  } else
    return false;
}

std::weak_ordering Product::static_lexicographic_compare(
    const Expr &that, const LexicographicCompareOptions &opts) const {
  const auto &that_cast = static_cast<const Product &>(that);
  if (auto cmp = std::lexicographical_compare_three_way(
          factors_.begin(), factors_.end(), that_cast.factors_.begin(),
          that_cast.factors_.end(),
          [&opts](const ExprPtr &a, const ExprPtr &b) {
            return a->lexicographic_compare(*b, opts);
          });
      cmp != 0)
    return cmp;
  if (auto cmp = std::compare_weak_order_fallback(scalar_.real(),
                                                  that_cast.scalar_.real());
      cmp != 0)
    return cmp;
  return std::compare_weak_order_fallback(scalar_.imag(),
                                          that_cast.scalar_.imag());
}

CProduct::CProduct(const Product &other) : Product(other) {}
CProduct::CProduct(Product &&other) : Product(std::move(other)) {}

bool CProduct::is_commutative() const { return true; }

ExprPtr CProduct::clone() const { return ex<CProduct>(this->deep_copy()); }

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

bool CProduct::static_commutativity() const { return true; }

NCProduct::NCProduct(const Product &other) : Product(other) {}
NCProduct::NCProduct(Product &&other) : Product(std::move(other)) {}

bool NCProduct::is_commutative() const { return false; }

ExprPtr NCProduct::clone() const { return ex<NCProduct>(this->deep_copy()); }

void NCProduct::adjoint() {
  auto adj_scalar = conj(scalar());
  using namespace ranges;
  // factors must be reversed since they do not commute
  auto adj_factors =
      factors() | views::reverse |
      views::transform([](auto &&expr) { return ::sequant::adjoint(expr); });
  *this = NCProduct(adj_scalar, ranges::begin(adj_factors),
                    ranges::end(adj_factors));
}

bool NCProduct::static_commutativity() const { return true; }

}  // namespace sequant
