#ifndef SEQUANT_EXPRESSIONS_PRODUCT_HPP
#define SEQUANT_EXPRESSIONS_PRODUCT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/meta.hpp>

#include <cassert>
#include <string>
#include <type_traits>

namespace sequant {

/// @brief generalized product, i.e. a scalar times a product of zero or more
/// terms.
///
/// Product is distributive over addition (see Sum). It is associative.
/// All constructors and mutating methods (append/prepend) will by default
/// flatten its factors recursively (but this can be fully controlled
/// by specifying the flatten tag). The commutativity of factors is checked at
/// runtime for each factor (see CProduct and NCProduct for statically
/// commutative and noncommutative Product, respectively)
class Product : public Expr {
 public:
  enum class Flatten { Once, Recursively, Yes = Recursively, No };

  using scalar_type = Constant::scalar_type;

  Product() = default;
  virtual ~Product() = default;
  Product(const Product &) = default;
  Product(Product &&) = default;
  Product &operator=(const Product &) = default;
  Product &operator=(Product &&) = default;

  /// construct a Product out of zero or more factors (multiplied by 1)
  /// @param factors the factors
  /// @param flatten_tag if Flatten::Yes, flatten the factors
  Product(ExprPtrList factors, Flatten flatten_tag = Flatten::Yes) {
    using std::begin;
    using std::end;
    for (auto it = begin(factors); it != end(factors); ++it)
      append(1, *it, flatten_tag);
  }

  /// construct a Product out of zero or more factors (multiplied by 1)
  /// @param rng a range of factors
  /// @param flatten_tag if Flatten::Yes, flatten the factors
  template <typename Range,
            typename = std::enable_if_t<meta::is_range_v<std::decay_t<Range>> &&
                                        !meta::is_same_v<Range, ExprPtrList> &&
                                        !meta::is_same_v<Range, Product>>>
  explicit Product(Range &&rng, Flatten flatten_tag = Flatten::Yes) {
    using ranges::begin;
    using ranges::end;
    for (auto &&v : rng) append(1, std::forward<decltype(v)>(v), flatten_tag);
  }

  /// construct a Product out of zero or more factors (multiplied by 1)
  /// @tparam T a numeric type; it must be able to multiply Product::scalar_type
  /// @param scalar a scalar of type T
  /// @param rng a range of factors
  /// @param flatten_tag if Flatten::Yes, flatten the factors
  template <typename T, typename Range,
            typename = std::enable_if_t<
                meta::is_range_v<std::decay_t<Range>> &&
                !std::is_same_v<std::remove_reference_t<Range>, ExprPtrList> &&
                !std::is_same_v<std::remove_reference_t<Range>, Product>>>
  explicit Product(T scalar, Range &&rng, Flatten flatten_tag = Flatten::Yes)
      : scalar_(std::move(scalar)) {
    using ranges::begin;
    using ranges::end;
    for (auto &&v : rng) append(1, std::forward<decltype(v)>(v), flatten_tag);
  }

  /// construct a Product out of zero or more factors multiplied by a scalar
  /// @tparam T a numeric type; it must be able to multiply Product::scalar_type
  /// @param scalar a scalar of type T
  /// @param factors an initializer list of factors
  /// @param flatten_tag if Flatten::Yes, flatten the factors
  template <typename T>
  Product(T scalar, ExprPtrList factors, Flatten flatten_tag = Flatten::Yes)
      : scalar_(std::move(scalar)) {
    using std::begin;
    using std::end;
    for (auto it = begin(factors); it != end(factors); ++it)
      append(1, *it, flatten_tag);
  }

  /// construct a Product out of a range of factors
  /// @param begin the begin iterator
  /// @param end the end iterator
  /// @param flatten_tag specifies whether (and how) to flatten the argument(s)
  template <typename Iterator>
  Product(Iterator begin, Iterator end, Flatten flatten_tag = Flatten::Yes) {
    for (auto it = begin; it != end; ++it) append(1, *it, flatten_tag);
  }

  /// construct a Product out of a range of factors
  /// @tparam T a numeric type; it must be able to multiply Product::scalar_type
  /// @param scalar a scalar of type T
  /// @param begin the begin iterator
  /// @param end the end iterator
  /// @param flatten_tag specifies whether (and how) to flatten the argument(s)
  template <typename T, typename Iterator>
  Product(T scalar, Iterator begin, Iterator end,
          Flatten flatten_tag = Flatten::Yes)
      : scalar_(std::move(scalar)) {
    for (auto it = begin; it != end; ++it) append(1, *it, flatten_tag);
  }

  /// multiplies the product by @c scalar
  /// @param scalar a scalar by which to multiply the product
  /// @return @c *this
  template <typename T>
  Product &scale(T scalar) {
    scalar_ *= scalar;
    return *this;
  }

  /// (post-)multiplies the product by @c scalar times @c factor
  /// @param scalar a scalar by which to multiply the product
  /// @param factor a factor by which to multiply the product
  /// @param flatten_tag specifies whether (and how) to flatten the argument(s)
  /// @return @c *this
  template <typename T>
  Product &append(T scalar, ExprPtr factor,
                  Flatten flatten_tag = Flatten::Yes) {
    assert(factor);
    scalar_ *= scalar;
    if (!factor->is<Product>()) {
      if (factor->is<Constant>()) {  // factor in Constant
        auto factor_constant = factor->as<Constant>();
        scalar_ *= factor_constant.value();
        // no need to reset the hash since scalar is not hashed!
      } else {
        factors_.push_back(factor->clone());
        reset_hash_value();
      }
    } else {                             // factor is a product also ..
      if (flatten_tag != Flatten::No) {  // flatten, once or recursively
        const auto &factor_product = factor->as<Product>();
        scalar_ *= factor_product.scalar_;
        for (auto &&subfactor : factor_product)
          this->append(1, subfactor,
                       flatten_tag == Flatten::Once ? Flatten::No
                                                    : Flatten::Recursively);
      } else {
        factors_.push_back(factor->clone());
        reset_hash_value();
      }
    }
    return *this;
  }

  /// (post-)multiplies the product by @c scalar times @c factor
  /// @param scalar a scalar by which to multiply the product
  /// @param factor a factor by which to multiply the product
  /// @param flatten_tag specifies whether (and how) to flatten the argument(s)
  /// @return @c *this
  /// @warning if @p factor is a Product, it is flattened recursively
  template <typename T, typename Factor,
            typename = std::enable_if_t<is_an_expr_v<Factor>>>
  Product &append(T scalar, Factor &&factor,
                  Flatten flatten_tag = Flatten::Yes) {
    return this->append(scalar,
                        std::static_pointer_cast<Expr>(
                            std::forward<Factor>(factor).shared_from_this()),
                        flatten_tag);
  }

  /// (pre-)multiplies the product by @c scalar times @c factor
  /// @param scalar a scalar by which to multiply the product
  /// @param factor a factor by which to multiply the product
  /// @param flatten_tag specifies whether (and how) to flatten the argument(s)
  /// @return @c *this
  /// @warning if @p factor is a Product, it is flattened recursively
  /// @note this is less efficient than append()
  template <typename T>
  Product &prepend(T scalar, ExprPtr factor,
                   Flatten flatten_tag = Flatten::Yes) {
    assert(factor);
    scalar_ *= scalar;
    if (!factor->is<Product>()) {
      if (factor->is<Constant>()) {  // factor in Constant
        auto factor_constant = std::static_pointer_cast<Constant>(factor);
        scalar_ *= factor_constant->value();
        // no need to reset the hash since scalar is not hashed!
      } else {
        factors_.insert(factors_.begin(), factor->clone());
        reset_hash_value();
      }
    } else {  // factor is a product also  ... flatten recursively
      const auto &factor_product = factor->as<Product>();
      scalar_ *= factor_product.scalar_;
      if (flatten_tag != Flatten::No) {  // flatten, once or recursively
        for (auto &&subfactor : factor_product)
          this->prepend(1, subfactor,
                        flatten_tag == Flatten::Once ? Flatten::No
                                                     : Flatten::Recursively);
      } else {
        factors_.insert(factors_.begin(), factor->clone());
        reset_hash_value();
      }
    }
    return *this;
  }

  /// (pre-)multiplies the product by @c scalar times @c factor
  /// @param scalar a scalar by which to multiply the product
  /// @param factor a factor by which to multiply the product
  /// @param flatten_tag specifies whether (and how) to flatten the argument(s)
  /// @return @c *this
  /// @warning if @p factor is a Product, it is flattened recursively
  /// @note this is less efficient than append()
  template <typename T, typename Factor,
            typename = std::enable_if_t<is_an_expr_v<Factor>>>
  Product &prepend(T scalar, Factor &&factor,
                   Flatten flatten_tag = Flatten::Yes) {
    return this->prepend(scalar,
                         std::static_pointer_cast<Expr>(
                             std::forward<Factor>(factor).shared_from_this()),
                         flatten_tag);
  }

  const auto &scalar() const { return scalar_; }

  /// @return `Constant::is_zero(this->scalar())`
  bool is_zero() const { return Constant::is_zero(this->scalar()); }

  const auto &factors() const { return factors_; }
  auto &factors() { return factors_; }

  /// Factor accessor
  /// @param i factor index
  /// @return ith factor
  const ExprPtr &factor(size_t i) const { return factors_.at(i); }

  /// @return true if the number of factors is zero
  bool empty() const { return factors_.empty(); }

  /// @brief checks commutativity recursively
  /// @return true if definitely commutative, false definitely not commutative
  /// @note this is memoizing
  /// @sa CProduct::is_commutative() and NCProduct::is_commutative()
  virtual bool is_commutative() const;

  /// @brief adjoint of a Product is a reversed product of adjoints of its
  /// factors, with complex-conjugated scalar
  virtual void adjoint() override;

 private:
  /// @return true if commutativity is decidable statically
  /// @sa CProduct::static_commutativity() and NCProduct::static_commutativity()
  virtual bool static_commutativity() const { return false; }

 public:
  std::wstring to_latex() const override { return to_latex(false); }

  /// just like Expr::to_latex() , but can negate before conversion
  /// @param[in] negate if true, scalar will be before conversion
  std::wstring to_latex(bool negate) const {
    std::wstring result;
    result = L"{";
    if (!scalar().is_zero()) {
      const auto scal = negate ? -scalar() : scalar();
      if (!scal.is_identity()) {
        // replace -1 prefactor by -
        if (!(negate ? scalar() : -scalar()).is_identity()) {
          result += sequant::to_latex(scal);
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

  std::wstring to_wolfram() const override {
    std::wstring result =
        is_commutative() ? L"Times[" : L"NonCommutativeMultiply[";
    if (scalar() != decltype(scalar_)(1)) {
      result += sequant::to_wolfram(scalar()) + L",";
    }
    const auto nfactors = factors().size();
    size_t factor_count = 1;
    for (const auto &i : factors()) {
      result += i->to_wolfram() + (factor_count == nfactors ? L"" : L",");
      ++factor_count;
    }
    result += L"]";
    return result;
  }

  type_id_type type_id() const override { return get_type_id<Product>(); };

  /// @return an identical clone of this Product (a deep copy allocated on the
  ///         heap)
  /// @note this does not flatten the product
  ExprPtr clone() const override { return ex<Product>(this->deep_copy()); }

  Product deep_copy() const {
    auto cloned_factors =
        factors() | ranges::views::transform([](const ExprPtr &ptr) {
          return ptr ? ptr->clone() : nullptr;
        });
    Product result(this->scalar(), ExprPtrList{});
    ranges::for_each(cloned_factors, [&](const auto &cloned_factor) {
      result.append(1, std::move(cloned_factor), Flatten::No);
    });
    return result;
  }

  virtual Expr &operator*=(const Expr &that) override {
    if (!that.is<Constant>()) {
      this->append(1, const_cast<Expr &>(that).shared_from_this());
    } else {
      scalar_ *= that.as<Constant>().value();
    }
    return *this;
  }

  void add_identical(const Product &other) {
    assert(this->hash_value() == other.hash_value());
    scalar_ += other.scalar_;
  }

  void add_identical(const std::shared_ptr<Product> &other) {
    assert(this->hash_value() == other->hash_value());
    scalar_ += other->scalar_;
  }

  void add_identical(const ExprPtr &other) {
    if (other.is<Product>()) return this->add_identical(other.as<Product>());

    // only makes sense if this has a single factor
    assert(this->factors_.size() == 1 &&
           this->factors_[0]->hash_value() == other->hash_value());
    scalar_ += 1;
  }

 private:
  scalar_type scalar_ = {1, 0};
  container::svector<ExprPtr, 2> factors_{};

  cursor begin_cursor() override {
    return factors_.empty() ? Expr::begin_cursor() : cursor{&factors_[0]};
  };
  cursor end_cursor() override {
    return factors_.empty() ? Expr::end_cursor()
                            : cursor{&factors_[0] + factors_.size()};
  };

  cursor begin_cursor() const override {
    return factors_.empty() ? Expr::begin_cursor() : cursor{&factors_[0]};
  };
  cursor end_cursor() const override {
    return factors_.empty() ? Expr::end_cursor()
                            : cursor{&factors_[0] + factors_.size()};
  };

  /// @return the hash of this object, by hashing only the factors,
  /// not the scalar to make possible rapid finding of Products that only
  /// differ by a factor
  /// @note this ensures that hash of a Product involving a single factor is
  /// identical to the hash of the factor itself.
  hash_type memoizing_hash() const override {
    auto compute_hash = [this]() {
      if (factors().size() == 1)
        return factors_[0]->hash_value();
      else {
        auto deref_factors =
            factors() |
            ranges::views::transform(
                [](const ExprPtr &ptr) -> const Expr & { return *ptr; });
        auto value = hash::range(ranges::begin(deref_factors),
                                 ranges::end(deref_factors));
        return value;
      }
    };

    if (!hash_value_) {
      hash_value_ = compute_hash();
    } else {
      assert(*hash_value_ == compute_hash());
    }

    return *hash_value_;
  }

  ExprPtr canonicalize_impl(CanonicalizeOptions);
  virtual ExprPtr canonicalize(
      CanonicalizeOptions opt =
          CanonicalizeOptions::default_options()) override;
  virtual ExprPtr rapid_canonicalize(
      CanonicalizeOptions opts =
          CanonicalizeOptions::default_options().copy_and_set_method(
              CanonicalizationMethod::Rapid)) override;

  bool static_equal(const Expr &that) const override {
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
};  // class Product

class CProduct : public Product {
 public:
  using Product::Product;
  CProduct(const Product &other) : Product(other) {}
  CProduct(Product &&other) : Product(other) {}

  bool is_commutative() const override { return true; }

  /// @brief adjoint of a CProduct is a product of adjoints of its factors, with
  /// complex-conjugated scalar
  /// @note factors are not reversed since the factors commute
  virtual void adjoint() override;

 private:
  bool static_commutativity() const override { return true; }
};  // class CProduct

class NCProduct : public Product {
 public:
  using Product::Product;
  NCProduct(const Product &other) : Product(other) {}
  NCProduct(Product &&other) : Product(other) {}

  bool is_commutative() const override { return false; }

  /// @brief adjoint of a NCProduct is a reserved product of adjoints of its
  /// factors, with complex-conjugated scalar
  virtual void adjoint() override;

 private:
  bool static_commutativity() const override { return true; }
};  // class NCProduct

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_PRODUCT_HPP
