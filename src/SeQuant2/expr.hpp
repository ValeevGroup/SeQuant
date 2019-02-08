//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_EXPR_HPP
#define SEQUANT2_EXPR_HPP

#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include <range/v3/all.hpp>

#include <boost/functional/hash.hpp>
#include <boost/callable_traits.hpp>

#include "latex.hpp"
#include "wolfram.hpp"
#include "vector.hpp"
#include "hash.hpp"

namespace sequant2 {

extern bool debug_canonicalize;

/// @brief Base expression class

/// Expr represents the interface needed to form expression trees. Classes that
/// represent expressions should publicly derive from this class. Each Expr on a tree has links
/// to its children Expr objects. The lifetime of Expr objects is expected to be managed by std::shared_ptr .
/// Expr is an Iterable over subexpressions (each of which is an Expr itself). More precisely,
/// Expr meets the SizedIterable concept (see https://raw.githubusercontent.com/ericniebler/range-v3/master/doc/std/D4128.md).
/// Specifically, iterators to subexpressions
/// dereference to std::shared_ptr<Expr>. Since Expr is a range, it provides begin/end/etc. that can participate in overloads
///       with other functions in the derived class. Consider a Container class derived from a BaseContainer class:
/// @code
///   template <typename T> class Container : public BaseContainer, public Expr {
///     // WARNING: BaseContainer::begin clashes with Expr::begin
///     // WARNING: BaseContainer::end clashes with Expr::end
///     // etc.
///   };
/// @endcode
/// There are two possible scenarios:
///   - if @c Container is a container of Expr objects, BaseContainer will iterate over std::shared_ptr<Expr> objects already
///     and both ranges will be equivalent; it is sufficient to add `using BaseContainer::begin`, etc. to Container's public API.
///   - if @c Container is a container of non-Expr objects, iteration over BaseContainer is likely to be more commonly used
///     in practice, hence again adding `using BaseContainer::begin`, etc. will suffice. To be able to iterate over subexpression
///     range (in this case it is empty) Expr provides Expr::expr member to cast to Expr:
/// @code
///    Container c(...);
///    for(const auto& e: c) {  // iterates over elements of BaseContainer
///    }
///    for(const auto& e: c.expr()) {  // iterates over subexpressions
///    }
/// @endcode
class Expr : public std::enable_shared_from_this<Expr>, public ranges::view_facade<Expr> {
 public:
  using range_type = ranges::view_facade<Expr>;
  using hash_type = std::size_t;
  using type_id_type = int;   // to speed up comparisons

  Expr() = default;
  virtual ~Expr() = default;

  /// @return true if this is a leaf
  bool is_atom() const {
    return ranges::empty(*this);
  }

  /// @return the string representation of @c this
  virtual std::wstring to_latex() const {
    throw std::logic_error("Expr::to_latex not implemented in this derived class");
  }

  /// @return the string representation of @c this in the Wolfram Language format
  virtual std::wstring to_wolfram() const {
    throw std::logic_error("Expr::to_wolfram not implemented in this derived class");
  }

  /// @return a clone of this object
  /// @note must be overridden in the derived class
  virtual std::shared_ptr<Expr> clone() const {
    throw std::logic_error("Expr::clone not implemented in this derived class");
  }

  /// Canonicalizes @c this and returns the biproduct of canonicalization (e.g. phase)
  /// @return the biproduct of canonicalization, or @c nullptr if no biproduct generated
  virtual std::shared_ptr<Expr> canonicalize() {
    return {};  // by default do nothing and return nullptr
  }

  /// Performs approximate, but fast, canonicalization of @c this and returns the biproduct of canonicalization (e.g. phase)
  /// The default is to use canonicalize(), unless overridden in the derived class.
  /// @return the biproduct of canonicalization, or @c nullptr if no biproduct generated
  virtual std::shared_ptr<Expr> rapid_canonicalize() {
    return this->canonicalize();
  }

  /// recursively visit the tree, i.e. call visitor on each subexpression in depth-first fashion
  /// @tparam Visitor a callable with (std::shared_ptr<Expr>&) or (const std::shared_ptr<Expr>&) signature
  /// @param visitor the visitor object
  /// @param atoms_only if true, will visit only the leafs; the default is to visit all nodes
  /// @return true if this object was visited
  /// @sa expr_range
  template <typename Visitor> bool visit(Visitor& visitor, const bool atoms_only = false) {
    for(auto& subexpr_ptr: expr()) {
      const auto subexpr_is_an_atom = subexpr_ptr->is_atom();
      const auto need_to_visit_subexpr = !atoms_only || subexpr_is_an_atom;
      bool visited = false;
      if (!subexpr_is_an_atom)  // if not a leaf, recur into it
        visited = subexpr_ptr->visit(visitor, atoms_only);
      // call on the subexpression itself, if not yet done so
      if (need_to_visit_subexpr && !visited)
        visitor(subexpr_ptr);
    }
    // can only visit itself here if visitor(const ExprPtr&) is valid
    bool this_visited = false;
    if constexpr(boost::callable_traits::is_invocable_r<void,Visitor,const std::shared_ptr<Expr>&>::value) {
      if (!atoms_only || this->is_atom()) {
        visitor(shared_from_this());
        this_visited = true;
      }
    }
    return this_visited;
  }

  auto begin_subexpr() {
    return range_type::begin();
  }

  auto end_subexpr() {
    return range_type::end();
  }

  auto begin_subexpr() const {
    return range_type::begin();
  }

  auto end_subexpr() const {
    return range_type::end();
  }

  Expr& expr() {
    return *this;
  }

  template <typename T, typename Enabler = void> struct is_shared_ptr_of_expr : std::false_type {};
  template <typename T>
  struct is_shared_ptr_of_expr<std::shared_ptr<T>, std::enable_if_t<std::is_same_v<Expr,T>> > : std::true_type {};
  template <typename T, typename Enabler = void> struct is_shared_ptr_of_expr_or_derived : std::false_type {};
  template <typename T>
  struct is_shared_ptr_of_expr_or_derived<std::shared_ptr<T>, std::enable_if_t<std::is_base_of<Expr,T>::value> > : std::true_type {};

  /// Computes and returns the memoized hash value.
  /// @note always returns 0 unless this derived class overrides Expr::memoizing_hash
  /// @return the hash value for this Expr
  hash_type hash_value() const {
    return memoizing_hash();
  };

  /// Computes and returns the derived type identifier
  /// @note this function must be overridden in the derived class
  /// @sa Expr::get_type_id
  /// @return the hash value for this Expr
  virtual type_id_type type_id() const
#if __GNUG__
  { abort(); }
#else
  =0;
#endif

  /// @tparam T Expr or a class derived from Expr
  /// @return true if @c *this is equal to @c that
  /// @note the derived class must implement Expr::static_equal
  template <typename T>
  std::enable_if_t<std::is_base_of<Expr,T>::value, bool>
      operator==(const T& that) const {
    if (this->type_id() != that.type_id())
      return false;
    else
      return this->static_equal(static_cast<const Expr &>(that));
  }

  /// @tparam T Expr or a class derived from Expr
  /// @return true if @c *this is equal to @c that
  /// @note the derived class must implement Expr::static_equal
  template <typename T>
  std::enable_if_t<std::is_base_of<Expr,T>::value, bool>
  operator!=(const T& that) const {
     return ! operator==(that);
  }

  /// @tparam T Expr or a class derived from Expr
  /// @return true if @c *this is less than @c that
  /// @note the derived class must implement Expr::static_less_than
  template<typename T>
  std::enable_if_t<std::is_base_of<Expr, T>::value, bool>
  operator<(const T &that) const {
    if (type_id() == that.type_id()) { // if same type, use generic (or type-specific, if available) comparison
      return static_less_than(static_cast<const Expr &>(that));
    } else {  // order types by type id
      return type_id() < that.type_id();
    }
  }
  /// @return (unique) type id of class T
  template <typename T> static  type_id_type get_type_id() {
    static type_id_type type_id = get_next_type_id();
    return type_id;
  };

  /// @tparam T an Expr type
  /// @return true if this object is of type @c T
  template<typename T>
  bool is() const { return this->type_id() == get_type_id<T>(); }

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template<typename T>
  const T &as() const {
    assert(this->is<T>());
    return static_cast<const T &>(*this);
  }

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template<typename T>
  T &as() {
    assert(this->is<T>());
    return static_cast<T &>(*this);
  }

/** @name in-place arithmetic operators
 *  Virtual in-place arithmetic operators to be overridden in expressions for which these make sense.
 */
///@{

  /// @brief in-place multiply @c *this by @c that
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be implemented for the particular @c that
  virtual Expr &operator*=(const Expr &that) {
    throw std::logic_error("Expr::operator*= not implemented in this derived class");
  }

  /// @brief in-place add @c that to @c *this
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be implemented for the particular @c that
  virtual Expr &operator+=(const Expr &that) {
    throw std::logic_error("Expr::operator+= not implemented in this derived class");
  }

  /// @brief in-place subtract @c that from @c *this
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be implemented for the particular @c that
  virtual Expr &operator-=(const Expr &that) {
    throw std::logic_error("Expr::operator-= not implemented in this derived class");
  }

///@}

 private:
  friend ranges::range_access;

 protected:
  Expr(Expr &&) = default;
  Expr(const Expr &) = default;
  Expr &operator=(Expr &&) = default;
  Expr &operator=(const Expr &) = default;

  struct cursor
  {
    using value_type = std::shared_ptr<Expr>;

    cursor() = default;
    constexpr explicit cursor(std::shared_ptr<Expr>* subexpr_ptr) noexcept
        : ptr_{subexpr_ptr}
    {}
    /// when take const ptr note runtime const flag
    constexpr explicit cursor(const std::shared_ptr<Expr>* subexpr_ptr) noexcept
        : ptr_{const_cast<std::shared_ptr<Expr>*>(subexpr_ptr)}, const_{true}
    {}
    bool equal(const cursor& that) const
    {
      return ptr_ == that.ptr_;
    }
    void next()
    {
      ++ptr_;
    }
    void prev()
    {
      --ptr_;
    }
    // TODO figure out why can't return const here if want to be able to assign to *begin(Expr&)
    std::shared_ptr<Expr>& read() const
    {
      RANGES_EXPECT(ptr_);
      return *ptr_;
    }
    std::shared_ptr<Expr>& read()
    {
      RANGES_EXPECT(const_ == false);
      RANGES_EXPECT(ptr_);
      return *ptr_;
    }
    void assign(const std::shared_ptr<Expr>& that_ptr)
    {
      RANGES_EXPECT(ptr_);
      *ptr_ = that_ptr;
    }
    std::ptrdiff_t distance_to(cursor const &that) const
    {
      return that.ptr_ - ptr_;
    }
    void advance(std::ptrdiff_t n)
    {
      ptr_ += n;
    }
   private:
    std::shared_ptr<Expr>* ptr_ = nullptr;  // both begin and end will be represented by this, so Expr without subexpressions begin() equals end() automatically
    bool const_ = false;  // assert in nonconst ops
  };

  /// @return the cursor for the beginning of the range (must override in a derived Expr that has subexpressions)
  virtual cursor begin_cursor()
  {
    return cursor{};
  }
  /// @return the cursor for the end of the range (must override in a derived Expr that has subexpressions)
  virtual cursor end_cursor()
  {
    return cursor{};
  }

  /// @return the cursor for the beginning of the range (must override in a derived Expr that has subexpressions)
  virtual cursor begin_cursor() const
  {
    return cursor{};
  }
  /// @return the cursor for the end of the range (must override in a derived Expr that has subexpressions)
  virtual cursor end_cursor() const
  {
    return cursor{};
  }

  mutable std::optional<hash_type> hash_value_;  // not initialized by default
  virtual hash_type memoizing_hash() const {
    static const hash_type default_hash_value = 0;
    if (hash_value_)
      return *hash_value_;
    else
      return default_hash_value;
  }
  void reset_hash_value() const {
    hash_value_.reset();
  }

  /// @param that an Expr object
  /// @note @c that is guaranteed to be of same type as @c *this, hence can be statically cast
  /// @return true if @c that is equivalent to *this
  virtual bool static_equal(const Expr &that) const
#if __GNUG__
  { abort(); }
#else
  =0;
#endif

  /// @param that an Expr object
  /// @note @c that is guaranteed to be of same type as @c *this, hence can be statically cast
  /// @note base comparison compares Expr::hash_value() , specialize to each type as needed
  /// @return true if @c *this is less than @c that
  virtual bool static_less_than(const Expr &that) const {
    return this->hash_value() < that.hash_value();
  }

 private:
  /// @return returns next type id in the grand class list
  static type_id_type get_next_type_id() {
    static std::atomic<type_id_type> grand_type_id = 0;
    return ++grand_type_id;
  };

};  // Expr

using ExprPtr = std::shared_ptr<Expr>;

/// make an ExprPtr to a new object of type T
/// @tparam T a class derived from Expr
/// @tparam Args a parameter pack type such that T(std::forward<Args>...) is well-formed
/// @param args a parameter pack such that T(args...) is well-formed
template <typename T, typename ... Args> ExprPtr make(Args&& ... args) {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

// this is needed when using std::make_shared<X>({ExprPtr,ExprPtr}), i.e. must std::make_shared<X>(ExprPtrList{ExprPtr,ExprPtr})
using ExprPtrList = std::initializer_list<ExprPtr>;
static auto expr_ptr_comparer = [](const auto& ptr1, const auto& ptr2) { return *ptr1 == *ptr2; };

class Constant : public Expr {
 public:
  Constant() = default;
  virtual ~Constant() = default;
  Constant(const Constant&) = default;
  Constant(Constant&&) = default;
  Constant& operator=(const Constant&) = default;
  Constant& operator=(Constant&&) = default;
  template <typename U> explicit Constant(U && value) : value_(std::forward<U>(value)) {}

  auto value() const { return value_; }

  std::wstring to_latex() const override {
    return L"{" + sequant2::to_latex(value()) + L"}";
  }

  std::wstring to_wolfram() const override {
    return sequant2::to_wolfram(value());
  }

  type_id_type type_id() const override {
    return get_type_id<Constant>();
  };

  std::shared_ptr<Expr> clone() const override {
    return make<Constant>(this->value());
  }

  virtual Expr &operator*=(const Expr &that) override {
    if (that.is<Constant>()) {
      value_ *= that.as<Constant>().value();
    } else {
      throw std::logic_error("Constant::operator*=(that): not valid for that");
    }
    return *this;
  }

  virtual Expr &operator+=(const Expr &that) override {
    if (that.is<Constant>()) {
      this->value() += that.as<Constant>().value();
    } else {
      throw std::logic_error("Constant::operator+=(that): not valid for that");
    }
    return *this;
  }

  virtual Expr &operator-=(const Expr &that) override {
    if (that.is<Constant>()) {
      this->value() -= that.as<Constant>().value();
    } else {
      throw std::logic_error("Constant::operator-=(that): not valid for that");
    }
    return *this;
  }

 private:
  std::complex<double> value_;

  hash_type memoizing_hash() const override {
    hash_value_ = boost::hash_value(value_);
    return *hash_value_;
  }

  bool static_equal(const Expr &that) const override {
    return value() == static_cast<const Constant&>(that).value();
  }

};

/// @brief generalized product, i.e. a scalar times a product of zero or more terms
class Product : public Expr {
 public:
  Product() = default;
  virtual ~Product() = default;
  Product(const Product&) = default;
  Product(Product&&) = default;
  Product& operator=(const Product&) = default;
  Product& operator=(Product&&) = default;

  /// construct a Product out of zero or more factors (multiplied by 1)
  /// @param factors the factors
  Product(ExprPtrList factors) : factors_(std::move(factors)) {}

  /// construct a Product out of zero or more factors multiplied by a scalar
  /// @tparam T a numeric type; it must be able to multiply std::complex<double>
  /// @param scalar a scalar of type T
  /// @param factors an initializer list of factors
  template<typename T>
  Product(T scalar, ExprPtrList factors) : scalar_(std::move(scalar)), factors_(std::move(factors)) {}

  /// construct a Product out of a range of factors
  /// @param begin the begin iterator
  /// @param end the end iterator
  template <typename Iterator> Product(Iterator begin, Iterator end) : factors_(begin, end) {}

  /// construct a Product out of a range of factors
  /// @tparam T a numeric type; it must be able to multiply std::complex<double>
  /// @param scalar a scalar of type T
  /// @param begin the begin iterator
  /// @param end the end iterator
  template <typename T, typename Iterator> Product(T scalar, Iterator begin, Iterator end) : scalar_(std::move(scalar)), factors_(begin, end) {}

  /// (post-)multiplies the product by @c scalar times @c factor
  template<typename T>
  Product &append(T scalar, ExprPtr factor) {
    scalar_ *= scalar;
    if (!factor->is<Product>()) {
      if (factor->is<Constant>()) {  // factor in Constant
        auto factor_constant = std::static_pointer_cast<Constant>(factor);
        scalar_ *= factor_constant->value();
        // no need to reset the hash since scalar is not hashed!
      } else {
        factors_.push_back(std::move(factor));
        reset_hash_value();
      }
    } else {  // factor is a product also ... flatten recursively
      auto factor_product = std::static_pointer_cast<Product>(factor);
      scalar_ *= factor_product->scalar_;
      for (auto &subfactor: *factor_product)
        this->append(1, subfactor);
//      using std::end;
//      using std::cbegin;
//      using std::cend;
//      factors_.insert(end(factors_), cbegin(factor_product->factors_), cend(factor_product->factors_));
    }
    return *this;
  }

  /// (pre-)multiplies the product by @c scalar times @c factor ; less efficient than append()
  template<typename T>
  Product &prepend(T scalar, ExprPtr factor) {
    scalar_ *= scalar;
    if (!factor->is<Product>()) {
      if (factor->is<Constant>()) {  // factor in Constant
        auto factor_constant = std::static_pointer_cast<Constant>(factor);
        scalar_ *= factor_constant->value();
        // no need to reset the hash since scalar is not hashed!
      } else {
        factors_.insert(factors_.begin(), std::move(factor));
        reset_hash_value();
      }
    } else {  // factor is a product also  ... flatten recursively
      auto factor_product = std::static_pointer_cast<Product>(factor);
      scalar_ *= factor_product->scalar_;
      for (auto &subfactor: *factor_product)
        this->prepend(1, subfactor);
//      using std::begin;
//      using std::cbegin;
//      using std::cend;
//      factors_.insert(begin(factors_), cbegin(factor_product->factors_), cend(factor_product->factors_));
    }
    return *this;
  }

  const std::complex<double> &scalar() const { return scalar_; }
  const auto& factors() const { return factors_; }
  auto &factors() { return factors_; }

  /// @return true if the number of factors is zero
  bool empty() const { return factors_.empty(); }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    if (scalar() != 0.) {
      if (scalar() != 1.) {
        result += sequant2::to_latex(scalar()) + L" \\times ";
      }
      for (const auto &i : factors())
        result += i->to_latex();
    }
    result += L"}";
    return result;
  }

  std::wstring to_wolfram() const override {
    std::wstring result = L"Times[";
    if (scalar() != decltype(scalar_)(1)) {
      result += sequant2::to_wolfram(scalar()) + L",";
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

  type_id_type type_id() const override {
    return get_type_id<Product>();
  };

  std::shared_ptr<Expr> clone() const override {
    auto cloned_factors =
        factors() | ranges::view::transform([](const ExprPtr &ptr) { return ptr ? ptr->clone() : nullptr; });
    return make<Product>(ranges::begin(cloned_factors), ranges::end(cloned_factors));
  }

  virtual Expr &operator*=(const Expr &that) override {
    this->append(1, const_cast<Expr &>(that).shared_from_this());
    return *this;
  }

  void add_identical(const std::shared_ptr<Product> &other) {
    assert(this->hash_value() == other->hash_value());
    scalar_ += other->scalar_;
  }

 private:
  std::complex<double> scalar_ = {1.0, 0.0};
  container::svector<ExprPtr> factors_{};

  cursor begin_cursor() override {
    return factors_.empty() ? Expr::begin_cursor() : cursor{&factors_[0]};
  };
  cursor end_cursor() override {
    return factors_.empty() ? Expr::end_cursor() : cursor{&factors_[0] + factors_.size()};
  };

  cursor begin_cursor() const override {
    return factors_.empty() ? Expr::begin_cursor() : cursor{&factors_[0]};
  };
  cursor end_cursor() const override {
    return factors_.empty() ? Expr::end_cursor() : cursor{&factors_[0] + factors_.size()};
  };

  /// @note this hashes only the factors, not the scalar to make possible rapid finding of identical factors
  hash_type memoizing_hash() const override {
    auto deref_factors = factors() | ranges::view::transform([](const ExprPtr &ptr) -> const Expr & { return *ptr; });
    hash_value_ = boost::hash_range(ranges::begin(deref_factors), ranges::end(deref_factors));
    return *hash_value_;
  }

  virtual std::shared_ptr<Expr> canonicalize() override;
  virtual std::shared_ptr<Expr> rapid_canonicalize() override;

  bool static_equal(const Expr &that) const override {
    const auto& that_cast = static_cast<const Product&>(that);
    if (scalar() == that_cast.scalar() && factors().size() == that_cast.factors().size()) {
      if (this->empty()) return true;
      // compare hash values first
      if (this->hash_value() == that.hash_value()) // hash values agree -> do full comparison
        return std::equal(begin_subexpr(), end_subexpr(), that.begin_subexpr(), expr_ptr_comparer);
      else
        return false;
    } else return false;
  }
};

/// @brief generalized sum of zero or more summands
class Sum : public Expr {
 public:
  Sum() = default;
  virtual ~Sum() = default;
  Sum(const Sum&) = default;
  Sum(Sum&&) = default;
  Sum& operator=(const Sum&) = default;
  Sum& operator=(Sum&&) = default;

  /// construct a Sum out of zero or more summands
  /// @param summands an initializer list of summands
  Sum(ExprPtrList summands) {
    // use append to flatten out Sum summands
    for(auto& summand: summands) {
      append(std::move(summand));
    }
  }

  /// construct a Sum out of a range of summands
  /// @param begin the begin iterator
  /// @param end the end iterator
  template<typename Iterator>
  Sum(Iterator begin, Iterator end) : summands_(begin, end) {}

  /// append a summand to the sum
  /// @param summand the summand
  Sum &append(ExprPtr summand) {
    if (!summand->is<Sum>()) {
      if (summand->is<Constant>()) {  // exclude zeros
        auto summand_constant = std::static_pointer_cast<Constant>(summand);
        if (summand_constant != 0) summands_.push_back(std::move(summand));
      } else {
        summands_.push_back(std::move(summand));
      }
      reset_hash_value();
    }
    else {  // this recursively flattens Sum summands
      for(auto& subsummand: *summand)
        this->append(subsummand);
    }
    return *this;
  }

  /// prepend a summand to the sum
  /// @param summand the summand
  Sum &prepend(ExprPtr summand) {
    if (!summand->is<Sum>()) {
      if (summand->is<Constant>()) {  // exclude zeros
        auto summand_constant = std::static_pointer_cast<Constant>(summand);
        if (summand_constant != 0) summands_.push_back(std::move(summand));
      } else {
        summands_.push_back(std::move(summand));
      }
      reset_hash_value();
    } else {  // this recursively flattens Sum summands
      for (auto &subsummand: *summand)
        this->prepend(subsummand);
    }
    return *this;
  }

  const auto& summands() const { return summands_; }

  /// @return true if the number of factors is zero
  bool empty() const { return summands_.empty(); }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{ \\left(";
    std::size_t counter = 0;
    for (const auto &i : summands()) {
      result += i->to_latex();
      ++counter;
      if (counter != summands().size())
        result += L" + ";
    }
    result += L"\\right) }";
    return result;
  }

  std::wstring to_wolfram() const override {
    std::wstring result;
    result = L"Plus[";
    std::size_t counter = 0;
    for (const auto &i : summands()) {
      result += i->to_wolfram();
      ++counter;
      if (counter != summands().size())
        result += L",";
    }
    result += L"]";
    return result;
  }

  Expr::type_id_type type_id() const override {
    return Expr::get_type_id<Sum>();
  };

  std::shared_ptr<Expr> clone() const override {
    auto cloned_summands = summands() | ranges::view::transform([](const ExprPtr &ptr) { return ptr->clone(); });
    return make<Sum>(ranges::begin(cloned_summands), ranges::end(cloned_summands));
  }

  virtual Expr &operator+=(const Expr &that) override {
    this->append(const_cast<Expr &>(that).shared_from_this());
    return *this;
  }

  virtual Expr &operator-=(const Expr &that) override {
    this->append(make<Product>(-1, ExprPtrList{const_cast<Expr &>(that).shared_from_this()}));
    return *this;
  }

 private:
  container::svector<ExprPtr> summands_{};

  cursor begin_cursor() override {
    return summands_.empty() ? Expr::begin_cursor() : cursor{&summands_[0]};
  };
  cursor end_cursor() override {
    return summands_.empty() ? Expr::end_cursor() : cursor{&summands_[0] + summands_.size()};
  };
  cursor begin_cursor() const override {
    return summands_.empty() ? Expr::begin_cursor() : cursor{&summands_[0]};
  };
  cursor end_cursor() const override {
    return summands_.empty() ? Expr::end_cursor() : cursor{&summands_[0] + summands_.size()};
  };

  hash_type memoizing_hash() const override {
    auto deref_summands = summands() | ranges::view::transform([](const ExprPtr &ptr) -> const Expr & { return *ptr; });
    hash_value_ = boost::hash_range(ranges::begin(deref_summands), ranges::end(deref_summands));
    return *hash_value_;
  }

  virtual std::shared_ptr<Expr> canonicalize() override {
    return canonicalize_<true>();
  }
  virtual std::shared_ptr<Expr> rapid_canonicalize() override {
    return canonicalize_<false>();
  }

  template<bool TwoPass>
  std::shared_ptr<Expr> canonicalize_() {

    const auto npasses = TwoPass ? 2 : 1;
    for (auto pass = 0; pass != npasses; ++pass) {
      // recursively canonicalize summands ...
      const auto nsubexpr = ranges::size(*this);
      for (std::size_t i = 0; i != nsubexpr; ++i) {
        auto bp = (pass == 0) ? summands_[i]->rapid_canonicalize() : summands_[i]->canonicalize();
        if (bp) {
          assert(bp->template is<Constant>());
          summands_[i] = make<Product>(std::static_pointer_cast<Constant>(bp)->value(), ExprPtrList{summands_[i]});
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

  bool static_equal(const Expr &that) const override {
    const auto& that_cast = static_cast<const Sum&>(that);
    if (summands().size() == that_cast.summands().size()) {
      if (this->empty()) return true;
      // compare hash values first
      if (this->hash_value() == that.hash_value()) // hash values agree -> do full comparison
        return std::equal(begin_subexpr(), end_subexpr(), that.begin_subexpr(), expr_ptr_comparer);
      else
        return false;
    } else return false;
  }

};

inline std::wstring to_latex(const Expr& expr) {
  return expr.to_latex();
}

inline std::wstring to_latex(const ExprPtr& exprptr) {
  return exprptr->to_latex();
}

/// splits long outer sum into a multiline align
inline std::wstring to_latex_align(const ExprPtr &exprptr, size_t max_terms_per_align = 0) {
  std::wstring result = to_latex(exprptr);
  if (exprptr->is<Sum>()) {
    result.erase(0, 7);  // remove leading  "{ \left"
    result.replace(result.size() - 9, 9, L")");  // replace trailing "\right) }" with ")"
    result = std::wstring(L"\\begin{align}\n") + result;
    // assume no inner sums
    int term_counter = 0;
    std::wstring::size_type pos = 0;
    while ((pos = result.find(L" + ", pos + 1)) != std::wstring::npos) {
      ++term_counter;
      if (max_terms_per_align > 0 && term_counter >= max_terms_per_align) {
        result.insert(pos + 3, L"\n\\end{align}\n\\begin{align}\n");
        term_counter = 1;
      } else {
        result.insert(pos + 3, L"\\\\\n");
      }
    }
  } else {
    result = std::wstring(L"\\begin{align}\n") + result;
  }
  result += L"\n\\end{align}";
  return result;
}

inline std::wstring to_wolfram(const Expr &expr) {
  return expr.to_wolfram();
}

inline std::wstring to_wolfram(const ExprPtr &exprptr) {
  return exprptr->to_wolfram();
}

template<typename Sequence>
std::decay_t<Sequence> clone(Sequence &&exprseq) {
  auto cloned_seq = exprseq | ranges::view::transform([](const ExprPtr &ptr) { return ptr ? ptr->clone() : nullptr; });
  return std::decay_t<Sequence>(ranges::begin(cloned_seq), ranges::end(cloned_seq));
}



};  // namespace sequant2

#include "expr_operator.hpp"
#include "expr_algorithm.hpp"

#endif //SEQUANT2_EXPR_HPP
