//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT_EXPR_HPP
#define SEQUANT_EXPR_HPP

#include <atomic>
#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include <range/v3/all.hpp>

#include <boost/core/demangle.hpp>
#include <boost/functional/hash.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "container.hpp"
#include "expr_fwd.hpp"
#include "hash.hpp"
#include "latex.hpp"
#include "meta.hpp"
#include "utility.hpp"
#include "wolfram.hpp"

namespace sequant {

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

  /// @return the string representation of @c this in the LaTeX format
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

  /// recursively visit this expression, i.e. call visitor on each subexpression
  /// in depth-first fashion.
  /// @warning this will only work for tree expressions; no checking is
  /// performed that each subexpression has only been visited once
  /// TODO make work for graphs
  /// @tparam Visitor a callable with (std::shared_ptr<Expr>&) or (const
  /// std::shared_ptr<Expr>&) signature
  /// @param visitor the visitor object
  /// @param atoms_only if true, will visit only the leaves; the default is to
  /// visit all nodes
  /// @return true if this object was visited
  /// @sa expr_range
  template<typename Visitor>
  bool visit(Visitor &&visitor, const bool atoms_only = false) {
    for(auto& subexpr_ptr: expr()) {
      const auto subexpr_is_an_atom = subexpr_ptr->is_atom();
      const auto need_to_visit_subexpr = !atoms_only || subexpr_is_an_atom;
      bool visited = false;
      if (!subexpr_is_an_atom)  // if not a leaf, recur into it
        visited = subexpr_ptr->visit(std::forward<Visitor>(visitor), atoms_only);
      // call on the subexpression itself, if not yet done so
      if (need_to_visit_subexpr && !visited)
        visitor(subexpr_ptr);
    }
    // can only visit itself here if visitor(const ExprPtr&) is valid
    bool this_visited = false;
    if constexpr(std::is_invocable_r_v<void,
                                       std::remove_reference_t<Visitor>,
                                       const std::shared_ptr<Expr> &>) {
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
  struct is_shared_ptr_of_expr_or_derived<
      std::shared_ptr<T>, std::enable_if_t<std::is_base_of<Expr, T>::value>>
      : std::true_type {};

  /// @brief Reports if this is a c-number
  /// (https://en.wikipedia.org/wiki/C-number), i.e. it commutes
  /// multiplicatively with c-numbers and q-numbers
  /// @return true if this is a c-number
  /// @warning this returns true for all leaves, hence must be overridden for
  /// leaf q-numbers
  /// @note for leaves this has O(1) cost, for non-leaves this involves checking
  /// subexpressions
  virtual bool is_cnumber() const {
    if (is_atom())
      return true;
    else {
      bool result = true;
      for (auto it = begin_subexpr(); result && it != end_subexpr(); ++it) {
        result &= (*it)->is_cnumber();
      }
      return result;
    }
  }

  /// @brief Checks if this commutes (wrt multiplication) with @c that
  /// @return true if this commutes with @c that
  /// @note the default implementation checks if either is c-number; if both are
  /// q-numbers
  ///       this checks commutativity of each subexpression with (each
  ///       subexpression of) @c that
  /// @note expressions are assumed to always commute with respect to additions
  /// since
  ///       it does not appear that +nonabelian near-rings
  ///       (https://en.wikipedia.org/wiki/Near-ring) are commonly needed.
  /// @note commutativity of leaves is checked by commutes_with_atom()
  bool commutes_with(const Expr &that) const {
    auto this_is_atom = is_atom();
    auto that_is_atom = that.is_atom();
    bool result = true;
    if (this_is_atom && that_is_atom) {
      result =
          this->is_cnumber() || that.is_cnumber() || commutes_with_atom(that);
    } else if (this_is_atom) {
      if (!this->is_cnumber()) {
        for (auto it = that.begin_subexpr(); result && it != that.end_subexpr();
             ++it) {
          result &= this->commutes_with(**it);
        }
      }
    } else {
      for (auto it = this->begin_subexpr(); result && it != this->end_subexpr();
           ++it) {
        result &= (*it)->commutes_with(that);
      }
    }
    return result;
  }

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
  template <typename T>
  static type_id_type get_type_id() {
    return type_id_accessor<T>();
  };

  /// sets (unique) type id of class T
  /// @param id the value of type id of class T
  /// @note since get_type_id does not check for duplicates, it's user's
  /// responsiblity to make sure that there are no collisions between type ids
  template <typename T>
  static void set_type_id(type_id_type id) {
    type_id_accessor<T>() = id;
  };

  /// @tparam T an Expr type
  /// @return true if this object is of type @c T
  template <typename T>
  bool is() const {
    return this->type_id() == get_type_id<std::decay_t<T>>();
  }

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template <typename T>
  const T &as() const {
    assert(this->is<std::decay_t<T>>());  // so that as<const T>() works fine
    return static_cast<const T &>(*this);
  }

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template <typename T>
  T &as() {
    assert(this->is<std::decay_t<T>>());  // so that as<const T>() works fine
    return static_cast<T &>(*this);
  }

  /// @return the (demangled) name of this type
  /// @note uses RTTI
  std::string type_name() const {
    return boost::core::demangle(typeid(*this).name());
  }

  /** @name in-place arithmetic operators
   *  Virtual in-place arithmetic operators to be overridden in expressions for
   * which these make sense.
   */
  ///@{

  /// @brief in-place multiply @c *this by @c that
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator*=(const Expr &that) {
    throw std::logic_error(
        "Expr::operator*= not implemented in this derived class");
  }

  /// @brief in-place non-commutatively-multiply @c *this by @c that
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator^=(const Expr &that) {
    throw std::logic_error(
        "Expr::operator^= not implemented in this derived class");
  }

  /// @brief in-place add @c that to @c *this
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
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
  virtual void reset_hash_value() const {
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
  /// @note @c that is guaranteed to be of same type as @c *this, hence can be
  /// statically cast
  /// @note base comparison compares Expr::hash_value() , specialize to each
  /// type as needed
  /// @return true if @c *this is less than @c that
  virtual bool static_less_than(const Expr &that) const {
    return this->hash_value() < that.hash_value();
  }

  /// @param that an Expr object
  /// @note @c *this and @c that are guaranteed to be leaves, and neither is a
  /// c-number , hence honest checking is needed
  /// @return true if @c *this multiplicatively commutes with @c that
  /// @note this returns true unless overridden in derived class
  virtual bool commutes_with_atom(const Expr& that) const {
    return true;
  }

 private:
  /// @return returns next type id in the grand class list
  static type_id_type get_next_type_id() {
    static std::atomic<type_id_type> grand_type_id = 0;
    return ++grand_type_id;
  };

  /// sets (unique) type id of class T
  /// @param id the value of type id of class T
  template <typename T>
  static type_id_type &type_id_accessor() {
    static type_id_type type_id = get_next_type_id();
    return type_id;
  };

};  // class Expr

/// make an ExprPtr to a new object of type T
/// @tparam T a class derived from Expr
/// @tparam Args a parameter pack type such that T(std::forward<Args>...) is well-formed
/// @param args a parameter pack such that T(args...) is well-formed
template<typename T, typename ... Args>
ExprPtr ex(Args &&... args) {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

// this is needed when using std::make_shared<X>({ExprPtr,ExprPtr}), i.e. must std::make_shared<X>(ExprPtrList{ExprPtr,ExprPtr})
using ExprPtrList = std::initializer_list<ExprPtr>;
static auto expr_ptr_comparer = [](const auto& ptr1, const auto& ptr2) { return *ptr1 == *ptr2; };

class Constant : public Expr {
 private:
  std::complex<double> value_;

 public:
  Constant() = default;
  virtual ~Constant() = default;
  Constant(const Constant&) = default;
  Constant(Constant&&) = default;
  Constant& operator=(const Constant&) = default;
  Constant& operator=(Constant&&) = default;
  template <typename U> explicit Constant(U && value) : value_(std::forward<U>(value)) {}

  /// @tparam T the result type; default to the type of value_
  /// @return the value cast to ResultType
  /// @throw std::invalid_argument if conversion to T is not possible
  /// @throw boost::numeric::positive_overflow or boost::numeric::negative_overflow if cast fails
  template <typename T = decltype(value_)>
  auto value() const {
    if constexpr (std::is_arithmetic_v<T>) {
      assert(value_.imag() == 0);
      return boost::numeric_cast<T>(value_.real());
    } else if constexpr (meta::is_complex_v<T>) {
      return T(boost::numeric_cast<typename T::value_type>(value_.real()),
               boost::numeric_cast<typename T::value_type>(value_.imag()));
    } else
      throw std::invalid_argument(
          "Constant::value<T>: cannot convert value to type T");
  }

  std::wstring to_latex() const override {
    return L"{" + sequant::to_latex(value()) + L"}";
  }

  std::wstring to_wolfram() const override {
    return sequant::to_wolfram(value());
  }

  type_id_type type_id() const override {
    return get_type_id<Constant>();
  }

  std::shared_ptr<Expr> clone() const override {
    return ex<Constant>(this->value());
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
      value_ += that.as<Constant>().value();
    } else {
      throw std::logic_error("Constant::operator+=(that): not valid for that");
    }
    return *this;
  }

  virtual Expr &operator-=(const Expr &that) override {
    if (that.is<Constant>()) {
      value_ -= that.as<Constant>().value();
    } else {
      throw std::logic_error("Constant::operator-=(that): not valid for that");
    }
    return *this;
  }

 private:
  hash_type memoizing_hash() const override {
    hash_value_ = boost::hash_value(value_);
    return *hash_value_;
  }

  bool static_equal(const Expr &that) const override {
    return value() == static_cast<const Constant &>(that).value();
  }
};  // class Constant

using ConstantPtr = std::shared_ptr<Constant>;

/// @brief generalized product, i.e. a scalar times a product of zero or more
/// terms.
///
/// Product is distributive over addition (see Sum). It is associative and is
/// flattened automatically. It's commutativity is checked at runtime for each
/// factor (see CProduct and NCProduct for statically commutative and
/// noncommutative Product, respectively)
class Product : public Expr {
 public:
  Product() = default;
  virtual ~Product() = default;
  Product(const Product &) = default;
  Product(Product &&) = default;
  Product &operator=(const Product &) = default;
  Product &operator=(Product &&) = default;

  /// construct a Product out of zero or more factors (multiplied by 1)
  /// @param factors the factors
  Product(ExprPtrList factors) {
    using std::begin;
    using std::end;
    for (auto it = begin(factors); it != end(factors); ++it) append(1, *it);
  }

  /// construct a Product out of zero or more factors multiplied by a scalar
  /// @tparam T a numeric type; it must be able to multiply std::complex<double>
  /// @param scalar a scalar of type T
  /// @param factors an initializer list of factors
  template <typename T>
  Product(T scalar, ExprPtrList factors) : scalar_(std::move(scalar)) {
    using std::begin;
    using std::end;
    for (auto it = begin(factors); it != end(factors); ++it) append(1, *it);
  }

  /// construct a Product out of a range of factors
  /// @param begin the begin iterator
  /// @param end the end iterator
  template <typename Iterator>
  Product(Iterator begin, Iterator end) {
    for (auto it = begin; it != end; ++it) append(1, *it);
  }

  /// construct a Product out of a range of factors
  /// @tparam T a numeric type; it must be able to multiply std::complex<double>
  /// @param scalar a scalar of type T
  /// @param begin the begin iterator
  /// @param end the end iterator
  template <typename T, typename Iterator>
  Product(T scalar, Iterator begin, Iterator end) : scalar_(std::move(scalar)) {
    for (auto it = begin; it != end; ++it) append(1, *it);
  }

  /// multiplies the product by @c scalar
  template <typename T>
  Product &scale(T scalar) {
    scalar_ *= scalar;
    return *this;
  }

  /// (post-)multiplies the product by @c scalar times @c factor
  template <typename T>
  Product &append(T scalar, ExprPtr factor) {
    assert(factor);
    scalar_ *= scalar;
    if (!factor->is<Product>()) {
      if (factor->is<Constant>()) {  // factor in Constant
        auto factor_constant = factor->as<Constant>();
        scalar_ *= factor_constant.value();
        // no need to reset the hash since scalar is not hashed!
      } else {
        factors_.push_back(std::move(factor));
        reset_hash_value();
      }
    } else {  // factor is a product also ... flatten recursively
      auto factor_product = factor->as<Product>();
      scalar_ *= factor_product.scalar_;
      for (auto &&subfactor : factor_product) this->append(1, subfactor);
      //      using std::end;
      //      using std::cbegin;
      //      using std::cend;
      //      factors_.insert(end(factors_), cbegin(factor_product->factors_),
      //      cend(factor_product->factors_));
    }
    return *this;
  }

  /// append @c factor without flattening
  Product &append(ExprPtr factor) {
    if (factor->is<Constant>()) {
      auto factor_constant = factor->as<Constant>();
      scalar_ *= factor_constant.value();
    } else {
    factors_.push_back(std::move(factor));
    reset_hash_value();
    }

    return *this;
  }

  /// (post-)multiplies the product by @c scalar times @c factor
  template <typename T, typename Factor, typename = std::enable_if_t<std::is_base_of_v<Expr, std::remove_reference_t<Factor>>>>
  Product &append(T scalar, Factor&& factor) {
    return this->append(scalar,
                        std::static_pointer_cast<Expr>(
                            std::forward<Factor>(factor).shared_from_this()));
  }

  /// (pre-)multiplies the product by @c scalar times @c factor ; less efficient
  /// than append()
  template <typename T>
  Product &prepend(T scalar, ExprPtr factor) {
    assert(factor);
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
      for (auto &subfactor : *factor_product) this->prepend(1, subfactor);
      //      using std::begin;
      //      using std::cbegin;
      //      using std::cend;
      //      factors_.insert(begin(factors_), cbegin(factor_product->factors_),
      //      cend(factor_product->factors_));
    }
    return *this;
  }

  /// (pre-)multiplies the product by @c scalar times @c factor
  template <typename T, typename Factor, typename = std::enable_if_t<std::is_base_of_v<Expr, std::remove_reference_t<Factor>>>>
  Product &prepend(T scalar, Factor&& factor) {
    return this->prepend(scalar,
                         std::static_pointer_cast<Expr>(
                             std::forward<Factor>(factor).shared_from_this()));
  }

  const std::complex<double> &scalar() const { return scalar_; }
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

 private:
  /// @return true if commutativity is decidable statically
  /// @sa CProduct::static_commutativity() and NCProduct::static_commutativity()
  virtual bool static_commutativity() const { return false; }

 public:
  std::wstring to_latex() const override {
    return to_latex(false);
  }

  /// just like Expr::to_latex() , but can negate before conversion
  /// @param[in] negate if true, scalar will be before conversion
  std::wstring to_latex(bool negate) const {
    std::wstring result;
    result = L"{";
    if (scalar() != 0.) {
      const auto scal = negate ? -scalar() : scalar();
      if (scal != 1.) {
        result += sequant::to_latex(scal);
      }
      for (const auto &i : factors()) {
        if (i->is<Product>())
          result += L"\\left(" + i->to_latex() + L"\\right)";
        else result += i->to_latex();
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

  std::shared_ptr<Expr> clone() const override {
    auto cloned_factors =
        factors() | ranges::views::transform([](const ExprPtr &ptr) {
          return ptr ? ptr->clone() : nullptr;
        });
    return ex<Product>(this->scalar(), ranges::begin(cloned_factors),
                       ranges::end(cloned_factors));
  }

  Product deep_copy() const {
    auto cloned_factors =
        factors() | ranges::views::transform([](const ExprPtr &ptr) {
          return ptr ? ptr->clone() : nullptr;
        });
    return Product(this->scalar(), ranges::begin(cloned_factors),
                   ranges::end(cloned_factors));
  }

  virtual Expr &operator*=(const Expr &that) override {
    if (!that.is<Constant>()) {
      this->append(1, const_cast<Expr &>(that).shared_from_this());
    }
    else {
      scalar_ *= that.as<Constant>().value();
    }
    return *this;
  }

  void add_identical(const std::shared_ptr<Product> &other) {
    assert(this->hash_value() == other->hash_value());
    scalar_ += other->scalar_;
  }

 private:
  std::complex<double> scalar_ = {1.0, 0.0};
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

  /// @note this hashes only the factors, not the scalar to make possible rapid
  /// finding of identical factors
  hash_type memoizing_hash() const override {
    auto deref_factors =
        factors() |
        ranges::views::transform(
            [](const ExprPtr &ptr) -> const Expr & { return *ptr; });
    hash_value_ = boost::hash_range(ranges::begin(deref_factors),
                                    ranges::end(deref_factors));
    return *hash_value_;
  }

  std::shared_ptr<Expr> canonicalize_impl(bool rapid = false);
  virtual std::shared_ptr<Expr> canonicalize() override;
  virtual std::shared_ptr<Expr> rapid_canonicalize() override;

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

using ProductPtr = std::shared_ptr<Product>;

class CProduct : public Product {
 public:
  using Product::Product;
  CProduct(const Product &other) : Product(other) {}
  CProduct(Product &&other) : Product(other) {}

  bool is_commutative() const override { return true; }

 private:
  bool static_commutativity() const override { return true; }
};  // class CProduct

using CProductPtr = std::shared_ptr<CProduct>;

class NCProduct : public Product {
 public:
  using Product::Product;
  NCProduct(const Product &other) : Product(other) {}
  NCProduct(Product &&other) : Product(other) {}

  bool is_commutative() const override { return false; }

 private:
  bool static_commutativity() const override { return true; }
};  // class NCProduct

using NCProductPtr = std::shared_ptr<NCProduct>;

/// @brief sum of zero or more summands

/// Sum is associative and is flattened automatically.
class Sum : public Expr {
 public:
  Sum() = default;
  virtual ~Sum() = default;
  Sum(const Sum &) = default;
  Sum(Sum &&) = default;
  Sum &operator=(const Sum &) = default;
  Sum &operator=(Sum &&) = default;

  /// construct a Sum out of zero or more summands
  /// @param summands an initializer list of summands
  Sum(ExprPtrList summands) {
    // use append to flatten out Sum summands
    for (auto &summand : summands) {
      append(std::move(summand));
    }
  }

  /// construct a Sum out of a range of summands
  /// @param begin the begin iterator
  /// @param end the end iterator
  template <typename Iterator>
  Sum(Iterator begin, Iterator end) {
    // use append to flatten out Sum summands
    for (auto it = begin; it != end; ++it) {
      append(*it);
    }
  }

  /// append a summand to the sum
  /// @param summand the summand
  Sum &append(ExprPtr summand) {
    assert(summand);
    if (!summand->is<Sum>()) {
      if (summand->is<Constant>()) {  // exclude zeros, add up constants
                                      // immediately, if possible
        auto summand_constant = std::static_pointer_cast<Constant>(summand);
        if (constant_summand_idx_) {
          assert(summands_.at(*constant_summand_idx_)->is<Constant>());
          *(summands_[*constant_summand_idx_]) += *summand;
        } else {
          if (*summand_constant != Constant(0)) {
            summands_.push_back(std::move(summand));
            constant_summand_idx_ = summands_.size() - 1;
          }
        }
      } else {
        summands_.push_back(std::move(summand));
      }
      reset_hash_value();
    } else {  // this recursively flattens Sum summands
      for (auto &subsummand : *summand) this->append(subsummand);
    }
    return *this;
  }

  /// prepend a summand to the sum
  /// @param summand the summand
  Sum &prepend(ExprPtr summand) {
    assert(summand);
    if (!summand->is<Sum>()) {
      if (summand->is<Constant>()) {  // exclude zeros
        auto summand_constant = std::static_pointer_cast<Constant>(summand);
        if (constant_summand_idx_) {  // add up to the existing constant ...
          assert(summands_.at(*constant_summand_idx_)->is<Constant>());
          *summands_[*constant_summand_idx_] += *summand_constant;
        } else {  // or include the nonzero constant and update
                  // constant_summand_idx_
          if (summand_constant != 0) {
            summands_.insert(summands_.begin(), std::move(summand));
            constant_summand_idx_ = 0;
          }
        }
      } else {
        summands_.insert(summands_.begin(), std::move(summand));
        if (constant_summand_idx_)  // if have a constant, update its position
          ++*constant_summand_idx_;
      }
      reset_hash_value();
    } else {  // this recursively flattens Sum summands
      for (auto &subsummand : *summand) this->prepend(subsummand);
    }
    return *this;
  }

  /// Summands accessor
  const auto &summands() const { return summands_; }

  /// Summand accessor
  /// @param i summand index
  /// @return ith summand
  const ExprPtr &summand(size_t i) const { return summands_.at(i); }

  /// Takes the first @c count elements of the sum
  ExprPtr take_n(size_t count) const {
    const auto e = (count >= summands_.size()? summands_.end() : (summands_.begin() + count));
    return ex<Sum>(summands_.begin(), e);
  }

  /// Takes the first @c count elements of the sum starting with element @c offset
  ExprPtr take_n(size_t offset, size_t count) const {
    const auto offset_plus_count = offset + count;
    const auto b = (offset >= summands_.size() ? summands_.end() : (summands_.begin() + offset));
    const auto e = (offset_plus_count >= summands_.size() ? summands_.end() : (summands_.begin() + offset_plus_count));
    return ex<Sum>(b, e);
  }

  /// @return true if the number of factors is zero
  bool empty() const { return summands_.empty(); }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{ \\left(";
    std::size_t counter = 0;
    for (const auto &i : summands()) {
      const auto i_is_product = i->is<Product>();
      if (!i_is_product) {
        result += (counter == 0) ? i->to_latex() : (L" + " + i->to_latex());
      } else {  // i_is_product
        const auto i_prod = i->as<Product>();
        const auto scalar = i_prod.scalar();
        if (scalar.real() < 0 || (scalar.real() == 0 && scalar.imag() < 0)) {
          result += L" - " + i_prod.to_latex(true);
        } else {
          result += (counter == 0) ? i->to_latex() : (L" + " + i->to_latex());
        }
      }
      ++counter;
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
      if (counter != summands().size()) result += L",";
    }
    result += L"]";
    return result;
  }

  Expr::type_id_type type_id() const override {
    return Expr::get_type_id<Sum>();
  };

  std::shared_ptr<Expr> clone() const override {
    auto cloned_summands =
        summands() | ranges::views::transform(
                         [](const ExprPtr &ptr) { return ptr->clone(); });
    return ex<Sum>(ranges::begin(cloned_summands),
                   ranges::end(cloned_summands));
  }

  virtual Expr &operator+=(const Expr &that) override {
    this->append(const_cast<Expr &>(that).shared_from_this());
    return *this;
  }

  virtual Expr &operator-=(const Expr &that) override {
    if (that.is<Constant>())
      this->append(ex<Constant>(-that.as<Constant>().value()));
    else
      this->append(ex<Product>(
          -1, ExprPtrList{const_cast<Expr &>(that).shared_from_this()}));
    return *this;
  }

 private:
  container::svector<ExprPtr, 2> summands_{};
  std::optional<size_t>
      constant_summand_idx_{};  // points to the constant summand, if any; used
                                // to sum up constants in append/prepend

  cursor begin_cursor() override {
    return summands_.empty() ? Expr::begin_cursor() : cursor{&summands_[0]};
  };
  cursor end_cursor() override {
    return summands_.empty() ? Expr::end_cursor()
                             : cursor{&summands_[0] + summands_.size()};
  };
  cursor begin_cursor() const override {
    return summands_.empty() ? Expr::begin_cursor() : cursor{&summands_[0]};
  };
  cursor end_cursor() const override {
    return summands_.empty() ? Expr::end_cursor()
                             : cursor{&summands_[0] + summands_.size()};
  };

  hash_type memoizing_hash() const override {
    auto deref_summands =
        summands() |
        ranges::views::transform(
            [](const ExprPtr &ptr) -> const Expr & { return *ptr; });
    hash_value_ = boost::hash_range(ranges::begin(deref_summands),
                                    ranges::end(deref_summands));
    return *hash_value_;
  }

  /// @param multipass if true, will do a multipass canonicalization, with extra
  /// cleanup pass after the deep canonization pass
  ExprPtr canonicalize_impl(bool multipass);

  virtual ExprPtr canonicalize() override { return canonicalize_impl(true); }
  virtual ExprPtr rapid_canonicalize() override {
    return canonicalize_impl(false);
  }

  bool static_equal(const Expr &that) const override {
    const auto &that_cast = static_cast<const Sum &>(that);
    if (summands().size() == that_cast.summands().size()) {
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
};  // class Sum

using SumPtr = std::shared_ptr<Sum>;

inline std::wstring to_latex(const ExprPtr &exprptr) {
  return exprptr->to_latex();
}

/// splits long outer sum into a multiline align
/// @param exprptr the expression to be converted to a string
/// @param max_lines_per_align the maximum number of lines in the align before
/// starting new align block (if zero, will produce single align block)
/// @param max_terms_per_line the maximum number of terms per line
inline std::wstring to_latex_align(const ExprPtr &exprptr,
                                   size_t max_lines_per_align = 0,
                                   size_t max_terms_per_line = 1) {
  std::wstring result = to_latex(exprptr);
  if (exprptr->is<Sum>()) {
    result.erase(0, 7);  // remove leading  "{ \left"
    result.replace(result.size() - 9, 9,
                   L")");  // replace trailing "\right) }" with ")"
    result = std::wstring(L"\\begin{align}\n& ") + result;
    // assume no inner sums
    int line_counter = 0;
    int term_counter = 0;
    std::wstring::size_type pos = 0;
    std::wstring::size_type plus_pos = 0;
    std::wstring::size_type minus_pos = 0;
    bool last_pos_has_plus = false;
    bool have_next_term = true;
    auto insert_into_result_at = [&](std::wstring::size_type at, const auto& str) {
      assert(pos != std::wstring::npos);
      result.insert(at, str);
      const auto str_nchar = std::size(str) - 1;  // neglect end-of-string
      pos += str_nchar;
      if (plus_pos != std::wstring::npos)
        plus_pos += str_nchar;
      if (minus_pos != std::wstring::npos)
        minus_pos += str_nchar;
      if (pos != plus_pos)
        assert(plus_pos == result.find(L" + ", plus_pos));
      if (pos != minus_pos)
        assert(minus_pos == result.find(L" - ", minus_pos));
    };
    while (have_next_term) {
      if (max_lines_per_align > 0 &&
          line_counter == max_lines_per_align) {  // start new align block?
        insert_into_result_at(pos + 1, L"\n\\end{align}\n\\begin{align}\n& ");
        line_counter = 0;
      } else {
        // break the line if needed
        if (term_counter != 0 && term_counter % max_terms_per_line == 0) {
          insert_into_result_at(pos + 1, L"\\\\\n& ");
          ++line_counter;
        }
      }
      // next term, plz
      if (plus_pos == 0 || last_pos_has_plus)
        plus_pos = result.find(L" + ", plus_pos + 1);
      if (minus_pos == 0 || !last_pos_has_plus)
        minus_pos = result.find(L" - ", minus_pos + 1);
      pos = std::min(plus_pos, minus_pos);
      last_pos_has_plus = (pos == plus_pos);
      if (pos != std::wstring::npos)
        ++term_counter;
      else
        have_next_term = false;
    }
  } else {
    result = std::wstring(L"\\begin{align}\n& ") + result;
  }
  result += L"\n\\end{align}";
  return result;
}

inline std::wstring to_wolfram(const ExprPtr &exprptr) {
  return exprptr->to_wolfram();
}

template<typename Sequence>
std::decay_t<Sequence> clone(Sequence &&exprseq) {
  auto cloned_seq = exprseq | ranges::views::transform([](const ExprPtr &ptr) { return ptr ? ptr->clone() : nullptr; });
  return std::decay_t<Sequence>(ranges::begin(cloned_seq), ranges::end(cloned_seq));
}

};  // namespace sequant

#include "expr_operator.hpp"

#include "expr_algorithm.hpp"

#endif  // SEQUANT_EXPR_HPP
