//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT_EXPR_HPP
#define SEQUANT_EXPR_HPP

#include "SeQuant/core/expr_fwd.hpp"

#include "SeQuant/core/complex.hpp"
#include "SeQuant/core/container.hpp"
#include "SeQuant/core/hash.hpp"
#include "SeQuant/core/latex.hpp"
#include "SeQuant/core/logger.hpp"
#include "SeQuant/core/meta.hpp"
#include "SeQuant/core/rational.hpp"
#include "SeQuant/core/wolfram.hpp"

#include <range/v3/all.hpp>

#include <boost/core/demangle.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <atomic>
#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

namespace sequant {

/// @brief ExprPtr is a multiple-owner smart pointer to Expr

/// It can be used mostly interchangeably with `std::shared_ptr<Expr>`, but
/// also provides convenient mathematical operators (`+=`, etc.)
class ExprPtr : public std::shared_ptr<Expr> {
 public:
  using base_type = std::shared_ptr<Expr>;
  using base_type::operator->;
  using base_type::base_type;

  ExprPtr() = default;
  template <typename E, typename = std::enable_if_t<
                            std::is_same_v<std::remove_const_t<E>, Expr> ||
                            std::is_base_of_v<Expr, std::remove_const_t<E>>>>
  ExprPtr(const std::shared_ptr<E> &other_sptr) : base_type(other_sptr) {}
  template <typename E, typename = std::enable_if_t<
                            std::is_same_v<std::remove_const_t<E>, Expr> ||
                            std::is_base_of_v<Expr, std::remove_const_t<E>>>>
  ExprPtr(std::shared_ptr<E> &&other_sptr) : base_type(std::move(other_sptr)) {}
  template <typename E, typename = std::enable_if_t<
                            std::is_same_v<std::remove_const_t<E>, Expr> ||
                            std::is_base_of_v<Expr, std::remove_const_t<E>>>>
  ExprPtr &operator=(const std::shared_ptr<E> &other_sptr) {
    as_shared_ptr() = other_sptr;
    return *this;
  }
  template <typename E, typename = std::enable_if_t<
                            std::is_same_v<std::remove_const_t<E>, Expr> ||
                            std::is_base_of_v<Expr, std::remove_const_t<E>>>>
  ExprPtr &operator=(std::shared_ptr<E> &&other_sptr) {
    as_shared_ptr() = std::move(other_sptr);
    return *this;
  }

  ~ExprPtr() = default;

  base_type &as_shared_ptr() &;
  const base_type &as_shared_ptr() const &;
  base_type &&as_shared_ptr() &&;

  template <typename E, typename = std::enable_if_t<!std::is_same_v<E, Expr>>>
  std::shared_ptr<E> as_shared_ptr() const {
    assert(this->is<E>());
    return std::static_pointer_cast<E>(this->as_shared_ptr());
  }

  ExprPtr &operator+=(const ExprPtr &);
  ExprPtr &operator-=(const ExprPtr &);
  ExprPtr &operator*=(const ExprPtr &);

  /// @tparam T an Expr type
  /// @return true if this object is of type @c T
  template <typename T>
  bool is() const;

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template <typename T>
  const T &as() const;

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template <typename T>
  T &as();

  std::wstring to_latex() const;
};  // class ExprPtr

/// ExprPtr is equal to a null pointer if it's uninitialized
inline bool operator==(const ExprPtr &x, std::nullptr_t) {
  return x.get() == nullptr;
}

/// ExprPtr is equal to a null pointer if it's uninitialized
inline bool operator==(std::nullptr_t, const ExprPtr &x) {
  return x.get() == nullptr;
}

/// @brief Base expression class

/// Expr represents the interface needed to form expression trees. Classes that
/// represent expressions should publicly derive from this class. Each Expr on a
/// tree has links to its children Expr objects. The lifetime of Expr objects is
/// expected to be managed by std::shared_ptr . Expr is an Iterable over
/// subexpressions (each of which is an Expr itself). More precisely, Expr meets
/// the SizedIterable concept (see
/// https://raw.githubusercontent.com/ericniebler/range-v3/master/doc/std/D4128.md).
/// Specifically, iterators to subexpressions
/// dereference to ExprPtr. Since Expr is a range, it provides begin/end/etc.
/// that can participate in overloads
///       with other functions in the derived class. Consider a Container class
///       derived from a BaseContainer class:
/// @code
///   template <typename T> class Container : public BaseContainer, public Expr
///   {
///     // WARNING: BaseContainer::begin clashes with Expr::begin
///     // WARNING: BaseContainer::end clashes with Expr::end
///     // etc.
///   };
/// @endcode
/// There are two possible scenarios:
///   - if @c Container is a container of Expr objects, BaseContainer will
///   iterate over ExprPtr objects already
///     and both ranges will be equivalent; it is sufficient to add `using
///     BaseContainer::begin`, etc. to Container's public API.
///   - if @c Container is a container of non-Expr objects, iteration over
///   BaseContainer is likely to be more commonly used
///     in practice, hence again adding `using BaseContainer::begin`, etc. will
///     suffice. To be able to iterate over subexpression range (in this case it
///     is empty) Expr provides Expr::expr member to cast to Expr:
/// @code
///    Container c(...);
///    for(const auto& e: c) {  // iterates over elements of BaseContainer
///    }
///    for(const auto& e: c.expr()) {  // iterates over subexpressions
///    }
/// @endcode
class Expr : public std::enable_shared_from_this<Expr>,
             public ranges::view_facade<Expr> {
 public:
  using range_type = ranges::view_facade<Expr>;
  using hash_type = std::size_t;
  using type_id_type = int;  // to speed up comparisons

  Expr() = default;
  virtual ~Expr() = default;

  /// @return true if this is a leaf
  bool is_atom() const { return ranges::empty(*this); }

  /// @return the string representation of @c this in the LaTeX format
  virtual std::wstring to_latex() const;

  /// @return the string representation of @c this in the Wolfram Language
  /// format
  virtual std::wstring to_wolfram() const;

  /// @return a clone of this object, i.e. an object that is equal to @c this
  /// @note - must be overridden in the derived class.
  ///       - the default implementation throws an exception
  virtual ExprPtr clone() const;

  /// Canonicalizes @c this and returns the biproduct of canonicalization (e.g.
  /// phase)
  /// @return the biproduct of canonicalization, or @c nullptr if no biproduct
  /// generated
  virtual ExprPtr canonicalize() {
    return {};  // by default do nothing and return nullptr
  }

  /// Performs approximate, but fast, canonicalization of @c this and returns
  /// the biproduct of canonicalization (e.g. phase) The default is to use
  /// canonicalize(), unless overridden in the derived class.
  /// @return the biproduct of canonicalization, or @c nullptr if no biproduct
  /// generated
  virtual ExprPtr rapid_canonicalize() { return this->canonicalize(); }

  // clang-format off
  /// recursively visit this expression, i.e. call visitor on each subexpression
  /// in depth-first fashion.
  /// @warning this will only work for tree expressions; no checking is
  /// performed that each subexpression has only been visited once
  /// TODO make work for graphs
  /// @tparam Visitor a callable of type void(ExprPtr&) or void(const ExprPtr&)
  /// @param visitor the visitor object
  /// @param atoms_only if true, will visit only the leaves; the default is to
  /// visit all nodes
  /// @return true if this object was visited
  /// @sa expr_range
  // clang-format on
  template <typename Visitor>
  bool visit(Visitor &&visitor, const bool atoms_only = false) {
    for (auto &subexpr_ptr : expr()) {
      const auto subexpr_is_an_atom = subexpr_ptr->is_atom();
      const auto need_to_visit_subexpr = !atoms_only || subexpr_is_an_atom;
      bool visited = false;
      if (!subexpr_is_an_atom)  // if not a leaf, recur into it
        visited =
            subexpr_ptr->visit(std::forward<Visitor>(visitor), atoms_only);
      // call on the subexpression itself, if not yet done so
      if (need_to_visit_subexpr && !visited) visitor(subexpr_ptr);
    }
    // can only visit itself here if visitor(const ExprPtr&) is valid
    bool this_visited = false;
    if constexpr (std::is_invocable_r_v<void, std::remove_reference_t<Visitor>,
                                        const ExprPtr &>) {
      if (!atoms_only || this->is_atom()) {
        visitor(shared_from_this());
        this_visited = true;
      }
    }
    return this_visited;
  }

  auto begin_subexpr() { return range_type::begin(); }

  auto end_subexpr() { return range_type::end(); }

  auto begin_subexpr() const { return range_type::begin(); }

  auto end_subexpr() const { return range_type::end(); }

  Expr &expr() { return *this; }

  template <typename T, typename Enabler = void>
  struct is_shared_ptr_of_expr : std::false_type {};
  template <typename T>
  struct is_shared_ptr_of_expr<std::shared_ptr<T>,
                               std::enable_if_t<std::is_same_v<Expr, T>>>
      : std::true_type {};
  template <typename T, typename Enabler = void>
  struct is_shared_ptr_of_expr_or_derived : std::false_type {};
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

  /// @brief changes this to its adjoint
  /// @note base implementation throws, must be reimplemented in the derived
  /// class
  virtual void adjoint();

  /// Computes and returns the hash value. If default @p hasher is used then the
  /// value will be memoized, otherwise @p hasher will be used to compute the
  /// hash every time.
  /// @param hasher the hasher object, if omitted, the default is used (@sa
  /// Expr::memoizing_hash )
  /// @note always returns 0 unless this derived class overrides
  /// Expr::memoizing_hash
  /// @return the hash value for this Expr
  hash_type hash_value(
      std::function<hash_type(const std::shared_ptr<const Expr> &)> hasher = {})
      const {
    return hasher ? hasher(shared_from_this()) : memoizing_hash();
  }

  /// Computes and returns the derived type identifier
  /// @note this function must be overridden in the derived class
  /// @sa Expr::get_type_id
  /// @return the hash value for this Expr
  virtual type_id_type type_id() const
#if __GNUG__
  {
    abort();
  }
#else
      = 0;
#endif

  friend inline bool operator==(const Expr &a, const Expr &b);

  /// @tparam T Expr or a class derived from Expr
  /// @return true if @c *this is less than @c that
  /// @note the derived class must implement Expr::static_less_than
  template <typename T>
  std::enable_if_t<std::is_base_of<Expr, T>::value, bool> operator<(
      const T &that) const {
    if (type_id() ==
        that.type_id()) {  // if same type, use generic (or type-specific, if
                           // available) comparison
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
    if constexpr (std::is_same_v<std::decay_t<T>, Expr>)
      return true;
    else
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
  virtual Expr &operator*=(const Expr &that);

  /// @brief in-place non-commutatively-multiply @c *this by @c that
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator^=(const Expr &that);

  /// @brief in-place add @c that to @c *this
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator+=(const Expr &that);

  /// @brief in-place subtract @c that from @c *this
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator-=(const Expr &that);

  ///@}

 private:
  friend ranges::range_access;

 protected:
  Expr(Expr &&) = default;
  Expr(const Expr &) = default;
  Expr &operator=(Expr &&) = default;
  Expr &operator=(const Expr &) = default;

  struct cursor {
    using value_type = ExprPtr;

    cursor() = default;
    constexpr explicit cursor(ExprPtr *subexpr_ptr) noexcept
        : ptr_{subexpr_ptr} {}
    /// when take const ptr note runtime const flag
    constexpr explicit cursor(const ExprPtr *subexpr_ptr) noexcept
        : ptr_{const_cast<ExprPtr *>(subexpr_ptr)}, const_{true} {}
    bool equal(const cursor &that) const { return ptr_ == that.ptr_; }
    void next() { ++ptr_; }
    void prev() { --ptr_; }
    // TODO figure out why can't return const here if want to be able to assign
    // to *begin(Expr&)
    ExprPtr &read() const {
      RANGES_EXPECT(ptr_);
      return *ptr_;
    }
    ExprPtr &read() {
      RANGES_EXPECT(const_ == false);
      RANGES_EXPECT(ptr_);
      return *ptr_;
    }
    void assign(const ExprPtr &that_ptr) {
      RANGES_EXPECT(ptr_);
      *ptr_ = that_ptr;
    }
    std::ptrdiff_t distance_to(cursor const &that) const {
      return that.ptr_ - ptr_;
    }
    void advance(std::ptrdiff_t n) { ptr_ += n; }

   private:
    ExprPtr *ptr_ =
        nullptr;  // both begin and end will be represented by this, so Expr
                  // without subexpressions begin() equals end() automatically
    bool const_ = false;  // assert in nonconst ops
  };

  /// @return the cursor for the beginning of the range (must override in a
  /// derived Expr that has subexpressions)
  virtual cursor begin_cursor() { return cursor{}; }
  /// @return the cursor for the end of the range (must override in a derived
  /// Expr that has subexpressions)
  virtual cursor end_cursor() { return cursor{}; }

  /// @return the cursor for the beginning of the range (must override in a
  /// derived Expr that has subexpressions)
  virtual cursor begin_cursor() const { return cursor{}; }
  /// @return the cursor for the end of the range (must override in a derived
  /// Expr that has subexpressions)
  virtual cursor end_cursor() const { return cursor{}; }

  mutable std::optional<hash_type> hash_value_;  // not initialized by default
  virtual hash_type memoizing_hash() const {
    static const hash_type default_hash_value = 0;
    if (hash_value_)
      return *hash_value_;
    else
      return default_hash_value;
  }
  virtual void reset_hash_value() const { hash_value_.reset(); }

  /// @param that an Expr object
  /// @note @c that is guaranteed to be of same type as @c *this, hence can be
  /// statically cast
  /// @return true if @c that is equivalent to *this
  virtual bool static_equal(const Expr &that) const
#if __GNUG__
  {
    abort();
  }
#else
      = 0;
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
  virtual bool commutes_with_atom(const Expr &that) const { return true; }

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

 private:
  /// @input[in] fn the name of function that is missing in this class
  /// @return a std::logic_error object containing a message describing that @p
  /// fn is missing from this type
  std::logic_error not_implemented(const char *fn) const;
};  // class Expr

template <>
struct Expr::is_shared_ptr_of_expr<ExprPtr, void> : std::true_type {};
template <>
struct Expr::is_shared_ptr_of_expr_or_derived<ExprPtr, void> : std::true_type {
};

/// @return true if @c a is equal to @c b
inline bool operator==(const Expr &a, const Expr &b) {
  if (a.type_id() != b.type_id())
    return false;
  else
    return a.static_equal(b);
}

#if __cplusplus < 202002L
/// @return true if @c a is not equal to @c b
inline bool operator!=(const Expr &a, const Expr &b) { return !(a == b); }
#endif  // __cplusplus < 202002L

/// make an ExprPtr to a new object of type T
/// @tparam T a class derived from Expr
/// @tparam Args a parameter pack type such that T(std::forward<Args>...) is
/// well-formed
/// @param args a parameter pack such that T(args...) is well-formed
template <typename T, typename... Args>
ExprPtr ex(Args &&...args) {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

// this is needed when using std::make_shared<X>({ExprPtr,ExprPtr}), i.e. must
// std::make_shared<X>(ExprPtrList{ExprPtr,ExprPtr})
using ExprPtrList = std::initializer_list<ExprPtr>;
static auto expr_ptr_comparer = [](const auto &ptr1, const auto &ptr2) {
  return *ptr1 == *ptr2;
};

using ExprPtrVector = container::svector<ExprPtr>;

/// @brief computes the adjoint of @p expr
/// @param[in] expr an Expr object
/// @return the adjoint of @p expr
ExprPtr adjoint(const ExprPtr &expr);

/// An object with a string label that be used for defining a canonical order of
/// expressions (defined at runtime)
class Labeled {
 public:
  Labeled() = default;
  virtual ~Labeled() = default;

  virtual std::wstring_view label() const = 0;
};

/// @brief a constant number

/// This is represented as a "compile-time" complex rational number
class Constant : public Expr {
 public:
  using scalar_type = Complex<sequant::rational>;

 private:
  scalar_type value_;

 public:
  Constant() = delete;
  virtual ~Constant() = default;
  Constant(const Constant &) = default;
  Constant(Constant &&) = default;
  Constant &operator=(const Constant &) = default;
  Constant &operator=(Constant &&) = default;
  template <typename U, typename = std::enable_if_t<
                            !std::is_same_v<std::decay_t<U>, Constant>>>
  explicit Constant(U &&value) : value_(std::forward<U>(value)) {}

 private:
  template <typename X>
  static X numeric_cast(const sequant::rational &r) {
    if constexpr (std::is_integral_v<X>) {
      assert(denominator(r) == 1);
      return boost::numeric_cast<X>(numerator(r));
    } else {
      return boost::numeric_cast<X>(numerator(r)) /
             boost::numeric_cast<X>(denominator(r));
    }
  };

 public:
  /// @tparam T the result type; default to the type of value_
  /// @return the value cast to ResultType
  /// @throw std::invalid_argument if conversion to T is not possible
  /// @throw boost::numeric::positive_overflow or
  /// boost::numeric::negative_overflow if cast fails
  template <typename T = decltype(value_)>
  auto value() const {
    if constexpr (std::is_arithmetic_v<T>) {
      assert(value_.imag() == 0);
      return numeric_cast<T>(value_.real());
    } else if constexpr (meta::is_complex_v<T>) {
      return T(numeric_cast<typename T::value_type>(value_.real()),
               numeric_cast<typename T::value_type>(value_.imag()));
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

  type_id_type type_id() const override { return get_type_id<Constant>(); }

  ExprPtr clone() const override { return ex<Constant>(this->value()); }

  /// @brief adjoint of a Constant is its complex conjugate
  virtual void adjoint() override;

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

  /// @param[in] v a scalar
  /// @return true if this is a soft zero, i.e. its magnitude is less than
  /// `std::sqrt(std::numeric_limits<float>::epsilon())`
  static bool is_zero(scalar_type v) { return v.is_zero(); }

  /// @return `Constant::is_zero(this->value())`
  bool is_zero() const { return is_zero(this->value()); }

 private:
  hash_type memoizing_hash() const override {
    hash_value_ = hash::value(value_);
    return *hash_value_;
  }

  bool static_equal(const Expr &that) const override {
    return value() == static_cast<const Constant &>(that).value();
  }
};  // class Constant

/// This is represented as a "run-time" complex rational number
class Variable : public Expr, public Labeled {
 public:
  Variable() = delete;
  virtual ~Variable() = default;
  Variable(const Variable &) = default;
  Variable(Variable &&) = default;
  Variable &operator=(const Variable &) = default;
  Variable &operator=(Variable &&) = default;
  template <typename U, typename = std::enable_if_t<
                            !std::is_same_v<std::decay_t<U>, Variable>>>
  explicit Variable(U &&label) : label_(std::forward<U>(label)) {}

  Variable(std::wstring label, bool conjugated)
      : label_(std::move(label)), conjugated_(conjugated) {}


  std::wstring_view label() const override;

  bool conjugated() const;

  std::wstring to_latex() const override;

  type_id_type type_id() const override { return get_type_id<Variable>(); }

  ExprPtr clone() const override;

  /// @brief adjoint of a Variable is its complex conjugate
  virtual void adjoint() override;

 private:
  std::wstring label_;
  bool conjugated_ = false;

  hash_type memoizing_hash() const override {
    hash_value_ = hash::value(label_);
    hash::combine(hash_value_.value(), conjugated_);
    return *hash_value_;
  }

  bool static_equal(const Expr &that) const override {
    return label_ == static_cast<const Variable &>(that).label_ &&
           conjugated_ == static_cast<const Variable &>(that).conjugated_;
  }
};  // class Variable

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
            typename = std::enable_if_t<
                meta::is_range_v<std::decay_t<Range>> &&
                !std::is_same_v<std::remove_reference_t<Range>, ExprPtrList> &&
                !std::is_same_v<std::remove_reference_t<Range>, Product>>>
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
            typename = std::enable_if_t<
                std::is_base_of_v<Expr, std::remove_reference_t<Factor>>>>
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
            typename = std::enable_if_t<
                std::is_base_of_v<Expr, std::remove_reference_t<Factor>>>>
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
        result += sequant::to_latex(scal);
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
      if (cloned_factor.template is<Product>()) std::cout << "";
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

  /// @note this hashes only the factors, not the scalar to make possible rapid
  /// finding of identical factors
  hash_type memoizing_hash() const override {
    auto deref_factors =
        factors() |
        ranges::views::transform(
            [](const ExprPtr &ptr) -> const Expr & { return *ptr; });
    hash_value_ =
        hash::range(ranges::begin(deref_factors), ranges::end(deref_factors));
    return *hash_value_;
  }

  ExprPtr canonicalize_impl(bool rapid = false);
  virtual ExprPtr canonicalize() override;
  virtual ExprPtr rapid_canonicalize() override;

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

  /// construct a Sum out of a range of summands
  /// @param rng a range
  template <typename Range,
            typename = std::enable_if_t<
                meta::is_range_v<std::decay_t<Range>> &&
                !std::is_same_v<std::remove_reference_t<Range>, ExprPtrList>>>
  explicit Sum(Range &&rng) {
    // use append to flatten out Sum summands
    for (auto &&v : rng) {
      append(std::forward<decltype(v)>(v));
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
          if (!summand_constant->is_zero()) {
            summands_.push_back(summand->clone());
            constant_summand_idx_ = summands_.size() - 1;
          }
        }
      } else {
        summands_.push_back(summand->clone());
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
          if (!summand_constant->is_zero()) {
            summands_.insert(summands_.begin(), summand->clone());
            constant_summand_idx_ = 0;
          }
        }
      } else {
        summands_.insert(summands_.begin(), summand->clone());
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
    const auto e = (count >= summands_.size() ? summands_.end()
                                              : (summands_.begin() + count));
    return ex<Sum>(summands_.begin(), e);
  }

  /// Takes the first @c count elements of the sum starting with element @c
  /// offset
  ExprPtr take_n(size_t offset, size_t count) const {
    const auto offset_plus_count = offset + count;
    const auto b = (offset >= summands_.size() ? summands_.end()
                                               : (summands_.begin() + offset));
    const auto e = (offset_plus_count >= summands_.size()
                        ? summands_.end()
                        : (summands_.begin() + offset_plus_count));
    return ex<Sum>(b, e);
  }

  /// @tparam Filter a boolean predicate type, such `Filter(const ExprPtr&)`
  /// evaluates to true
  /// @param f an object of Filter type
  /// Selects elements {`e`} for which `f(e)` is true
  template <typename Filter>
  ExprPtr filter(Filter &&f) const {
    return ex<Sum>(summands_ | ranges::views::filter(f));
  }

  /// @return true if the number of factors is zero
  bool empty() const { return summands_.empty(); }

  /// @return the number of summands in a Sum
  std::size_t size() const { return summands_.size(); }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{ \\bigl(";
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
    result += L"\\bigr) }";
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

  ExprPtr clone() const override {
    auto cloned_summands =
        summands() | ranges::views::transform(
                         [](const ExprPtr &ptr) { return ptr->clone(); });
    return ex<Sum>(ranges::begin(cloned_summands),
                   ranges::end(cloned_summands));
  }

  /// @brief adjoint of a Sum is a sum of adjoints of its factors
  virtual void adjoint() override;

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
    hash_value_ =
        hash::range(ranges::begin(deref_summands), ranges::end(deref_summands));
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
    result.erase(0, 7);  // remove leading  "{ \bigl"
    result.replace(result.size() - 8, 8,
                   L")");  // replace trailing "\bigr) }" with ")"
    result = std::wstring(L"\\begin{align}\n& ") + result;
    // assume no inner sums
    size_t line_counter = 0;
    size_t term_counter = 0;
    std::wstring::size_type pos = 0;
    std::wstring::size_type plus_pos = 0;
    std::wstring::size_type minus_pos = 0;
    bool last_pos_has_plus = false;
    bool have_next_term = true;
    auto insert_into_result_at = [&](std::wstring::size_type at,
                                     const auto &str) {
      assert(pos != std::wstring::npos);
      result.insert(at, str);
      const auto str_nchar = std::size(str) - 1;  // neglect end-of-string
      pos += str_nchar;
      if (plus_pos != std::wstring::npos) plus_pos += str_nchar;
      if (minus_pos != std::wstring::npos) minus_pos += str_nchar;
      if (pos != plus_pos) assert(plus_pos == result.find(L" + ", plus_pos));
      if (pos != minus_pos) assert(minus_pos == result.find(L" - ", minus_pos));
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

template <typename Sequence>
std::decay_t<Sequence> clone(Sequence &&exprseq) {
  auto cloned_seq = exprseq | ranges::views::transform([](const ExprPtr &ptr) {
                      return ptr ? ptr->clone() : nullptr;
                    });
  return std::decay_t<Sequence>(ranges::begin(cloned_seq),
                                ranges::end(cloned_seq));
}

// finish off ExprPtr members that depend on Expr

template <typename T>
bool ExprPtr::is() const {
  return as_shared_ptr()->is<T>();
}

template <typename T>
const T &ExprPtr::as() const {
  return as_shared_ptr()->as<T>();
}

template <typename T>
T &ExprPtr::as() {
  return as_shared_ptr()->as<T>();
}

}  // namespace sequant

#include "expr_operator.hpp"

#include "expr_algorithm.hpp"

#endif  // SEQUANT_EXPR_HPP
