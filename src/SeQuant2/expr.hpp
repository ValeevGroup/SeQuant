//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_EXPR_HPP
#define SEQUANT2_EXPR_HPP

#include <complex>
#include <iostream>
#include <memory>
#include <vector>

#include <range/v3/all.hpp>

#include <boost/optional.hpp>
#include <boost/functional/hash.hpp>
#include <boost/callable_traits.hpp>

#include "latex.hpp"
#include "vector.hpp"

namespace sequant2 {

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

  virtual ~Expr() = default;

  /// @return true if this is a leaf
  bool is_atom() const {
    return ranges::empty(*this);
  }

  /// @return the string representation of @c this
  virtual std::wstring to_latex() const {
    throw std::logic_error("Expr::to_latex not implemented in this derived class");
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
  /// @note always returns 0 unless this derived class overrides Expr::memoized_hash
  /// @return the hash value for this Expr
  hash_type hash_value() const {
    return memoizing_hash();
  };

  /// Computes and returns the derived type identifier
  /// @note this function must be overridden in the derived class
  /// @sa Expr::get_type_id
  /// @return the hash value for this Expr
  virtual type_id_type type_id() const =0;

  /// @tparam T Expr or a class derived from Expr
  /// @return true if this is equal to that
  /// @note the derived class must implement Expr::static_compare
  template <typename T>
  std::enable_if_t<std::is_base_of<Expr,T>::value, bool>
      operator==(const T& that) const {
    if (this->type_id() != that.type_id())
      return false;
    else
      return this->static_compare(static_cast<const Expr&>(that));
  }

  /// @tparam T Expr or a class derived from Expr
  /// @return true if this is equal to that
  /// @note the derived class must implement Expr::static_compare
  template <typename T>
  std::enable_if_t<std::is_base_of<Expr,T>::value, bool>
  operator!=(const T& that) const {
     return ! operator==(that);
  }

  /// @return (unique) type id of class T
  template <typename T> static  type_id_type get_type_id() {
    static type_id_type type_id = get_next_type_id();
    return type_id;
  };

 private:
  friend ranges::range_access;

 protected:
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

  /// @return the cursor for the beginning of the range (must overridden in a derived Expr that has subexpressions)
  virtual cursor begin_cursor()
  {
    return cursor{};
  }
  /// @return the cursor for the end of the range (must overridden in a derived Expr that has subexpressions)
  virtual cursor end_cursor()
  {
    return cursor{};
  }

  /// @return the cursor for the beginning of the range (must overridden in a derived Expr that has subexpressions)
  virtual cursor begin_cursor() const
  {
    return cursor{};
  }
  /// @return the cursor for the end of the range (must overridden in a derived Expr that has subexpressions)
  virtual cursor end_cursor() const
  {
    return cursor{};
  }

  mutable boost::optional<hash_type> hash_value_;  // not initialized by default
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

 private:
  /// @return returns next type id in the grand class list
  static type_id_type get_next_type_id() {
    static std::atomic<type_id_type> grand_type_id = 0;
    return ++grand_type_id;
  };

  /// @param that an Expr object
  /// @note @c that is guaranteed to be of same type as @c *this, hence can be statically cast
  /// @return true if @c that is equivalent to *this
  virtual bool static_compare(const Expr& that) const =0;
};

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

  virtual std::shared_ptr<Expr> canonicalize() override {
    return {};
  }

  type_id_type type_id() const override {
    return get_type_id<Constant>();
  };

  std::shared_ptr<Expr> clone() const override {
    return make<Constant>(this->value());
  }

 private:
  std::complex<double> value_;

  hash_type memoizing_hash() const override {
    hash_value_ = boost::hash_value(value_);
    return *hash_value_;
  }

  bool static_compare(const Expr& that) const override {
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
  Product(std::initializer_list<ExprPtr> factors) : factors_(std::move(factors)) {}

  /// construct a Product out of zero or more factors multiplied by a scalar
  /// @tparam T a numeric type; it must be able to multiply std::complex<double>
  /// @param scalar a scalar of type T
  /// @param factors an initializer list of factors
  template<typename T>
  Product(T scalar, std::initializer_list<ExprPtr> factors) : scalar_(std::move(scalar)), factors_(std::move(factors)) {}

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

  /// multiplies the product by @c scalar times @c factor
  template<typename T>
  Product &append(T scalar, ExprPtr factor) {
    scalar_ *= scalar;
    factors_.push_back(std::move(factor));
    reset_hash_value();
    return *this;
  }

  const std::complex<double> &scalar() const { return scalar_; }
  const auto& factors() const { return factors_; }

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

  type_id_type type_id() const override {
    return get_type_id<Product>();
  };

  std::shared_ptr<Expr> clone() const override {
    return make<Product>(factors().begin(), factors().end());
  }

 private:
  std::complex<double> scalar_ = {1.0, 0.0};
  container::vector<ExprPtr> factors_{};

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

  bool static_compare(const Expr& that) const override {
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
  Sum(std::initializer_list<ExprPtr> summands) {
    // use append to flatten out Sum summands
    for(auto& summand: summands) {
      append(std::move(summand));
    }
  }

  /// append a summand to the sum
  /// @param summand the summand
  Sum &append(ExprPtr summand) {
    if (summand->type_id() != Expr::get_type_id<Sum>())
      summands_.push_back(std::move(summand));
    else {  // this recursively flattens Sum summands
      for(auto& subsummand: *summand)
        this->append(subsummand);
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

  Expr::type_id_type type_id() const override {
    return Expr::get_type_id<Sum>();
  };


 private:
  container::vector<ExprPtr> summands_{};

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

  bool static_compare(const Expr& that) const override {
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

};  // namespace sequant2

#include "expr_operator.hpp"
#include "expr_algorithm.hpp"

#endif //SEQUANT2_EXPR_HPP
