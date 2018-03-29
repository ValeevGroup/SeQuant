//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_EXPR_HPP
#define SEQUANT2_EXPR_HPP

#include <complex>
#include <memory>
#include <vector>

#include <range/v3/all.hpp>

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

  virtual ~Expr() = default;

  /// @throw std::logic_error if Expr's type does match the type of this
  virtual bool operator==(const Expr &) const {
    throw std::logic_error("Expr::operator== not implemented in this derived class");
  }

  /// @return the string representation of @c this
  virtual std::wstring to_latex() const {
    throw std::logic_error("Expr::to_latex not implemented in this derived class");
  }

  /// Canonicalizes @c this and returns the biproduct of canonicalization (e.g. phase)
  /// @return the biproduct of canonicalization
  virtual std::shared_ptr<Expr> canonicalize() {
    throw std::logic_error("Expr::canonicalize not implemented in this derived class");
  }

  auto begin_subexpr() {
    return range_type::begin();
  }

  auto end_subexpr() {
    return range_type::end();
  }

  Expr& expr() {
    return *this;
  }

 private:
  friend ranges::range_access;

  struct cursor
  {
    using value_type = std::shared_ptr<Expr>;

    cursor() = default;
    constexpr explicit cursor(std::shared_ptr<Expr>* subexpr_ptr) noexcept
        : ptr_{subexpr_ptr}
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
    const std::shared_ptr<Expr>& read() const
    {
      RANGES_EXPECT(ptr_);
      return *ptr_;
    }
    std::shared_ptr<Expr>& read()
    {
      RANGES_EXPECT(ptr_);
      return *ptr_;
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
  };

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

};

class Constant : public Expr {
 public:
  Constant() = default;
  virtual ~Constant() = default;
  Constant(const Constant&) = default;
  Constant(Constant&&) = default;
  Constant& operator=(const Constant&) = default;
  Constant& operator=(Constant&&) = default;

  bool operator==(const Expr &expr) const override{
    const auto *expr_cast_ptr = dynamic_cast<const Constant *>(&expr);
    if (!expr_cast_ptr)
      throw std::bad_cast();
    return *this == *expr_cast_ptr;
  }

  std::wstring to_latex() const override {
    return L"{" + std::to_wstring(value_.real()) + L"}";
  }

  virtual std::shared_ptr<Expr> canonicalize() override {
    return {};
  }

  auto value() const { return value_; }

 private:
  std::complex<double> value_;
};

class ScaledProduct : public Expr {
 public:
  ScaledProduct() = default;
  virtual ~ScaledProduct() = default;
  ScaledProduct(const ScaledProduct&) = default;
  ScaledProduct(ScaledProduct&&) = default;
  ScaledProduct& operator=(const ScaledProduct&) = default;
  ScaledProduct& operator=(ScaledProduct&&) = default;

  template<typename T>
  ScaledProduct &append(T scalar, std::shared_ptr<Expr> factor) {
    scalar_ *= scalar;
    factors_.push_back(factor);
    return *this;
  }

  const std::complex<double> &scalar() const { return scalar_; }
  const auto& factors() const { return factors_; }

  /// @return true if the number of factors is zero
  bool empty() const { return factors_.empty(); }

  /// @param expr an Expr object
  /// @return true if @c expr is identical to @c *this
  /// @throw std::bad_cast if Expr's type does match the type of this
  bool operator==(const Expr &expr) const override {
    const auto *expr_cast_ptr = dynamic_cast<const ScaledProduct *>(&expr);
    if (!expr_cast_ptr)
      throw std::bad_cast();
    return *this == *expr_cast_ptr;
  }

  /// @param other a ScaledProduct object
  /// @return true if @c other is identical to @c *this
  bool operator==(const ScaledProduct &other) const {
    using std::begin;
    using std::end;
    if (this->scalar() == other.scalar()) {
      return std::equal(begin(factors_), end(factors_), begin(other.factors_),
                        [](const auto &i1, const auto &i2) { return i1->operator==(*i2); });
    } else
      return false;
  }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    assert(scalar().imag() == 0.0);
    result += std::to_wstring(scalar().real());
    result += L" \\times ";
    for (const auto &i : factors())
      result += i->to_latex();
    result += L"}";
    return result;
  }

 private:
  std::complex<double> scalar_ = {1.0, 0.0};
  container::vector<std::shared_ptr<Expr>> factors_{};
};

class ScaledSum : public Expr {
 public:
  ScaledSum() = default;
  virtual ~ScaledSum() = default;
  ScaledSum(const ScaledSum&) = default;
  ScaledSum(ScaledSum&&) = default;
  ScaledSum& operator=(const ScaledSum&) = default;
  ScaledSum& operator=(ScaledSum&&) = default;

  template<typename T>
  ScaledSum &append(std::shared_ptr<Expr> summand) {
    summands_.push_back(summand);
    return *this;
  }

  const std::complex<double> &scalar() const { return scalar_; }
  const auto& summands() const { return summands_; }

  /// @return true if the number of factors is zero
  bool empty() const { return summands_.empty(); }

  /// @param expr an Expr object
  /// @return true if @c expr is identical to @c *this
  /// @throw std::bad_cast if Expr's type does match the type of this
  bool operator==(const Expr &expr) const override {
    const auto *expr_cast_ptr = dynamic_cast<const ScaledSum *>(&expr);
    if (!expr_cast_ptr)
      throw std::bad_cast();
    return *this == *expr_cast_ptr;
  }

  /// @param other a ScaledProduct object
  /// @return true if @c other is identical to @c *this
  bool operator==(const ScaledSum &other) const {
    using std::begin;
    using std::end;
    if (this->scalar() == other.scalar()) {
      return std::equal(begin(summands_), end(summands_), begin(other.summands_),
                        [](const auto &i1, const auto &i2) { return i1->operator==(*i2); });
    } else
      return false;
  }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    assert(scalar().imag() == 0.0);
    result += std::to_wstring(scalar().real());
    result += L" \\times \\left(";
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

 private:
  std::complex<double> scalar_ = {1.0, 0.0};
  container::vector<std::shared_ptr<Expr>> summands_{};
};

inline std::wstring to_latex(const Expr& expr) {
  return expr.to_latex();
}

/// Canonicalizes an Expr and replaces it as needed
/// @param[in,out] expr expression to be canonicalized; will be replaced if canonicalization is impure
inline void canonicalize(std::shared_ptr<Expr>& expr) {
  const auto biproduct = expr->canonicalize();
  const auto& biproduct_value = *biproduct;
  const auto& biproduct_typeid = typeid(biproduct_value);
  if (biproduct_typeid == typeid(Constant)) {
    const auto constant_ptr = std::dynamic_pointer_cast<Constant>(biproduct);
    auto new_expr = std::make_shared<ScaledProduct>();
    new_expr->append(constant_ptr->value(), expr);
    expr = new_expr;
  }
}

}  // namespace sequant2

#endif //SEQUANT2_EXPR_HPP
