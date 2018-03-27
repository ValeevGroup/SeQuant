//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_EXPR_HPP
#define SEQUANT2_EXPR_HPP

#include <complex>
#include <vector>

namespace sequant2 {

/// Base expression class
class Expr {
 public:
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
  virtual std::shared_ptr<Expr> canonicalize() { abort(); }
};

class ScaledProduct : public Expr {
 public:
  ScaledProduct() = default;
  virtual ~ScaledProduct() = default;

  template<typename T>
  ScaledProduct &append(T scalar, std::shared_ptr<Expr> factor) {
    scalar_ *= scalar;
    factors_.push_back(factor);
    return *this;
  }

  const std::complex<double> &scalar() const { return scalar_; }
  const std::vector<std::shared_ptr<Expr>> &factors() const { return factors_; }

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
  std::vector<std::shared_ptr<Expr>> factors_{};
};

class ScaledSum : public Expr {
 public:
  ScaledSum() = default;
  virtual ~ScaledSum() = default;

  template<typename T>
  ScaledSum &append(std::shared_ptr<Expr> summand) {
    summands_.push_back(summand);
    return *this;
  }

  const std::complex<double> &scalar() const { return scalar_; }
  const std::vector<std::shared_ptr<Expr>> &summands() const { return summands_; }

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
  std::vector<std::shared_ptr<Expr>> summands_{};
};

inline std::wstring to_latex(const Expr& expr) {
  return expr.to_latex();
}

/// Canonicalizes @c expr and returns the biproduct of canonicalization (e.g. phase)
/// @return the biproduct of canonicalization
inline std::shared_ptr<Expr> canonicalize(Expr& expr) {
  return expr.canonicalize();
}

}  // namespace sequant2

#endif //SEQUANT2_EXPR_HPP
