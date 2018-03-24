//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_EXPR_HPP
#define SEQUANT2_EXPR_HPP

#include <complex>
#include <vector>

namespace sequant2 {

class Expr {
 public:
  virtual ~Expr() = default;

  /// @throw std::logic_error if Expr's type does match the type of this
  virtual bool operator==(const Expr &) const {
    throw std::logic_error("Expr::operator== not implemented in this derived class");
  }

  virtual std::wstring to_latex() const {
    throw std::logic_error("Expr::to_latex not implemented in this derived class");
  }
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

}  // namespace sequant2

#endif //SEQUANT2_EXPR_HPP
