#include <SeQuant/core/expressions/expr_holder.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <memory>

namespace sequant {

Power::Power(const ExprPtr& base, exponent_type exponent)
    : Power(ExprHolder(base), std::move(exponent)) {}

Power::Power(ExprHolder base, exponent_type exponent)
    : base_{std::move(base)}, exponent_{std::move(exponent)} {
  SEQUANT_ASSERT(base_->is<Constant>() || base_->is<Variable>());
  // 0^n is defined only for n >= 0 (0^0 = 1 by convention)
  SEQUANT_ASSERT(!base_->is<Constant>() || !base_->as<Constant>().is_zero() ||
                 exponent_ >= 0);
}

const ExprHolder& Power::base() const { return base_; }

const Power::exponent_type& Power::exponent() const { return exponent_; }

bool Power::conjugated() const { return conjugated_; }

void Power::conjugate() {
  conjugated_ = !conjugated_;
  reset_hash_value();
}

bool Power::is_zero() const {
  return exponent_ > 0 && base_->is<Constant>() &&
         base_->as<Constant>().is_zero();
}

std::wstring Power::to_latex() const { return io::latex::to_string(*this); }

template <expr_holder E>
void flatt_impl(E& expr) {
  auto create_constant = [](const auto& val) {
    if constexpr (std::same_as<std::remove_cvref_t<E>, ExprHolder>) {
      return Constant(val);
    } else {
      return ex<Constant>(val);
    }
  };

  const auto& pw = expr->template as<Power>();

  // b^1 = b and conjugate if needed
  if (pw.exponent() == 1) {
    auto lifted = pw.base().copy();
    if (pw.conjugated()) lifted->adjoint();
    expr = std::move(lifted);
    return;
  }
  // b^0 = 1 for any base (the ctor rejects 0^(negative)
  if (pw.exponent() == 0) {
    expr = create_constant(Constant::scalar_type{1});
    return;
  }
  if (!pw.base()->template is<Constant>()) return;

  using scalar_type = Constant::scalar_type;
  const auto& base_val = pw.base()->template as<Constant>().value();

  // 1^k = 1 for any rational k.
  if (base_val == scalar_type{1}) {
    expr = create_constant(scalar_type{1});
    return;
  }
  // 0^k for k < 0 is undefined; the ctor rejects it, and a disabled assert
  // must not let it fold to 0 or divide by zero below
  if (Constant::is_zero(base_val) && pw.exponent() < 0) return;

  scalar_type value{0};
  if (denominator(pw.exponent()) == 1) {
    // integer exponent, any Constant base: exp-by-squaring; for negative
    // exponents we power the magnitude and invert at the end
    auto exp_nr = numerator(pw.exponent());
    const bool negate = exp_nr < 0;
    if (negate) exp_nr = -exp_nr;
    // refuse a value too large to build: with the base written as
    // (n_re + i n_im) / d over a common denominator d, every part of the
    // value is bounded by (|n_re| + |n_im|)^k and every denominator by d^k
    const auto& d_re = denominator(base_val.real());
    const auto& d_im = denominator(base_val.imag());
    const intmax_t common_denominator = d_re * d_im;
    const intmax_t numerator_bound = abs(numerator(base_val.real())) * d_im +
                                     abs(numerator(base_val.imag())) * d_re;
    unsigned base_bits = msb(common_denominator) + 1;
    if (numerator_bound != 0)
      base_bits = std::max<unsigned>(base_bits, msb(numerator_bound) + 1);
    if (exp_nr * base_bits > default_max_folded_bits) return;
    value = scalar_type{1};
    scalar_type b = base_val;
    while (exp_nr > 0) {
      if (exp_nr % 2 != 0) value *= b;
      exp_nr /= 2;
      if (exp_nr > 0) b *= b;
    }
    if (negate) value = scalar_type{1} / value;
  } else if (Constant::is_zero(base_val)) {
    // 0^k = 0 for the k > 0 the ctor admits (k = 0 was handled above)
  } else if (base_val.imag() == 0 && base_val.real() > 0) {
    // a rational power of a positive rational folds iff it has no root
    const auto form = simplest_radical_form({{base_val.real(), pw.exponent()}});
    if (!form || form->root_order != 1) return;
    value = scalar_type{form->rational_part};
  } else {
    return;
  }

  if (pw.conjugated()) value = conj(value);
  expr = create_constant(std::move(value));
}

void Power::flatten(ExprPtr& expr) {
  if (!expr || !expr->is<Power>()) return;

  flatt_impl(expr);
}

void Power::flatten(ExprHolder& expr) {
  if (!expr->is<Power>()) return;

  flatt_impl(expr);
}

Expr::type_id_type Power::type_id() const { return get_type_id<Power>(); }

bool Power::is_scalar() const { return true; }

void Power::adjoint() { conjugate(); }

Power& Power::operator*=(const Expr& that) {
  // b^e1 *= b^e2  ->  b^(e1+e2)
  if (that.is<Power>()) {
    const auto& other = that.as<Power>();
    if (conjugated_ == other.conjugated_ && *base_ == *other.base_) {
      exponent_ += other.exponent_;
      reset_hash_value();
      return *this;
    }
  }
  // (b^e)* *= b*  ->  (b^(e+1))*
  else if (base_->is<Variable>() && that.is<Variable>()) {
    // check effective conjugation of Variable in this and that, if valid
    // operation iff they match
    const auto& base_var = base_->as<Variable>();
    const auto& that_var = that.as<Variable>();
    if (base_var.label() == that_var.label() &&
        (base_var.conjugated() ^ conjugated_) == that_var.conjugated()) {
      exponent_ += rational{1};
      reset_hash_value();
      return *this;
    }
  }
  // C^e *= C  ->  C^(e+1)
  else if (!conjugated_ && *base_ == that) {
    exponent_ += rational{1};
    reset_hash_value();
    return *this;
  }
  throw Exception("Power::operator*=(that): not valid for that");
}

std::unique_ptr<Expr> Power::unique_copy() const {
  auto copy = std::make_unique<Power>(base_.copy(), exponent_);
  if (conjugated_) copy->as<Power>().conjugate();
  copy->copy_canonical_mark(*this);
  return copy;
}

std::unique_ptr<Expr> Power::unique_move() && {
  return std::make_unique<Power>(std::move(*this));
}

Expr::hash_type Power::memoizing_hash() const {
  auto compute_hash = [this]() {
    if (exponent_ == 1 && !conjugated_) return hash::value(*base_);
    auto val = hash::value(*base_);
    hash::combine(val, hash::value(exponent_));
    hash::combine(val, conjugated_);
    return val;
  };

  if (!hash_value_) {
    hash_value_ = compute_hash();
  } else {
    SEQUANT_ASSERT(*hash_value_ == compute_hash());
  }
  return *hash_value_;
}

bool Power::static_equal(const Expr& that) const {
  const auto& other = static_cast<const Power&>(that);
  return exponent_ == other.exponent_ && conjugated_ == other.conjugated_ &&
         *base_ == *other.base_;
}

bool Power::static_less_than(const Expr& that) const {
  const auto& other = static_cast<const Power&>(that);
  if (*base_ != *other.base_) return *base_ < *other.base_;
  if (exponent_ != other.exponent_) return exponent_ < other.exponent_;
  return conjugated_ < other.conjugated_;
}

std::weak_ordering Power::static_lexicographic_compare(
    const Expr& that, const LexicographicCompareOptions& opts) const {
  const auto& other = static_cast<const Power&>(that);
  if (auto cmp = base_->lexicographic_compare(*other.base_, opts); cmp != 0)
    return cmp;
  if (auto cmp = std::compare_weak_order_fallback(exponent_, other.exponent_);
      cmp != 0)
    return cmp;
  return conjugated_ <=> other.conjugated_;
}
}  // namespace sequant
