#include <SeQuant/core/optimize/fusion.hpp>

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm.hpp>
#include <range/v3/iterator.hpp>
#include <range/v3/view.hpp>

namespace sequant::opt {

using ranges::views::drop;
using ranges::views::reverse;
using ranges::views::zip;

// convert a Product of single tensor and scalar == 1 into a tensor exprptr
auto lift_tensor = [](Product const& p) -> ExprPtr {
  return p.scalar() == 1 && p.size() == 1
             ? p.factor(0)
             : ex<Product>(Product{p.scalar(), p.factors().begin(),
                                   p.factors().end(), Product::Flatten::No});
};

Fusion::Fusion(Product const& lhs, Product const& rhs)
    : left_{fuse_left(lhs, rhs)}, right_{fuse_right(lhs, rhs)} {}

ExprPtr Fusion::left() const { return left_; }

ExprPtr Fusion::right() const { return right_; }

ExprPtr Fusion::fuse_left(Product const& lhs, Product const& rhs) {
  auto fac = container::svector<ExprPtr>{};

  for (auto&& [e1, e2] : zip(lhs.factors(), rhs.factors())) {
    if (e1 == e2)
      fac.push_back(e1);
    else
      break;
  }

  if (fac.empty()) return nullptr;

  auto lsmand = lhs.factors() | drop(fac.size());
  auto rsmand = rhs.factors() | drop(fac.size());

  auto fac_prod = Product{fac.begin(), fac.end()};
  auto lsmand_prod = Product{lsmand.begin(), lsmand.end()};
  auto rsmand_prod = Product{rsmand.begin(), rsmand.end()};

  SEQUANT_ASSERT(lhs.scalar().imag().is_zero() &&
                 rhs.scalar().imag().is_zero() &&
                 "Complex valued gcd not supported");
  auto scalars_fused = fuse_scalar(lhs.scalar().real(), rhs.scalar().real());

  fac_prod.scale(scalars_fused.at(0));
  lsmand_prod.scale(scalars_fused.at(1));
  rsmand_prod.scale(scalars_fused.at(2));

  // f (a + b)

  auto f = lift_tensor(fac_prod);
  auto a = lift_tensor(lsmand_prod);
  auto b = lift_tensor(rsmand_prod);

  return ex<Product>(ExprPtrList{f, ex<Sum>(ExprPtrList{a, b})});
}

ExprPtr Fusion::fuse_right(Product const& lhs, Product const& rhs) {
  auto fac = container::svector<ExprPtr>{};

  for (auto&& [e1, e2] :
       zip(lhs.factors() | reverse, rhs.factors() | reverse)) {
    if (e1 == e2)
      fac.push_back(e1);
    else
      break;
  }

  if (fac.empty()) return nullptr;

  ranges::reverse(fac);
  auto lsmand = lhs.factors() | reverse | drop(fac.size()) | reverse;
  auto rsmand = rhs.factors() | reverse | drop(fac.size()) | reverse;

  auto fac_prod = Product{fac.begin(), fac.end()};
  auto lsmand_prod = Product{lsmand.begin(), lsmand.end()};
  auto rsmand_prod = Product{rsmand.begin(), rsmand.end()};

  SEQUANT_ASSERT(lhs.scalar().imag().is_zero() &&
                 rhs.scalar().imag().is_zero() &&
                 "Complex valued gcd not supported");
  auto scalars_fused = fuse_scalar(lhs.scalar().real(), rhs.scalar().real());

  fac_prod.scale(scalars_fused.at(0));
  lsmand_prod.scale(scalars_fused.at(1));
  rsmand_prod.scale(scalars_fused.at(2));

  // (a + b) f

  auto a = lift_tensor(lsmand_prod);
  auto b = lift_tensor(rsmand_prod);
  auto f = lift_tensor(fac_prod);

  return ex<Product>(ExprPtrList{ex<Sum>(ExprPtrList{a, b}), f});
}

rational Fusion::gcd_rational(rational const& left, rational const& right) {
  auto&& r1 = left.real();
  auto&& r2 = right.real();
  auto&& n1 = numerator(r1);
  auto&& d1 = denominator(r1);
  auto&& n2 = numerator(r2);
  auto&& d2 = denominator(r2);

  auto num = gcd(n1 * d2, n2 * d1);
  return {num, d1 * d2};
}

std::array<rational, 3> Fusion::fuse_scalar(rational const& left,
                                            rational const& right) {
  auto fused = gcd_rational(left, right);
  rational left_fused = left / fused;
  rational right_fused = right / fused;
  if (left < 0 && right < 0) {
    fused *= -1;
    left_fused *= -1;
    right_fused *= -1;
  }
  return {fused, left_fused, right_fused};
}

}  // namespace sequant::opt
