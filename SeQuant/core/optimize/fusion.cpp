#include "fusion.hpp"
#include <range/v3/action.hpp>
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
  auto fac = container::vector<ExprPtr>{};

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

  if (lhs.scalar() == rhs.scalar())
    fac_prod.scale(lhs.scalar());
  else {
    lsmand_prod.scale(lhs.scalar());
    rsmand_prod.scale(rhs.scalar());
  }

  // f (a + b)

  auto f = lift_tensor(fac_prod);
  auto a = lift_tensor(lsmand_prod);
  auto b = lift_tensor(rsmand_prod);

  return ex<Product>(ExprPtrList{f, ex<Sum>(ExprPtrList{a, b})});
}

ExprPtr Fusion::fuse_right(Product const& lhs, Product const& rhs) {
  auto fac = container::vector<ExprPtr>{};

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

  if (lhs.scalar() == rhs.scalar())
    fac_prod.scale(lhs.scalar());
  else {
    lsmand_prod.scale(lhs.scalar());
    rsmand_prod.scale(rhs.scalar());
  }

  // (a + b) f

  auto a = lift_tensor(lsmand_prod);
  auto b = lift_tensor(rsmand_prod);
  auto f = lift_tensor(fac_prod);

  return ex<Product>(ExprPtrList{ex<Sum>(ExprPtrList{a, b}), f});
}

}  // namespace sequant::opt
