#include <SeQuant/core/optimize/fusion.hpp>

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/reverse.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/zip.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

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

namespace {

/// \return whether \p e is a Tensor, or a Product of such trees, none of whose
///         indices has proto indices
bool is_plain_tensor_tree(Expr const& e) {
  if (e.is<Tensor>())
    return ranges::all_of(
        e.as<Tensor>().const_braketaux_indices(),
        [](Index const& ix) { return !ix.has_proto_indices(); });
  return e.is<Product>() &&
         ranges::all_of(e.as<Product>().factors(), [](ExprPtr const& f) {
           return is_plain_tensor_tree(*f);
         });
}

/// \return whether \p e is a Product of two plain tensor trees
bool is_binary_contraction(Expr const& e) {
  return e.is<Product>() && e.as<Product>().size() == 2 &&
         is_plain_tensor_tree(e);
}

/// Calls \p f on each slot index of each tensor of the tensor tree \p e.
template <typename F>
void for_each_slot(Expr const& e, F&& f) {
  if (e.is<Tensor>()) {
    for (auto const& ix : e.as<Tensor>().const_braketaux_indices()) f(ix);
    return;
  }
  for (auto const& x : e.as<Product>().factors()) for_each_slot(*x, f);
}

container::set<Index> free_indices(ExprPtr const& e) {
  auto const g = get_unique_indices(e);
  container::set<Index> result;
  result.insert(g.bra.begin(), g.bra.end());
  result.insert(g.ket.begin(), g.ket.end());
  result.insert(g.aux.begin(), g.aux.end());
  return result;
}

/// Key of the shared factor F of a top contraction F X: F's structure with
/// the external indices kept and the other indices numbered in order of first
/// appearance, followed by X's free indices numbered the same way. Two
/// contractions with equal keys have the same F up to a renaming of its
/// non-external indices, and that renaming maps the X's free indices onto each
/// other.
struct FactorKey {
  std::wstring key;
  /// the numbered indices, in numbering order
  container::vector<Index> order;
};

class KeyBuilder {
 public:
  explicit KeyBuilder(container::set<Index> const& external)
      : external_(external) {}

  void tree(Expr const& e) {
    if (e.is<Tensor>()) {
      auto const& t = e.as<Tensor>();
      key_.key += t.label();
      key_.key += L'{' + std::to_wstring(static_cast<int>(t.symmetry())) +
                  L',' +
                  std::to_wstring(static_cast<int>(t.braket_symmetry())) +
                  L',' + std::to_wstring(static_cast<int>(t.column_symmetry()));
      auto slots = [this](auto const& indices) {
        key_.key += L';';
        for (auto const& ix : indices) index(ix);
      };
      slots(t.bra());
      slots(t.ket());
      slots(t.aux());
      key_.key += L'}';
      return;
    }
    auto const& p = e.as<Product>();
    key_.key += L'(' + ex<Constant>(p.scalar())->to_latex();
    for (auto const& f : p.factors()) {
      key_.key += L' ';
      tree(*f);
    }
    key_.key += L')';
  }

  void free_indices_of(ExprPtr const& x) {
    std::vector<std::wstring> tokens;
    for (auto const& ix : free_indices(x)) {
      std::wstring saved;
      std::swap(saved, key_.key);
      index(ix);
      tokens.push_back(std::move(key_.key));
      key_.key = std::move(saved);
    }
    std::ranges::sort(tokens);
    key_.key += L'|';
    for (auto const& t : tokens) key_.key += t + L',';
  }

  FactorKey release() { return std::move(key_); }

 private:
  container::set<Index> const& external_;
  container::map<Index, std::size_t> ids_;
  FactorKey key_;

  /// the basis of \p ix, everything Index equality compares but the ordinal
  /// and the proto indices
  static std::wstring basis_token(Index const& ix) {
    auto const& b = ix.basis();
    auto h = hash::value(ix.space());
    if (b.has_basis_instance()) {
      hash::combine(h, *b.basis_instance());
      if (b.has_name()) hash::combine(h, b.name());
    }
    return std::wstring(ix.space().base_key()) + L'<' + std::to_wstring(h) +
           L'>';
  }

  void index(Index const& ix) {
    if (external_.contains(ix)) {
      key_.key += L'E';
      key_.key += ix.full_label();
    } else {
      auto const [it, inserted] = ids_.try_emplace(ix, ids_.size());
      if (inserted) key_.order.push_back(ix);
      key_.key += basis_token(ix) + L'#' + std::to_wstring(it->second);
    }
    key_.key += L',';
  }
};

/// \pre is_binary_contraction(*summand)
FactorKey factor_key(ExprPtr const& summand, std::size_t side) {
  auto const& p = summand->as<Product>();
  auto const external = free_indices(summand);
  KeyBuilder b(external);
  b.tree(*p.factor(side));
  b.free_indices_of(p.factor(1 - side));
  return b.release();
}

/// One summand of a fusion group: its position and the side of its shared
/// factor.
struct Member {
  std::size_t pos;
  std::size_t side;
};

/// \return a copy of the tensor tree \p e with its indices replaced per \p map
ExprPtr renamed(ExprPtr const& e, container::map<Index, Index> const& map) {
  if (e->is<Tensor>()) {
    auto result = e->clone();
    auto& t = result->as<Tensor>();
    // transform_indices skips tagged indices and tags the ones it replaces
    t.reset_tags();
    t.transform_indices(map);
    t.reset_tags();
    return result;
  }
  auto const& p = e->as<Product>();
  Product result;
  result.scale(p.scalar());
  for (auto const& f : p.factors())
    result.append(1, renamed(f, map), Product::Flatten::No);
  return ex<Product>(std::move(result));
}

Sum::summands_type fuse(Sum::summands_type const& summands,
                        container::set<Index> const& outer);

/// Fuses \p group (two or more members with equal factor keys \p keys).
/// \p outer holds the indices of the enclosing expression, which renamed
/// dummies must avoid.
ExprPtr fuse_group(Sum::summands_type const& summands,
                   std::vector<std::array<FactorKey, 2>> const& keys,
                   std::vector<Member> const& group,
                   container::set<Index> const& outer) {
  auto const [pos0, side0] = group.front();
  auto const& shared = summands[pos0]->as<Product>().factor(side0);

  // a renamed X must not reuse an index of the shared factor or of the
  // external ones for its own dummies
  container::set<Index> taken = free_indices(summands[pos0]);
  for_each_slot(*shared, [&](Index const& ix) { taken.insert(ix); });
  taken.insert(outer.begin(), outer.end());
  // fresh indices avoid every index of the group and of the enclosing
  // expression
  container::set<Index> used = outer;
  for (auto const& m : group)
    for_each_slot(*summands[m.pos], [&](Index const& ix) { used.insert(ix); });
  IndexFactory fresh([&used](Index const& ix) { return !used.contains(ix); },
                     1);

  auto const& order0 = keys[pos0][side0].order;
  Sum::summands_type xs;
  for (auto const& m : group) {
    auto const& p = summands[m.pos]->as<Product>();
    auto const& x0 = p.factor(1 - m.side);
    container::map<Index, Index> map;
    if (m.pos != pos0) {
      auto const& order = keys[m.pos][m.side].order;
      SEQUANT_ASSERT(order.size() == order0.size());
      container::set<Index> mapped(order.begin(), order.end());
      for (std::size_t i = 0; i < order.size(); ++i)
        if (order[i] != order0[i]) map.emplace(order[i], order0[i]);
      auto const external = free_indices(summands[m.pos]);
      for_each_slot(*x0, [&](Index const& ix) {
        if (mapped.contains(ix) || external.contains(ix) ||
            !taken.contains(ix) || map.contains(ix))
          return;
        auto const f = fresh.make(ix);
        used.insert(f);
        map.emplace(ix, f);
      });
    }
    auto x = renamed(x0, map);
    if (p.scalar() != 1) {
      if (x->is<Product>())
        x->as<Product>().scale(p.scalar());
      else
        x = ex<Product>(p.scalar(), ExprPtrList{x}, Product::Flatten::No);
    }
    xs.push_back(std::move(x));
  }

  auto sum = ex<Sum>(fuse(xs, taken), Sum::move_only_tag{});
  auto result = ex<Product>();
  result->as<Product>().append(1, shared->clone(), Product::Flatten::No);
  if (side0 == 0)
    result->as<Product>().append(1, sum, Product::Flatten::No);
  else
    result->as<Product>().prepend(1, sum, Product::Flatten::No);
  return result;
}

Sum::summands_type fuse(Sum::summands_type const& summands,
                        container::set<Index> const& outer) {
  auto const n = summands.size();
  std::vector<bool> eligible(n);
  std::vector<std::array<FactorKey, 2>> keys(n);
  container::map<std::wstring, std::size_t> count;
  for (std::size_t i = 0; i < n; ++i) {
    eligible[i] = is_binary_contraction(*summands[i]);
    if (!eligible[i]) continue;
    for (std::size_t side = 0; side < 2; ++side)
      keys[i][side] = factor_key(summands[i], side);
    ++count[keys[i][0].key];
    if (keys[i][1].key != keys[i][0].key) ++count[keys[i][1].key];
  }

  // each summand joins the group of its side with the most candidates
  container::map<std::wstring, std::vector<Member>> groups;
  std::vector<std::optional<std::wstring>> group_of(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (!eligible[i]) continue;
    std::size_t const side =
        count[keys[i][1].key] > count[keys[i][0].key] ? 1 : 0;
    groups[keys[i][side].key].push_back({i, side});
    group_of[i] = keys[i][side].key;
  }

  Sum::summands_type result;
  for (std::size_t i = 0; i < n; ++i) {
    auto const* g = group_of[i] ? &groups.at(*group_of[i]) : nullptr;
    if (!g || g->size() < 2)
      result.push_back(summands[i]);
    else if (g->front().pos == i)
      result.push_back(fuse_group(summands, keys, *g, outer));
  }
  return result;
}

}  // namespace

Sum::summands_type fuse_top_factors(Sum::summands_type const& summands) {
  return fuse(summands, {});
}

}  // namespace sequant::opt
