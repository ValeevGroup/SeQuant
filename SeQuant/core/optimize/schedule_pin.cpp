#include <SeQuant/core/optimize/schedule_pin.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace sequant::opt::schedule {

namespace {

namespace v1 = sequant::io::serialization::v1;
using eval::schedule::AxisRecord;
using eval::schedule::ModeKind;
using eval::schedule::NodeAnnotationRecord;
using eval::schedule::TermRecord;
using eval::schedule::TreeShape;

[[noreturn]] void fail(std::string const& key, std::string const& what) {
  throw Exception("schedule pin: term " + key + ": " + what);
}

std::string expr_string(ExprPtr const& e) { return toUtf8(v1::to_string(e)); }

/// the scalar spelled the way the V1 serializer prints a Constant
std::string scalar_string(Constant::scalar_type const& s) {
  return expr_string(ex<Constant>(s));
}

Constant::scalar_type parse_scalar(std::string const& key,
                                   std::string const& text) {
  ExprPtr e;
  try {
    e = v1::from_string<ExprPtr>(std::string_view{text});
  } catch (std::exception const&) {
    fail(key, "malformed scalar '" + text + "'");
  }
  if (!e || !e->is<Constant>()) fail(key, "malformed scalar '" + text + "'");
  return e->as<Constant>().value();
}

/// only a Constant, a lone Tensor, or a Product of Tensors can be pinned
void require_pure(ExprPtr const& canonical, std::string const& key) {
  if (canonical->is<Constant>() || canonical->is<Tensor>()) return;
  if (canonical->is<Product>()) {
    bool pure = true;
    for (auto const& f : canonical->as<Product>().factors())
      pure = pure && f->is<Tensor>();
    if (pure) return;
  }
  fail(key, "only tensor products can be pinned");
}

/// slot indices of a tensor factor and their proto indices
template <typename F>
void for_each_index(ExprPtr const& factor, F&& f) {
  for (auto const& ix : factor->as<Tensor>().const_slots()) {
    if (!ix.nonnull()) continue;
    f(ix);
    for (auto const& p : ix.proto_indices()) f(p);
  }
}

/// Left-first post-order walk over \p tree that mirrors \c binarize: an inner
/// node consumes an annotation (calls \p on_node with the ordinals of its
/// subtree's leaves) iff neither child evaluates to a scalar. A leaf is a
/// tensor; an inner node is a scalar iff none of its indices is named, i.e.
/// external to the term or used outside the subtree (the same named-index set
/// \c binarize canonicalizes the node's network with).
class NodeWalker {
 public:
  explicit NodeWalker(std::vector<ExprPtr> const& factors) : factors_(factors) {
    for (auto const& f : factors_)
      for_each_index(f, [this](Index const& ix) { ++total_[ix]; });
  }

  template <typename F>
  void walk(TreeShape const& t, F&& on_node) {
    std::vector<std::size_t> leaves;
    walk_impl(t, on_node, leaves);
  }

 private:
  std::vector<ExprPtr> const& factors_;
  std::map<Index, std::size_t> total_;

  bool is_scalar(std::vector<std::size_t> const& leaves) const {
    std::map<Index, std::size_t> here;
    for (auto o : leaves)
      for_each_index(factors_.at(o), [&here](Index const& ix) { ++here[ix]; });
    for (auto const& [ix, n] : here) {
      auto const N = total_.at(ix);
      if (n < N || N == 1) return false;
    }
    return true;
  }

  /// \return whether the subtree evaluates to a scalar
  template <typename F>
  bool walk_impl(TreeShape const& t, F& on_node,
                 std::vector<std::size_t>& leaves) {
    if (t.leaf) {
      leaves.push_back(*t.leaf);
      return false;
    }
    std::vector<std::size_t> sub;
    bool const ls = walk_impl(*t.left, on_node, sub);
    bool const rs = walk_impl(*t.right, on_node, sub);
    if (!ls && !rs) on_node(sub);
    leaves.insert(leaves.end(), sub.begin(), sub.end());
    return is_scalar(sub);
  }
};

ModeKind to_kind(BatchModeType t) {
  return t == BatchModeType::Contracted ? ModeKind::Contracted
                                        : ModeKind::External;
}

BatchModeType from_kind(ModeKind k) {
  return k == ModeKind::Contracted ? BatchModeType::Contracted
                                   : BatchModeType::External;
}

std::vector<AxisRecord> to_records(
    container::svector<std::pair<Index, BatchModeType>> const& axes) {
  std::vector<AxisRecord> result;
  result.reserve(axes.size());
  for (auto const& [ix, t] : axes)
    result.push_back({toUtf8(ix.full_label()), to_kind(t)});
  return result;
}

/// \throws if the leaves of \p tree are not a permutation of [0, n)
void require_permutation(TreeShape const& tree, std::size_t n,
                         std::string const& key) {
  std::vector<bool> seen(n, false);
  std::size_t count = 0;
  auto visit = [&](auto&& self, TreeShape const& t) -> void {
    if (t.leaf) {
      if (*t.leaf >= n || seen[*t.leaf])
        fail(key, "tree leaf " + std::to_string(*t.leaf) +
                      " is out of range or repeated");
      seen[*t.leaf] = true;
      ++count;
      return;
    }
    if (!t.left || !t.right) fail(key, "tree node lacks a child");
    self(self, *t.left);
    self(self, *t.right);
  };
  visit(visit, tree);
  if (count != n) fail(key, "tree does not use every factor");
}

}  // namespace

ExprPtr canonical_term(ExprPtr const& term) {
  return canonicalize(flatten(term->clone()));
}

std::string term_key(ExprPtr const& canonical) {
  return expr_string(canonical);
}

std::vector<ExprPtr> canonical_factors(ExprPtr const& canonical) {
  if (canonical->is<Product>()) {
    auto const& fs = canonical->as<Product>().factors();
    return std::vector<ExprPtr>(fs.begin(), fs.end());
  }
  if (canonical->is<Constant>()) return {};
  return {canonical};
}

TermRecord record_term(ExprPtr const& canonical, ExprPtr const& optimized,
                       container::vector<NodeBatchAnnotation> const& nodes) {
  TermRecord rec;
  rec.key = term_key(canonical);
  require_pure(canonical, rec.key);
  auto const factors = canonical_factors(canonical);
  for (auto const& f : factors) rec.factors.push_back(expr_string(f));

  if (canonical->is<Constant>()) {
    if (!nodes.empty()) fail(rec.key, "a scalar term has no contraction nodes");
    rec.scalar = scalar_string(canonical->as<Constant>().value());
    return rec;
  }

  // factor -> ordinal: by ExprPtr identity first, then by == for a factor the
  // optimizer cloned; each ordinal is used once (identical factors)
  std::vector<bool> used(factors.size(), false);
  auto ordinal_of = [&](ExprPtr const& f) -> std::size_t {
    for (std::size_t i = 0; i < factors.size(); ++i)
      if (!used[i] && factors[i].get() == f.get()) return used[i] = true, i;
    for (std::size_t i = 0; i < factors.size(); ++i)
      if (!used[i] && *factors[i] == *f) return used[i] = true, i;
    fail(rec.key,
         "optimized factor " + expr_string(f) + " is not a factor of the term");
  };
  auto to_shape = [&](auto&& self, ExprPtr const& e, bool top) -> TreeShape {
    if (e->is<Tensor>()) return TreeShape{ordinal_of(e)};
    if (!e->is<Product>())
      fail(rec.key, "optimized output contains a non-tensor " + expr_string(e));
    auto const& p = e->as<Product>();
    if (!top && p.scalar() != 1)
      fail(rec.key, "optimized output has a scalar on an inner node");
    if (p.factors().size() == 1) return self(self, p.factor(0), false);
    if (p.factors().size() != 2)
      fail(rec.key, "optimized output is not a binary tree");
    TreeShape node;
    node.left = std::make_unique<TreeShape>(self(self, p.factor(0), false));
    node.right = std::make_unique<TreeShape>(self(self, p.factor(1), false));
    return node;
  };
  rec.tree = to_shape(to_shape, optimized, true);
  rec.scalar = optimized->is<Product>()
                   ? scalar_string(optimized->as<Product>().scalar())
                   : scalar_string(Constant::scalar_type{1});
  for (bool u : used)
    if (!u) fail(rec.key, "optimized output does not use every factor");

  std::size_t k = 0;
  NodeWalker(factors).walk(*rec.tree, [&](std::vector<std::size_t> const&) {
    if (k >= nodes.size())
      fail(rec.key, "fewer annotations than contraction nodes");
    auto const& n = nodes[k++];
    rec.nodes.push_back(NodeAnnotationRecord{
        to_records(n.axes), to_records(n.opened_here), n.order_aware});
  });
  if (k != nodes.size())
    fail(rec.key, "more annotations than contraction nodes");
  return rec;
}

PinnedTerm rebuild_pinned_term(TermRecord const& rec,
                               ExprPtr const& canonical) {
  require_pure(canonical, rec.key);
  auto const factors = canonical_factors(canonical);
  std::vector<std::string> strings;
  for (auto const& f : factors) strings.push_back(expr_string(f));
  if (strings != rec.factors)
    fail(rec.key, "factors differ from those of " + term_key(canonical));

  auto const scalar = parse_scalar(rec.key, rec.scalar);
  PinnedTerm result;
  if (!rec.tree) {
    if (!factors.empty()) fail(rec.key, "a tensor term needs a tree");
    if (!rec.nodes.empty())
      fail(rec.key, "a scalar term has no contraction nodes");
    result.expr = ex<Constant>(scalar);
    return result;
  }
  require_permutation(*rec.tree, factors.size(), rec.key);

  auto build = [&](auto&& self, TreeShape const& t) -> ExprPtr {
    if (t.leaf) return factors[*t.leaf];
    return ex<Product>(1,
                       ExprPtrList{self(self, *t.left), self(self, *t.right)},
                       Product::Flatten::No);
  };
  if (rec.tree->leaf) {
    if (canonical->is<Tensor>()) {
      if (scalar != 1) fail(rec.key, "a lone tensor carries no scalar");
      result.expr = factors[*rec.tree->leaf];
    } else {
      result.expr = ex<Product>(scalar, ExprPtrList{factors[*rec.tree->leaf]},
                                Product::Flatten::No);
    }
  } else {
    result.expr = build(build, *rec.tree);
    result.expr->as<Product>().scale(scalar);
  }

  auto resolve = [&](std::vector<AxisRecord> const& axes,
                     std::map<std::string, Index> const& by_label,
                     std::size_t k) {
    container::svector<std::pair<Index, BatchModeType>> out;
    for (auto const& a : axes) {
      auto it = by_label.find(a.label);
      if (it == by_label.end())
        fail(rec.key, "axis label " + a.label +
                          " is not an index of contraction node " +
                          std::to_string(k));
      out.emplace_back(it->second, from_kind(a.kind));
    }
    return out;
  };
  std::size_t k = 0;
  NodeWalker(factors).walk(
      *rec.tree, [&](std::vector<std::size_t> const& leaves) {
        if (k >= rec.nodes.size())
          fail(rec.key, "fewer node records than contraction nodes");
        std::map<std::string, Index> by_label;
        for (auto o : leaves)
          for_each_index(factors[o], [&by_label](Index const& ix) {
            by_label.emplace(toUtf8(ix.full_label()), ix);
          });
        auto const& n = rec.nodes[k];
        NodeBatchAnnotation ann;
        ann.axes = resolve(n.axes, by_label, k);
        ann.opened_here = resolve(n.opened_here, by_label, k);
        ann.order_aware = n.order_aware;
        result.nodes.push_back(std::move(ann));
        ++k;
      });
  if (k != rec.nodes.size())
    fail(rec.key, "more node records than contraction nodes");
  return result;
}

}  // namespace sequant::opt::schedule
