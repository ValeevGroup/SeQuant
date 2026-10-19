#include <SeQuant/core/optimize/schedule_pin.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <mutex>
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

/// Left-first post-order walk over \p tree calling \p on_node with the
/// ordinals of each inner node's subtree leaves; one call per inner node.
///
/// Every inner node of a pure tensor product tree is a tensor x tensor node
/// and owns one annotation. This mirrors the producer exactly: the optimizer
/// emits #tensors - 1 entries, one per inner node, in this post-order
/// (single_term.hpp, opt_pure_product and single_term_opt). binarize stamps an
/// entry only on a node neither of whose children is a scalar
/// (eval_expr.cpp, binarize(Product) make_prod), and a child is a scalar only
/// when it is fully contracted with respect to the binarize-time external set
/// (left_to_right_binarization_indices, utility/indices.hpp, whose
/// `uncontract` comes from the ResultExpr LHS that a term does not carry).
/// When that happens with a non-empty node_batch_axes, binarize's
/// node_counter == node_batch_axes.size() assertion fires; so a valid
/// annotation vector always has one entry per inner node, and that is the
/// rule used here.
template <typename F>
void walk_nodes(TreeShape const& t, F&& on_node,
                std::vector<std::size_t>& leaves) {
  if (t.leaf) {
    leaves.push_back(*t.leaf);
    return;
  }
  std::vector<std::size_t> sub;
  walk_nodes(*t.left, on_node, sub);
  walk_nodes(*t.right, on_node, sub);
  on_node(sub);
  leaves.insert(leaves.end(), sub.begin(), sub.end());
}

/// number of inner (tensor x tensor) nodes of \p t
std::size_t count_nodes(TreeShape const& t) {
  return t.leaf ? 0 : 1 + count_nodes(*t.left) + count_nodes(*t.right);
}

[[noreturn]] void fail_count(std::string const& key, std::size_t n,
                             std::size_t k) {
  fail(key, std::to_string(n) + " annotations supplied but the tree has " +
                std::to_string(k) +
                " tensor-times-tensor nodes (one annotation per node, or "
                "none at all)");
}

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

  // an empty annotation vector (non-batched objective) records the tree only
  if (nodes.empty()) return rec;
  if (nodes.size() != count_nodes(*rec.tree))
    fail_count(rec.key, nodes.size(), count_nodes(*rec.tree));
  std::size_t k = 0;
  std::vector<std::size_t> all;
  walk_nodes(
      *rec.tree,
      [&](std::vector<std::size_t> const&) {
        auto const& n = nodes[k++];
        rec.nodes.push_back(NodeAnnotationRecord{
            to_records(n.axes), to_records(n.opened_here), n.order_aware});
      },
      all);
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
  // a record without node annotations (non-batched objective) rebuilds the
  // tree only
  if (rec.nodes.empty()) return result;
  if (rec.nodes.size() != count_nodes(*rec.tree))
    fail_count(rec.key, rec.nodes.size(), count_nodes(*rec.tree));
  std::size_t k = 0;
  std::vector<std::size_t> all;
  walk_nodes(
      *rec.tree,
      [&](std::vector<std::size_t> const& leaves) {
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
      },
      all);
  return result;
}

namespace {

/// appends "\n  <what>: <key>" to \p msg for each of \p keys
void list_keys(std::string& msg, std::string const& what,
               std::vector<std::string> const& keys) {
  for (auto const& k : keys) msg += "\n  " + what + ": " + k;
}

/// sorted, without repeats
std::vector<std::string> sorted_unique(std::vector<std::string> keys) {
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  return keys;
}

}  // namespace

PinSession::PinSession(eval::schedule::RootRecord const& root)
    : root_(root), uses_(root.terms.size(), 0) {
  for (std::size_t i = 0; i < root_.terms.size(); ++i)
    if (!position_.emplace(root_.terms[i].key, i).second)
      record_duplicates_.push_back(root_.terms[i].key);
}

TermRecord const* PinSession::take(std::string const& key) {
  auto it = position_.find(key);
  std::lock_guard<std::mutex> lock(mutex_);
  if (it == position_.end()) {
    missing_.push_back(key);
    return nullptr;
  }
  if (++uses_[it->second] > 1) duplicates_.push_back(key);
  return &root_.terms[it->second];
}

std::optional<std::size_t> PinSession::document_position(
    std::string const& key) const {
  auto it = position_.find(key);
  if (it == position_.end()) return std::nullopt;
  return it->second;
}

void PinSession::finish() const {
  std::lock_guard<std::mutex> lock(mutex_);
  std::vector<std::string> unused;
  for (std::size_t i = 0; i < root_.terms.size(); ++i)
    if (uses_[i] == 0 && position_.at(root_.terms[i].key) == i)
      unused.push_back(root_.terms[i].key);
  auto duplicates = duplicates_;
  duplicates.insert(duplicates.end(), record_duplicates_.begin(),
                    record_duplicates_.end());
  if (missing_.empty() && unused.empty() && duplicates.empty()) return;
  std::string msg = "schedule pin: root " + root_.label +
                    " does not match its record (no partial pin):";
  list_keys(msg, "missing", sorted_unique(missing_));
  list_keys(msg, "unused", unused);
  list_keys(msg, "duplicate", sorted_unique(std::move(duplicates)));
  throw Exception(msg);
}

void PinSession::bind(Expr const* optimized, std::string key) {
  std::lock_guard<std::mutex> lock(mutex_);
  // newest wins: an address freed by a per-summand caller may be reused by a
  // later result, and only the live result's key is meaningful
  keys_[optimized] = std::move(key);
}

std::optional<std::string> PinSession::key_of(Expr const* optimized) const {
  std::lock_guard<std::mutex> lock(mutex_);
  auto it = keys_.find(optimized);
  if (it == keys_.end()) return std::nullopt;
  return it->second;
}

RecordSession::RecordSession(std::string root_label)
    : label_(std::move(root_label)) {}

void RecordSession::add(Expr const* optimized, TermRecord rec) {
  std::lock_guard<std::mutex> lock(mutex_);
  // append unconditionally: a per-summand caller that drops each result frees
  // its address, which the next result may reuse; the record must survive
  // that. The pointer map follows the newest record for that address.
  index_[optimized] = records_.size();
  records_.push_back(std::move(rec));
}

void RecordSession::relocate(
    std::vector<std::pair<Expr const*, Expr const*>> const& moves) {
  std::lock_guard<std::mutex> lock(mutex_);
  // look all sources up first: a target may equal another move's source
  std::vector<std::pair<std::size_t, Expr const*>> found;
  for (auto const& [from, to] : moves) {
    auto it = index_.find(from);
    if (it == index_.end()) continue;
    found.emplace_back(it->second, to);
    index_.erase(it);
  }
  for (auto const& [i, to] : found) index_[to] = i;
}

namespace {

eval::schedule::RootRecord checked_root(std::string const& label,
                                        std::vector<TermRecord> terms) {
  std::vector<std::string> keys, duplicates;
  for (auto const& t : terms) keys.push_back(t.key);
  std::sort(keys.begin(), keys.end());
  for (std::size_t i = 1; i < keys.size(); ++i)
    if (keys[i] == keys[i - 1]) duplicates.push_back(keys[i]);
  if (!duplicates.empty()) {
    std::string msg = "schedule record: root " + label +
                      " has the same term more than once (simplify the Sum):";
    list_keys(msg, "duplicate", sorted_unique(std::move(duplicates)));
    throw Exception(msg);
  }
  return {label, std::move(terms)};
}

}  // namespace

eval::schedule::RootRecord RecordSession::finish(
    std::vector<Expr const*> const& order) const {
  std::lock_guard<std::mutex> lock(mutex_);
  std::vector<TermRecord> terms;
  terms.reserve(order.size());
  for (auto const* p : order) {
    auto it = index_.find(p);
    if (it == index_.end())
      throw Exception("schedule record: root " + label_ +
                      ": a summand of the final Sum has no record");
    terms.push_back(records_[it->second]);
  }
  return checked_root(label_, std::move(terms));
}

eval::schedule::RootRecord RecordSession::finish() const {
  std::lock_guard<std::mutex> lock(mutex_);
  std::vector<TermRecord> terms;
  terms.reserve(records_.size());
  for (auto const& rec : records_) terms.push_back(rec);
  return checked_root(label_, std::move(terms));
}

RecordedRoot record_root(std::string label, ExprPtr const& root_sum,
                         OptimizeOptions opts) {
  RecordSession session(std::move(label));
  opts.recording = &session;
  RecordedRoot result;
  result.optimized = optimize_result(root_sum, opts).expr;
  std::vector<Expr const*> order;
  if (result.optimized->is<Sum>())
    for (auto const& s : result.optimized->as<Sum>().summands())
      order.push_back(s.get());
  else
    order.push_back(result.optimized.get());
  result.record = session.finish(order);
  return result;
}

namespace detail {

namespace {

std::string objective_name(ObjectiveFunction f) {
  switch (f) {
    case ObjectiveFunction::DenseFLOPs:
      return "DenseFLOPs";
    case ObjectiveFunction::DenseSize:
      return "DenseSize";
    case ObjectiveFunction::DenseSpaceTime:
      return "DenseSpaceTime";
    case ObjectiveFunction::DenseTimeSpace:
      return "DenseTimeSpace";
    case ObjectiveFunction::DenseSpaceTimeBatched:
      return "DenseSpaceTimeBatched";
    case ObjectiveFunction::DenseTimeSpaceBatched:
      return "DenseTimeSpaceBatched";
  }
  throw Exception("fingerprint: unknown ObjectiveFunction");
}

std::string reorder_name(ReorderSum r) {
  return r == ReorderSum::Reorder ? "Reorder" : "NoReorder";
}

std::string scheduler_name(BatchScheduler s) {
  return s == BatchScheduler::ordered ? "ordered" : "forest_descent";
}

std::string vacuum_name(Vacuum v) {
  switch (v) {
    case Vacuum::Physical:
      return "Physical";
    case Vacuum::SingleProduct:
      return "SingleProduct";
    case Vacuum::MultiProduct:
      return "MultiProduct";
  }
  throw Exception("fingerprint: unknown Vacuum");
}

std::string spbasis_name(SPBasis b) {
  return b == SPBasis::Spinor ? "Spinor" : "Spinfree";
}

}  // namespace

eval::schedule::Fingerprint make_fingerprint(
    std::vector<std::string> root_labels, OptimizeOptions const& opts,
    std::vector<std::wstring> const& mode_order,
    std::map<std::wstring, Index> const& representatives,
    std::set<std::string> const& volatile_leaf_labels) {
  eval::schedule::Fingerprint fp;
  fp.roots = std::move(root_labels);
  auto const& policy = opts.batch_policy;
  for (auto const& [space, ix] : representatives) {
    fp.batchable[toUtf8(space)] = {policy.is_batchable_contracted_index &&
                                       policy.is_batchable_contracted_index(ix),
                                   policy.is_batchable_external_index &&
                                       policy.is_batchable_external_index(ix)};
    fp.context.spaces.push_back(toUtf8(space));
  }
  std::sort(fp.context.spaces.begin(), fp.context.spaces.end());
  fp.volatile_leaves.assign(volatile_leaf_labels.begin(),
                            volatile_leaf_labels.end());
  fp.persistent_only = policy.persistent_only;
  fp.batch_spectator_indices = policy.batch_spectator_indices;
  fp.accumulation_factor = policy.accumulation_factor;
  fp.scheduler = scheduler_name(policy.scheduler);
  fp.objective_function = objective_name(opts.objective_function);
  fp.reorder = reorder_name(opts.reorder);
  fp.cse["subnet"] = opts.CSE.subnet ? "true" : "false";
  for (auto const& m : mode_order) fp.mode_order.push_back(toUtf8(m));
  auto const& ctx = get_default_context();
  fp.context.vacuum = vacuum_name(ctx.vacuum());
  fp.context.spbasis = spbasis_name(ctx.spbasis());
  fp.context.strict_braket_symmetry = ctx.assert_strict_braket_symmetry();
  for (auto const& l : ctx.cardinal_tensor_labels())
    fp.context.cardinal_tensor_labels.push_back(toUtf8(l));
  return fp;
}

}  // namespace detail

namespace {

std::string join(std::vector<std::string> const& v) {
  std::string s = "[";
  for (std::size_t i = 0; i != v.size(); ++i) s += (i ? "," : "") + v[i];
  return s + "]";
}

std::string str(std::string const& s) { return s; }
std::string str(bool b) { return b ? "true" : "false"; }
std::string str(double d) { return std::to_string(d); }
std::string str(std::vector<std::string> const& v) { return join(v); }
std::string str(std::pair<bool, bool> const& p) {
  return "(" + str(p.first) + "," + str(p.second) + ")";
}
std::string str(std::map<std::string, std::string> const& m) {
  std::string s = "{";
  bool first = true;
  for (auto const& [k, v] : m) {
    s += (first ? "" : ",") + k + "=" + v;
    first = false;
  }
  return s + "}";
}

template <typename T>
void diff(std::vector<std::string>& out, std::string const& field, T const& a,
          T const& b) {
  if (!(a == b))
    out.push_back("fingerprint." + field + ": " + str(a) + " vs " + str(b));
}

}  // namespace

void check_fingerprint(eval::schedule::Fingerprint const& loaded,
                       eval::schedule::Fingerprint const& live) {
  std::vector<std::string> d;
  diff(d, "roots", loaded.roots, live.roots);
  // per-space predicate results; a space absent on one side is reported too
  std::set<std::string> spaces;
  for (auto const& [k, v] : loaded.batchable) spaces.insert(k);
  for (auto const& [k, v] : live.batchable) spaces.insert(k);
  for (auto const& sp : spaces) {
    auto a = loaded.batchable.find(sp);
    auto b = live.batchable.find(sp);
    if (a == loaded.batchable.end() || b == live.batchable.end())
      d.push_back("fingerprint.batchable." + sp + ": " +
                  (a == loaded.batchable.end() ? "absent" : str(a->second)) +
                  " vs " +
                  (b == live.batchable.end() ? "absent" : str(b->second)));
    else
      diff(d, "batchable." + sp, a->second, b->second);
  }
  diff(d, "volatile_leaves", loaded.volatile_leaves, live.volatile_leaves);
  diff(d, "persistent_only", loaded.persistent_only, live.persistent_only);
  diff(d, "batch_spectator_indices", loaded.batch_spectator_indices,
       live.batch_spectator_indices);
  diff(d, "accumulation_factor", loaded.accumulation_factor,
       live.accumulation_factor);
  diff(d, "scheduler", loaded.scheduler, live.scheduler);
  diff(d, "objective_function", loaded.objective_function,
       live.objective_function);
  diff(d, "reorder", loaded.reorder, live.reorder);
  diff(d, "cse", loaded.cse, live.cse);
  diff(d, "mode_order", loaded.mode_order, live.mode_order);
  auto const& a = loaded.context;
  auto const& b = live.context;
  diff(d, "context.vacuum", a.vacuum, b.vacuum);
  diff(d, "context.spbasis", a.spbasis, b.spbasis);
  diff(d, "context.strict_braket_symmetry", a.strict_braket_symmetry,
       b.strict_braket_symmetry);
  diff(d, "context.cardinal_tensor_labels", a.cardinal_tensor_labels,
       b.cardinal_tensor_labels);
  diff(d, "context.spaces", a.spaces, b.spaces);
  if (d.empty()) return;
  std::string msg = "schedule fingerprint mismatch:";
  for (auto const& line : d) msg += "\n" + line;
  throw Exception(msg);
}

std::map<std::string, std::string> planner_record(OptimizeOptions const& opts) {
  auto set = [](auto const& f) -> std::string { return f ? "set" : "unset"; };
  std::map<std::string, std::string> r;
  r["idx_to_extent"] = set(opts.idx_to_extent);
  r["inner_pow"] = set(opts.inner_pow);
  r["batch_target_size"] = set(opts.batch_policy.batch_target_size);
  r["peak_threshold"] = std::to_string(opts.batch_policy.peak_threshold);
  r["volatile_weight"] = std::to_string(opts.volatile_weight);
  r["footprint_weight"] = std::to_string(opts.footprint_weight);
  r["peak_flops_tolerance"] = std::to_string(opts.peak_flops_tolerance);
  r["prune_outer_products"] = opts.prune_outer_products ? "true" : "false";
  r["roofline.machine_balance"] = std::to_string(opts.roofline.machine_balance);
  r["roofline.fast_mem_elems"] = std::to_string(opts.roofline.fast_mem_elems);
  r["roofline.block_tiles"] = std::to_string(opts.roofline.block_tiles);
  r["roofline.block_prefactor"] = std::to_string(opts.roofline.block_prefactor);
  return r;
}

}  // namespace sequant::opt::schedule
