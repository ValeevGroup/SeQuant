#include <SeQuant/core/optimize/multiterm.hpp>

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize/single_term.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sequant::opt {

namespace {

using Node = FullBinaryNode<EvalExpr>;
using index_vector = EvalExpr::index_vector;
using scalar_type = Constant::scalar_type;    ///< Complex<rational>
using Signature = std::vector<std::int64_t>;  ///< per-bucket signature key

// ===========================================================================
// Cost model
// ===========================================================================

/// Cost model for biclique scoring: owns the active OptimizeOptions and answers
/// the three size/cost questions the search asks. Constructed once per
/// factorize_multiterm call so the metric and extent map need not be threaded
/// through every search signature.
///
/// \note To be reconciled later with the existing cost counters in
/// \c single_term.hpp (\ref detail::flops_counter, \ref
/// detail::memsize_counter, \ref detail::footprint_counter) and the way
/// single-term optimization consumes them (compile-time \c ObjectiveFunction
/// dispatch, plus the volatile- and footprint-weighting in \ref
/// OptimizeOptions). \c contraction_cost already delegates to those counters;
/// \c tensor_size re-derives \c footprint_counter's element-count with a
/// different scalar/proto-index convention rather than calling it. We knowingly
/// keep it this way for now: a single cost model that every optimization pass
/// (single-term, multi-term, and whatever follows) shares is a longer-term
/// goal, and forcing premature unification here would either change single-term
/// behavior or freeze an interface before its requirements are clear. This
/// class is the first options-owning cost facade and the intended seam for that
/// consolidation.
class CostModel {
  OptimizeOptions const& opts_;

 public:
  explicit CostModel(OptimizeOptions const& o) : opts_(o) {}

  /// Product of index extents -- the number of elements of a tensor with the
  /// given result modes (1 for a scalar / empty index list).
  double tensor_size(index_vector const& idxs) const {
    double s = 1.;
    for (Index const& ix : idxs)
      s *= static_cast<double>(opts_.idx_to_extent(ix));
    return s;
  }

  /// Cost of the single binary contraction L*R -> result, under the active
  /// metric. The same counter that drives single-term optimization, so the two
  /// passes agree on "cheaper".
  double contraction_cost(index_vector const& l, index_vector const& r,
                          index_vector const& res) const {
    if (opts_.objective_function == ObjectiveFunction::DenseSize)
      return detail::memsize_counter(opts_.idx_to_extent)(l, r, res);
    return detail::flops_counter(opts_.idx_to_extent)(l, r, res);
  }

  /// Saving of an m x n biclique that rewrites m*n contractions L_i*R_j as one
  /// (sum_i L_i)*(sum_j R_j) plus the two factor-sums.
  ///
  /// \note The coefficient is (m*n - 1): the fold replaces m*n contractions
  /// with a single (sum L)*(sum R), so m*n - 1 contractions are avoided; from
  /// that we subtract the two factor-sum build costs. Result-add savings in the
  /// original sum are ignored, i.e. this is a conservative lower bound.
  double saving(std::size_t m, std::size_t n, double c_final, double l_size,
                double r_size) const {
    double const avoided = (static_cast<double>(m * n) - 1.) * c_final;
    double build = 0.;
    if (m > 1) build += static_cast<double>(m - 1) * l_size;
    if (n > 1) build += static_cast<double>(n - 1) * r_size;
    return avoided - build;
  }
};

// ===========================================================================
// Forest -> bipartite graph
// ===========================================================================

/// Label-exact, phase-relaxed factor equality. Two factor subtrees match iff
/// the structural/connectivity comparator agrees and their canonical result
/// indices are label-identical; a difference in canonical phase does not block
/// the match. Two factors that differ only by a canonicalization sign share one
/// vertex; the relative sign is folded into the edge coefficient (and thence
/// into a partner) at emission. canon_indices equality is retained: every
/// partner still exposes the same labeled free indices, so the factor-sums
/// remain summable with no relabeling (dummy-invariant matching is not handled
/// here).
struct FactorHash {
  std::size_t operator()(Node const* n) const {
    return TreeNodeHasher<Node>{}(*n);
  }
};
struct FactorEq {
  bool operator()(Node const* a, Node const* b) const {
    if (!TreeNodeEqualityComparator<Node>{}(*a, *b)) return false;
    index_vector const& ia = (*a)->canon_indices();
    index_vector const& ib = (*b)->canon_indices();
    if (ia.size() != ib.size()) return false;
    for (std::size_t k = 0; k < ia.size(); ++k)
      if (!(ia[k] == ib[k])) return false;
    return true;  // phase relaxed
  }
};

/// Interns factor subtrees into dense vertex ids 0..V-1 under the label-exact,
/// phase-relaxed predicate: a structurally identical factor gets
/// the same id whether it sits on the left or the right of a contraction and
/// regardless of its canonicalization sign. The first factor interned for a
/// vertex is its representative; \ref phase records that representative's
/// canonicalization phase so the per-edge relative sign can be recovered (see
/// \ref build_buckets). Ids are dense and assignment-ordered, so the reps and
/// phases are kept in flat vectors indexed by id.
class Interner {
  std::unordered_map<Node const*, int, FactorHash, FactorEq> vmap_;
  std::vector<Node const*> reps_;    ///< dense: id -> representative node
  std::vector<std::int8_t> phases_;  ///< dense: id -> representative phase

 public:
  /// Vertex id of \p n, assigning a fresh dense id (and recording \p n as the
  /// representative, with its canonicalization phase) on first sight.
  int intern(Node const& n) {
    auto [it, inserted] = vmap_.try_emplace(&n, static_cast<int>(vmap_.size()));
    if (inserted) {
      reps_.push_back(&n);
      phases_.push_back((*n).canon_phase());
    }
    return it->second;
  }
  Node const* rep(int id) const { return reps_[id]; }
  std::int8_t phase(int id) const { return phases_[id]; }
};

/// True if \p node is a top-level binary contraction of two tensor subtrees.
/// Only a tensor*tensor contraction carries a connectivity graph.
bool is_splittable(Node const& node) {
  return !node.leaf() && node->op_type() == EvalOp::Product &&
         node->has_connectivity_graph();
}

/// The contraction at the heart of a summand together with its scalar
/// prefactor. A plain tensor*tensor summand binarizes straight to a splittable
/// contraction (coeff == 1). A scalar-prefactored summand binarizes to
/// Product(contraction, Constant) with a null connectivity graph at the root;
/// we peel the Constant off and recurse to the contraction one level down so
/// that \c A*B - A*C and \c 2*A*B + 3*A*C participate in folding. Returns
/// \c std::nullopt for leaves, scalar*single-tensor, and
/// non-Constant (e.g. Variable) prefactors -- those pass through untouched.
struct PrefactoredContraction {
  Node const* node = nullptr;  ///< the splittable contraction
  scalar_type coeff{1};        ///< summand = coeff * to_expr(node)
};
std::optional<PrefactoredContraction> extract_core(Node const& root) {
  if (is_splittable(root)) return PrefactoredContraction{&root, scalar_type{1}};
  if (root.leaf() || root->op_type() != EvalOp::Product) return std::nullopt;
  // Product(contraction, scalar) -- scalar rides as a Constant child, the
  // contraction as the other (order is binarize's choice, so try both).
  auto try_pair =
      [](Node const& core,
         Node const& scal) -> std::optional<PrefactoredContraction> {
    if (scal->is_scalar() && scal->expr() && scal->expr()->is<Constant>() &&
        is_splittable(core))
      return PrefactoredContraction{&core,
                                    scal->expr()->as<Constant>().value()};
    return std::nullopt;
  };
  if (auto c = try_pair(root.left(), root.right())) return c;
  if (auto c = try_pair(root.right(), root.left())) return c;
  return std::nullopt;
}

/// Signature bucket key: the sorted multiset of index-space identifiers over
/// the final contraction's external and contracted indices (spaces only,
/// label-agnostic). Terms only factor together within a bucket; cross-bucket
/// bicliques are impossible by construction.
Signature signature(Node const& root) {
  index_vector const& ext = root->canon_indices();
  index_vector const& l = root.left()->canon_indices();

  auto push_space = [](Signature& v, Index const& ix) {
    v.push_back(static_cast<std::int64_t>(ix.space().attr()));
  };

  Signature key;
  for (Index const& ix : ext) push_space(key, ix);
  auto is_external = [&ext](Index const& ix) {
    for (Index const& e : ext)
      if (e == ix) return true;
    return false;
  };
  // Each contracted index appears on both L and R; scan L only to count once.
  for (Index const& ix : l)
    if (!is_external(ix)) push_space(key, ix);
  std::sort(key.begin(), key.end());
  return key;
}

/// One splittable summand: an edge between its left and right factor vertices,
/// carrying the per-summand coefficient that relates the summand's value to the
/// product of the two vertex representatives' expressions:
///   summand = coeff * to_expr(left_rep) * to_expr(right_rep).
/// coeff absorbs the scalar prefactor (\ref PrefactoredContraction) and the
/// relative canonicalization sign of each factor versus its vertex
/// representative.
struct Edge {
  std::size_t pos;    ///< position in sum.summands()
  int lvid, rvid;     ///< vertex ids of the left / right child factors
  scalar_type coeff;  ///< summand = coeff * L_rep * R_rep
};

struct Bucket {
  double c_final = 0.;  ///< memoized final-contraction cost (signature-only)
  std::vector<Edge> edges;
};

/// The live bipartite view of a bucket over not-yet-consumed edges. Parallel
/// edges (same vertex pair, e.g. duplicate summands) have their coefficients
/// summed and their positions concatenated. Built once per greedy round via
/// \ref build; the search reads it only through the accessors below, so the
/// underlying maps stay an implementation detail.
class Live {
  std::map<int, std::vector<int>> left_adj_;  ///< lvid -> sorted rvids
  std::map<int, Node const*> left_rep_, right_rep_;
  std::map<std::pair<int, int>, std::vector<std::size_t>> edge_pos_;
  std::map<std::pair<int, int>, scalar_type> edge_coeff_;

 public:
  static Live build(Bucket const& bucket, std::vector<bool> const& consumed,
                    Interner const& interner);

  /// lvid -> sorted rvids, over the live (not-yet-consumed) edges.
  std::map<int, std::vector<int>> const& left_adj() const { return left_adj_; }
  Node const* left_rep(int l) const { return left_rep_.at(l); }
  Node const* right_rep(int r) const { return right_rep_.at(r); }
  /// Summand positions covered by the (l, r) edge. Defined only for present
  /// edges -- throws on an absent vertex pair.
  std::vector<std::size_t> const& positions(int l, int r) const {
    return edge_pos_.at({l, r});
  }
  /// Folded coefficient of the (l, r) edge. Defined only for present edges.
  scalar_type coeff(int l, int r) const { return edge_coeff_.at({l, r}); }
};

Live Live::build(Bucket const& bucket, std::vector<bool> const& consumed,
                 Interner const& interner) {
  Live live;
  std::map<int, std::set<int>> adj;
  for (Edge const& e : bucket.edges) {
    if (consumed[e.pos]) continue;
    adj[e.lvid].insert(e.rvid);
    live.left_rep_.emplace(e.lvid, interner.rep(e.lvid));
    live.right_rep_.emplace(e.rvid, interner.rep(e.rvid));
    auto key = std::pair{e.lvid, e.rvid};
    live.edge_pos_[key].push_back(e.pos);
    auto [it, inserted] = live.edge_coeff_.try_emplace(key, e.coeff);
    if (!inserted) it->second = it->second + e.coeff;
  }
  for (auto const& [l, rs] : adj)
    live.left_adj_.emplace(l, std::vector<int>(rs.begin(), rs.end()));
  return live;
}

/// Build the per-signature buckets of factorable edges from the binarized
/// summands, interning factor vertices into \p interner along the way.
///
/// Each splittable summand contributes one edge between its left and right
/// factor vertices. The edge coefficient relates the summand's value to the
/// product of the two vertex representatives' expressions:
///   summand = coeff * to_expr(L_rep) * to_expr(R_rep),
///   coeff   = prefactor * sigma_L * sigma_R,
///   sigma   = canon_phase(factor) * canon_phase(representative)   (each +/-1),
/// so the relative canonicalization sign between a factor and its vertex
/// representative rides on the edge. Non-splittable
/// summands (leaf, scalar*leaf, opaque prefactor) contribute no edge and pass
/// through untouched. \c c_final is memoized per bucket on its first edge
/// (signature-only, hence equal for every member).
std::map<Signature, Bucket> build_buckets(container::vector<Node> const& nodes,
                                          Interner& interner,
                                          CostModel const& cost) {
  std::map<Signature, Bucket> buckets;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    auto core = extract_core(nodes[i]);
    if (!core) continue;  // leaf / scalar*leaf / opaque prefactor: untouched
    Node const& c = *core->node;
    Node const& lnode = c.left();
    Node const& rnode = c.right();
    int const lv = interner.intern(lnode), rv = interner.intern(rnode);
    scalar_type const sigma_l(static_cast<int>(lnode->canon_phase()) *
                              static_cast<int>(interner.phase(lv)));
    scalar_type const sigma_r(static_cast<int>(rnode->canon_phase()) *
                              static_cast<int>(interner.phase(rv)));
    scalar_type const coeff = core->coeff * sigma_l * sigma_r;

    auto& bucket = buckets[signature(c)];
    if (bucket.edges.empty())
      bucket.c_final = cost.contraction_cost(
          lnode->canon_indices(), rnode->canon_indices(), c->canon_indices());
    bucket.edges.push_back(Edge{i, lv, rv, coeff});
  }
  return buckets;
}

// ===========================================================================
// Biclique search
// ===========================================================================

std::vector<int> intersect_sorted(std::vector<int> const& a,
                                  std::vector<int> const& b) {
  std::vector<int> out;
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::back_inserter(out));
  return out;
}

/// A factorable biclique S_L x S_R within one bucket's live bipartite graph.
/// The coefficient of each covered summand factors as left_coeffs[i] *
/// right_coeffs[j] (rank-1), which is what lets it be emitted as the single
/// product (sum_i left_coeffs[i] L_i) * (sum_j right_coeffs[j] R_j).
struct Biclique {
  std::vector<int> left, right;        ///< vertex ids
  std::vector<std::size_t> positions;  ///< covered summand positions
  std::vector<Node const*> left_reps, right_reps;
  std::vector<scalar_type> left_coeffs, right_coeffs;
  double saving = 0.;
  std::size_t min_pos = 0;  ///< for deterministic tie-breaking
};

/// Strict "candidate beats incumbent" order: higher saving wins; ties broken by
/// lowest covered position (deterministic). A null incumbent is always beaten.
bool supersedes(Biclique const& cand, std::optional<Biclique> const& best) {
  return !best || cand.saving > best->saving ||
         (cand.saving == best->saving && cand.min_pos < best->min_pos);
}

/// Try to factor the m x n coefficient matrix c[i][j] = edge_coeff(left[i],
/// right[j]) as an outer product alpha[i] * beta[j] (multiplicative rank 1).
/// Returns nullopt if it does not factor (then the biclique must be reduced to
/// a one-sided fold, which is always factorable). For a complete biclique every
/// (left[i], right[j]) edge is present, so every lookup succeeds.
std::optional<std::pair<std::vector<scalar_type>, std::vector<scalar_type>>>
factor_coeffs(std::vector<int> const& left, std::vector<int> const& right,
              Live const& live) {
  auto at = [&](int l, int r) { return live.coeff(l, r); };
  scalar_type const c00 = at(left.front(), right.front());
  if (c00.is_zero()) return std::nullopt;
  // rank-1 iff c[i][j] * c[0][0] == c[i][0] * c[0][j] for all i, j.
  for (std::size_t i = 0; i < left.size(); ++i)
    for (std::size_t j = 0; j < right.size(); ++j)
      if (!(at(left[i], right[j]) * c00 ==
            at(left[i], right.front()) * at(left.front(), right[j])))
        return std::nullopt;
  std::vector<scalar_type> alpha, beta;
  alpha.reserve(left.size());
  beta.reserve(right.size());
  for (int l : left) alpha.push_back(at(l, right.front()));
  for (int r : right) beta.push_back(at(left.front(), r) / c00);
  return std::pair{std::move(alpha), std::move(beta)};
}

/// Fill the geometry (positions, reps, min_pos) of a biclique over the given
/// vertex sets, then score it. Coefficients (\p alpha, \p beta) are supplied by
/// the caller (already validated factorable).
///
/// \pre \p left x \p right is a complete biclique: every (l, r) pair is a live
/// edge. This walks the full product through \ref Live::positions, which throws
/// on an absent pair. Both callers honor it -- \ref best_biclique passes
/// maximal complete vertex sets, and \ref best_fold's one-sided slices are a
/// single row or column of one.
Biclique make_biclique(std::vector<int> left, std::vector<int> right,
                       std::vector<scalar_type> alpha,
                       std::vector<scalar_type> beta, double c_final,
                       Live const& live, CostModel const& cost) {
  Biclique bc;
  std::size_t const m = left.size(), n = right.size();
  double const l_size =
      cost.tensor_size((*live.left_rep(left.front()))->canon_indices());
  double const r_size =
      cost.tensor_size((*live.right_rep(right.front()))->canon_indices());
  bc.saving = cost.saving(m, n, c_final, l_size, r_size);

  std::size_t min_pos = std::numeric_limits<std::size_t>::max();
  for (int l : left)
    for (int r : right)
      for (std::size_t p : live.positions(l, r)) {
        bc.positions.push_back(p);
        min_pos = std::min(min_pos, p);
      }
  bc.min_pos = min_pos;
  for (int l : left) bc.left_reps.push_back(live.left_rep(l));
  for (int r : right) bc.right_reps.push_back(live.right_rep(r));
  bc.left = std::move(left);
  bc.right = std::move(right);
  bc.left_coeffs = std::move(alpha);
  bc.right_coeffs = std::move(beta);
  return bc;
}

/// The best factorable fold derivable from the complete biclique \p left x
/// \p right. If the coefficient matrix is rank-1 the full m x n fold is used;
/// otherwise it is reduced to the more profitable of a one-sided slice -- a
/// single left vertex against all of \p right, or all of \p left against a
/// single right vertex -- which is always factorable (the relative signs ride
/// entirely on the partner side). Returns nullopt if nothing is profitable.
std::optional<Biclique> best_fold(std::vector<int> const& left,
                                  std::vector<int> const& right, double c_final,
                                  Live const& live, CostModel const& cost) {
  std::optional<Biclique> best;
  auto consider = [&](Biclique bc) {
    if (bc.saving > 0. && supersedes(bc, best)) best = std::move(bc);
  };

  if (auto ab = factor_coeffs(left, right, live)) {
    consider(make_biclique(left, right, std::move(ab->first),
                           std::move(ab->second), c_final, live, cost));
  } else {
    // One-sided reductions (always rank-1). Coefficients fold onto the
    // many-side partner: the single side keeps coefficient 1.
    int const l0 = left.front();
    std::vector<scalar_type> beta;
    for (int r : right) beta.push_back(live.coeff(l0, r));
    consider(make_biclique({l0}, right, {scalar_type{1}}, std::move(beta),
                           c_final, live, cost));

    int const r0 = right.front();
    std::vector<scalar_type> alpha;
    for (int l : left) alpha.push_back(live.coeff(l, r0));
    consider(make_biclique(left, {r0}, std::move(alpha), {scalar_type{1}},
                           c_final, live, cost));
  }
  return best;
}

/// Internal safety backstop on the intersection-closure size -- deliberately a
/// file-local constant and not a field of \ref OptimizeOptions. It is not an
/// optimization tuning knob: the enumeration below is a first-cut, worst-case-
/// exponential closure (acceleration is deferred), and this merely caps its
/// memory/time on pathological inputs. Hitting it cannot change correctness:
/// the greedy driver re-enumerates every round, so a capped round only risks a
/// sub-optimal (never invalid) fold. There's no quality-for-speed trade a user
/// would sensibly dial. It sits far above any realistic problem size and so
/// never fires in practice; were it ever to fire, the remedy is the deferred
/// search acceleration, not a larger user-set bound.
constexpr std::size_t max_closure_size = 50000;

/// Enumerate maximal bicliques of \p live and return the highest-saving fold
/// with positive saving (deterministic tie-break: lowest covered position).
/// Returns nullopt if none is profitable. First cut: full enumeration via
/// intersection-closure of the left neighborhoods (worst-case exponential;
/// acceleration is deferred). \p cap bounds the
/// closure; if hit we stop growing it and keep the best biclique found so far
/// (still correct -- the greedy driver re-enumerates each round, so a bounded
/// round only risks sub-optimality, never an invalid fold).
std::optional<Biclique> best_biclique(Live const& live, double c_final,
                                      CostModel const& cost,
                                      std::size_t cap = max_closure_size) {
  // Closure of the left neighborhoods under intersection -> candidate right
  // vertex sets. Each closed set induces a maximal biclique.
  std::set<std::vector<int>> closed;
  std::vector<std::vector<int>> seeds;
  for (auto const& [l, nb] : live.left_adj())
    if (!nb.empty() && closed.insert(nb).second) seeds.push_back(nb);

  std::vector<std::vector<int>> queue = seeds;
  bool overflow = false;
  for (std::size_t qi = 0; qi < queue.size() && !overflow; ++qi) {
    for (std::vector<int> const& sd : seeds) {
      std::vector<int> inter = intersect_sorted(queue[qi], sd);
      if (!inter.empty() && closed.insert(inter).second) {
        queue.push_back(inter);
        if (closed.size() > cap) {
          overflow = true;
          break;
        }
      }
    }
  }

  std::optional<Biclique> best;
  for (std::vector<int> const& r0 : closed) {
    // S_L = every left vertex adjacent to all of r0.
    std::vector<int> left;
    for (auto const& [l, nb] : live.left_adj())
      if (std::includes(nb.begin(), nb.end(), r0.begin(), r0.end()))
        left.push_back(l);
    if (left.empty()) continue;
    // S_R = intersection of the neighborhoods of S_L (>= r0): makes it maximal.
    std::vector<int> right = live.left_adj().at(left.front());
    for (std::size_t k = 1; k < left.size(); ++k)
      right = intersect_sorted(right, live.left_adj().at(left[k]));
    if (right.empty()) continue;

    std::optional<Biclique> bc = best_fold(left, right, c_final, live, cost);
    if (!bc) continue;
    if (supersedes(*bc, best)) best = std::move(bc);
  }
  return best;
}

// ===========================================================================
// Emission
// ===========================================================================

/// (sum of) the given factor subtrees, each scaled by its coefficient, as a
/// symbolic ExprPtr. A unit coefficient contributes the bare factor; a
/// non-unit one wraps it in a scalar Product.
ExprPtr side_expr(std::vector<Node const*> const& reps,
                  std::vector<scalar_type> const& coeffs) {
  auto term = [](Node const* n, scalar_type c) -> ExprPtr {
    ExprPtr e = to_expr(*n);
    if (c.is_identity()) return e;
    return ex<Product>(c, ExprPtrList{e}, Product::Flatten::No);
  };
  if (reps.size() == 1) return term(reps.front(), coeffs.front());
  Sum::summands_type parts;
  parts.reserve(reps.size());
  for (std::size_t i = 0; i < reps.size(); ++i)
    parts.push_back(term(reps[i], coeffs[i]));
  return ex<Sum>(Sum(parts.begin(), parts.end()));
}

ExprPtr emit_biclique(Biclique const& bc) {
  return ex<Product>(ExprPtrList{side_expr(bc.left_reps, bc.left_coeffs),
                                 side_expr(bc.right_reps, bc.right_coeffs)},
                     Product::Flatten::No);
}

}  // namespace

ExprPtr factorize_multiterm(
    Sum const& sum, container::vector<FullBinaryNode<EvalExpr>> const& nodes,
    OptimizeOptions const& opts) {
  // Cost-driven selection needs index extents; guaranteed by construction in
  // the normal optimize() flow. No structural fallback.
  SEQUANT_ASSERT(opts.idx_to_extent);
  SEQUANT_ASSERT(nodes.size() == sum.size());

  CostModel cost{opts};
  std::size_t const N = sum.size();

  // Build: intern factor vertices and group splittable summands into
  // per-signature buckets of factorable edges.
  Interner interner;
  std::map<Signature, Bucket> buckets = build_buckets(nodes, interner, cost);

  std::vector<bool> consumed(N, false);
  std::vector<ExprPtr> folds;

  // Greedy cost-driven driver: repeatedly apply the single highest-saving
  // maximal biclique across all buckets, until none has positive saving.
  while (true) {
    std::optional<Biclique> best;
    for (auto const& [key, bucket] : buckets) {
      Live live = Live::build(bucket, consumed, interner);
      std::optional<Biclique> bc = best_biclique(live, bucket.c_final, cost);
      if (!bc) continue;
      if (supersedes(*bc, best)) best = std::move(bc);
    }
    if (!best) break;

    for (std::size_t p : best->positions) consumed[p] = true;
    folds.push_back(emit_biclique(*best));
  }

  if (folds.empty()) return ex<Sum>(sum);  // nothing profitable

  // Reassemble: untouched summands (original order) followed by folds.
  Sum::summands_type out;
  out.reserve(N + folds.size());
  for (std::size_t i = 0; i < N; ++i)
    if (!consumed[i]) out.push_back(sum.summand(i));
  for (ExprPtr const& f : folds) out.push_back(f);

  return ex<Sum>(Sum(out.begin(), out.end()));
}

}  // namespace sequant::opt
