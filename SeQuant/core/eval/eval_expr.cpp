#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/functional/not_fn.hpp>
#include <range/v3/range/operations.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/move.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <cmath>
#include <ranges>
#include <string_view>
#include <type_traits>
#include <utility>

namespace sequant {

using EvalExprNode = FullBinaryNode<EvalExpr>;

namespace {

size_t hash_terminal_tensor(Tensor const&) noexcept;

bool is_tot(Tensor const& t) noexcept {
  return ranges::any_of(t.const_indices(), &Index::has_proto_indices);
}

// Slot-derived metadata -- an intermediate's bra/ket partition, which fixes
// its result-column grouping -- must be computed on the unfolded spelling;
// see sequant::value_oriented (core/expressions/tensor.hpp).

}  // namespace

namespace detail {
inline constexpr std::wstring_view label_tensor{L"I"};
inline constexpr std::wstring_view label_scalar{L"Z"};

template <std::ranges::range Bra, std::ranges::range Ket,
          std::ranges::range Aux>
ExprPtr make_tensor(const BinarizationOptions& opts, bra<Bra> b, ket<Ket> k,
                    aux<Aux> a, Symmetry symm, BraKetSymmetry bksymm,
                    ColumnSymmetry csymm) {
  // This function is creating intermediate tensors, which don't come with
  // an externally provided "correct"/canonical order of its indices.
  // Hence, we are free to define our own canonical order, which we
  // conveniently set to the indices being sorted in each group.
  if (opts.merge_indices) {
    using std::ranges::begin;
    using std::ranges::end;

    Index::index_vector indices;
    indices.insert(indices.end(), begin(b), end(b));
    indices.insert(indices.end(), begin(k), end(k));
    indices.insert(indices.end(), begin(a), end(a));

    std::ranges::sort(indices);
    return ex<Tensor>(label_tensor, bra(), ket(), aux(std::move(indices)), symm,
                      bksymm, csymm);
  } else {
    std::ranges::sort(b);
    std::ranges::sort(k);
    std::ranges::sort(a);

    return ex<Tensor>(label_tensor, std::move(b), std::move(k), std::move(a),
                      symm, bksymm, csymm);
  }
}

template <std::ranges::range Bra, std::ranges::range Ket,
          std::ranges::range Aux>
ExprPtr make_tensor_wo_symmetries(const BinarizationOptions& opts, bra<Bra>&& b,
                                  ket<Ket>&& k, aux<Aux>&& a) {
  return make_tensor<Bra, Ket, Aux>(opts, b, k, a, Symmetry::Nonsymm,
                                    BraKetSymmetry::Nonsymm,
                                    ColumnSymmetry::Nonsymm);
}

ExprPtr make_tensor(Tensor const& t, bool with_symm,
                    const BinarizationOptions& opts) {
  Symmetry symm = with_symm ? t.symmetry() : Symmetry::Nonsymm;
  BraKetSymmetry bksymm =
      with_symm ? t.braket_symmetry() : BraKetSymmetry::Nonsymm;
  ColumnSymmetry csymm =
      with_symm ? t.column_symmetry() : ColumnSymmetry::Nonsymm;

  return make_tensor(opts, bra(t.bra()), ket(t.ket()), aux(t.aux()), symm,
                     bksymm, csymm);
}

ExprPtr make_variable() { return ex<Variable>(label_scalar); }

}  // namespace detail

std::string to_label_annotation(const Index& idx) {
  using namespace ranges::views;
  using ranges::to;

  return toUtf8(idx.label()) +
         (idx.proto_indices() | transform(&Index::label) |
          transform([](auto&& str) { return toUtf8(str); }) |
          ranges::views::join | to<std::string>);
}

std::string EvalExpr::indices_annot() const noexcept {
  using ranges::views::filter;
  using ranges::views::join;
  using ranges::views::transform;

  if (!is_tensor()) return {};
  auto outer = csv_labels(canon_indices_  //
                          | filter(ranges::not_fn(&Index::has_proto_indices)));

  auto inner = csv_labels(canon_indices_  //
                          | filter(&Index::has_proto_indices));

  return outer + (inner.empty() ? "" : (";" + inner));
}

EvalExpr::index_vector const& EvalExpr::canon_indices() const noexcept {
  return canon_indices_;
}

EvalExpr::EvalExpr(Tensor const& tnsr)
    : op_type_{std::nullopt},
      result_type_{ResultType::Tensor},
      expr_{tnsr.clone()} {
  SEQUANT_ASSERT(!tnsr.indices().empty());
  if (is_tot(tnsr)) {
    ExprPtrList tlist{expr_};
    auto tn = TensorNetwork(tlist);
    // N.B. pass default_idxptr_slottype_lesscompare{} explicitly, NOT {}: an
    // empty named_index_compare selects canonicalize_slots' internal fallback,
    // which orders named indices by space() ALONE, whereas the declared default
    // (default_idxptr_slottype_lesscompare) orders by proto-index count first.
    // The latter is what makes a proto-indexed (ToT) leaf's canon_indices
    // put occupieds first (canonicals.hpp) -- a layout downstream
    // coefficient-shape detectors rely on. Passing {} here silently broke
    // that, so name
    // it explicitly.
    auto md =
        tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels(),
                              nullptr, default_idxptr_slottype_lesscompare{});
    hash_value_ = md.hash_value();
    canon_phase_ = md.phase;
    // The graph hash is orientation-shared (bra/ket of a Conjugate tensor are
    // colored identically), so both orientations land on one cache slot. When
    // the canonical orientation is the swapped one (md.conj), rewrite expr_
    // to the canonical spelling: swap (the Conjugate adjoint) + toggle the
    // elementwise-conjugation marker; binarize(Tensor) serves the marker via
    // an EvalOp::Adjoint wrapper over the shared operand.
    if (md.conj) {
      auto& tt = expr_->as<Tensor>();
      tt.adjoint();
      tt.conjugate();
    }
    canon_indices_ = md.get_indices<index_vector>();
    connectivity_ = std::move(md.graph);
  } else {
    // Single (protoindex-free) tensor: block-canonicalize it in place. This is
    // a lightweight per-tensor canonicalization (no deep tensor-network
    // canonicalization is needed for a tensor that is not itself a network),
    // and it normalizes bra<->ket orientation for braket-symmetric tensors so
    // that equivalent half-tensor forms (e.g. X{a;;x} and X{;a;x}) fold.
    auto& t = expr_->as<Tensor>();
    // apply() folds the two bra<->ket orientations of a flat Conjugate
    // tensor onto the canonical one, toggling the tensor's
    // elementwise-conjugation marker when it swaps (fold_conjugate_braket);
    // the hash below is that of the unconjugated spelling so both
    // orientations share a cache slot, and binarize(Tensor) serves the
    // marker via an EvalOp::Adjoint on retrieval.
    auto phase = TensorBlockCanonicalizer{}.apply(t);
    canon_phase_ = phase ? -1 : 1;
    // Leaf-hash invariant: the hash is always that of the UNSTARRED spelling,
    // so the two orientations of a Conjugate tensor share one cache slot; the
    // conjugation marker stays on expr_ (its symbolic spelling) and is served
    // by binarize's Adjoint wrapper on retrieval.
    if (t.conjugated()) {
      Tensor bare{t};
      bare.conjugate();
      hash_value_ = hash_terminal_tensor(bare);
    } else {
      hash_value_ = hash_terminal_tensor(t);
    }
    canon_indices_ = t.const_indices() | ranges::to<index_vector>;
  }
}

EvalExpr::EvalExpr(Constant const& c)
    : op_type_{std::nullopt},
      result_type_{ResultType::Scalar},
      expr_{c.clone()},
      hash_value_{hash::value(c)} {}

EvalExpr::EvalExpr(Variable const& v)
    : op_type_{std::nullopt},
      result_type_{ResultType::Scalar},
      expr_{v.clone()},
      hash_value_{hash::value(v)} {}

EvalExpr::EvalExpr(Power const& p)
    : op_type_{std::nullopt},
      result_type_{ResultType::Scalar},
      expr_{p.clone()},
      hash_value_{hash::value(p)} {}

EvalExpr::EvalExpr(EvalOp op, ResultType res, ExprPtr const& ex,
                   index_vector ixs, std::int8_t p, size_t h,
                   std::shared_ptr<bliss::Graph> connectivity)
    : op_type_{op},
      result_type_{res},
      expr_{ex.clone()},
      canon_indices_{std::move(ixs)},
      canon_phase_{p},
      hash_value_{h},
      connectivity_{std::move(connectivity)} {
  if (connectivity_ != nullptr) {
    // Note: The non-const cmp function performs some internal cleanup that the
    // comparison depends on. However, we want to be able to do const
    // comparisons and hence we have to assume fully cleaned-up graphs which we
    // achieve by causing a self-cleanup of the graph via the non-const cmp
    // function.
    connectivity_->cmp(*connectivity_);
  }

  // Using Tensor objects to represent scalar results is just confusing
  SEQUANT_ASSERT(ex->is<Tensor>() == (res == ResultType::Tensor));
}

const std::optional<EvalOp>& EvalExpr::op_type() const noexcept {
  return op_type_;
}

ResultType EvalExpr::result_type() const noexcept { return result_type_; }

size_t EvalExpr::hash_value() const noexcept {
  // canon_phase (+1/-1) is part of the node's *value* identity: two nodes that
  // share a canonical graph/leaf but differ in antisymmetric-reorder parity
  // evaluate to negatives of each other (+T vs -T), so they must not share a
  // CSE cache slot. Folding it in here (rather than special-casing the
  // comparator) makes every node hash carry the phase. It is a no-op for real
  // closed-shell paths (every phase is +1, so all hashes shift uniformly and
  // the equality structure is unchanged); complex/Kramers paths, which do
  // produce -1 phases, are thereby kept apart.
  auto h = hash_value_;
  hash::combine(h, canon_phase_);
  return h;
}

ExprPtr EvalExpr::expr() const noexcept { return expr_; }

bool EvalExpr::tot() const noexcept {
  return ranges::any_of(canon_indices(), &Index::has_proto_indices);
}

std::wstring EvalExpr::to_latex() const noexcept { return expr_->to_latex(); }

Expr::type_id_type EvalExpr::type_id() const noexcept {
  return expr_->type_id();
}

bool EvalExpr::is_tensor() const noexcept {
  return expr().is<Tensor>() && result_type() == ResultType::Tensor;
}

bool EvalExpr::is_scalar() const noexcept { return !is_tensor(); }

bool EvalExpr::is_constant() const noexcept {
  return expr().is<Constant>() && result_type() == ResultType::Scalar;
}

bool EvalExpr::is_variable() const noexcept {
  return expr().is<Variable>() && result_type() == ResultType::Scalar;
}

bool EvalExpr::is_power() const noexcept {
  return expr().is<Power>() && result_type() == ResultType::Scalar;
}

bool EvalExpr::is_primary() const noexcept { return !op_type(); }

bool EvalExpr::is_sum() const noexcept { return op_type() == EvalOp::Sum; }

bool EvalExpr::is_product() const noexcept {
  return op_type() == EvalOp::Product;
}

bool EvalExpr::is_adjoint() const noexcept {
  return op_type() == EvalOp::Adjoint;
}

Tensor const& EvalExpr::as_tensor() const { return expr().as<Tensor>(); }

Constant const& EvalExpr::as_constant() const { return expr().as<Constant>(); }

Variable const& EvalExpr::as_variable() const { return expr().as<Variable>(); }

Power const& EvalExpr::as_power() const { return expr().as<Power>(); }

std::string EvalExpr::label() const noexcept {
  if (is_tensor())
    return toUtf8(as_tensor().label()) + "(" + indices_annot() + ")";
  else if (is_constant()) {
    return toUtf8(io::serialization::to_string(as_constant()));
  } else if (is_power()) {
    return toUtf8(io::serialization::to_string(as_power()));
  } else if (is_variable()) {
    return toUtf8(as_variable().label());
  } else {
    SEQUANT_ABORT("EvalExpr::label: unhandled expression type");
  }
}

std::int8_t EvalExpr::canon_phase() const noexcept { return canon_phase_; }

bool EvalExpr::has_connectivity_graph() const noexcept {
  return connectivity_ != nullptr;
}

const bliss::Graph& EvalExpr::connectivity_graph() const noexcept {
  SEQUANT_ASSERT(connectivity_ != nullptr);
  return *connectivity_;
}

std::shared_ptr<bliss::Graph> EvalExpr::copy_connectivity_graph()
    const noexcept {
  return connectivity_;
}

namespace {

///
/// \param bk iterable of sequant Index
/// \return combined hash values of the elements.
///
/// @note An Index object's IndexSpace type and quantum numbers contribute to
///       the hash.
///
template <typename T>
size_t hash_indices(T const& indices) noexcept {
  size_t h = 0;
  for (auto const& idx : indices) {
    hash::combine(h, hash::value(idx.space().type().to_int32()));
    hash::combine(h, hash::value(idx.space().qns().to_int32()));
    if (idx.has_proto_indices()) {
      hash::combine(h, hash::value(idx.proto_indices().size()));
      for (auto&& i : idx.proto_indices())
        hash::combine(h, hash::value(i.label()));
    }
  }
  return h;
}

size_t hash_terminal_tensor(Tensor const& tnsr) noexcept {
  size_t h = 0;
  hash::combine(h, hash::value(tnsr.label()));
  hash::combine(h, hash_indices(tnsr.const_slots()));
  return h;
}
}  // namespace

///
/// \brief Calls canon_hash on all inits subranges.
/// \see inits
/// \see canon_hash
///
template <typename Rng>
auto imed_hashes(Rng const& rng) {
  using std::views::transform;
  return inits(rng) | transform([](auto&& v) {
           return hash::range_unordered(ranges::begin(v), ranges::end(v));
         });
}

struct ExprWithHash {
  ExprPtr expr;
  size_t hash;
};

void all_indices(IndexSet& result, ExprPtr const& expr) {
  if (!expr) return;
  if (expr->is<Tensor>())
    for (auto&& ix : expr->as<Tensor>().const_indices()) result.emplace(ix);
  else if (expr->is<Sum>() && !expr->empty())
    all_indices(result, expr->front());
  else if (expr->is<Product>())
    for (auto&& fac : *expr) all_indices(result, fac);
}

IndexSet all_indices(ExprPtr const& expr) {
  IndexSet result;
  all_indices(result, expr);
  return result;
}

///
/// \brief Collect tensors appearing as a factor at the leaf node of a product
///        sub-tree, or, at the root node of a sum sub-tree.
///
template <typename Rng>
void collect_tensor_factors(EvalExprNode const& node,  //
                            Rng& collect) {
  static_assert(std::is_same_v<ranges::range_value_t<Rng>, ExprWithHash>);

  if (auto op = node->op_type();
      node->is_tensor() &&
      (!op || *op == EvalOp::Sum || *op == EvalOp::Adjoint))
    // Treat Adjoint the same as Sum here: it produces a tensor result that
    // enters a parent Product as a single factor — the parent shouldn't
    // try to recurse past the Adjoint boundary, just collect the adjointed
    // tensor (held in node->expr()) and move on.
    collect.emplace_back(
        ExprWithHash{.expr = node->expr(), .hash = node->hash_value()});
  else if (node->op_type() == EvalOp::Product && !node.leaf()) {
    collect_tensor_factors(node.left(), collect);
    collect_tensor_factors(node.right(), collect);
  }
}

EvalExprNode binarize(Constant const& c) { return EvalExprNode{EvalExpr{c}}; }

EvalExprNode binarize(Variable const& v) { return EvalExprNode{EvalExpr{v}}; }

EvalExprNode binarize(Power const& p) { return EvalExprNode{EvalExpr{p}}; }

namespace {
// Assemble the Adjoint(bare_leaf, Constant{1}) IR over the tensor `orig`,
// carrying the given slot order and phase. The right child is a sentinel
// (FullBinaryNode invariant; evaluate ignores it for EvalOp::Adjoint).
// Wrapper hash = bare-leaf hash ⊕ EvalOp::Adjoint, so the wrapped
// orientation gets its own cache slot layered over the shared operand.
// Shared by the '⁺'-marked-adjoint and Conjugate-fold paths of
// binarize(Tensor).
EvalExprNode make_adjoint_over(Tensor const& orig, EvalExprNode bare_leaf,
                               EvalExpr::index_vector idxs, std::int8_t phase) {
  EvalExprNode sentinel{EvalExpr{Constant{1}}};
  auto h = bare_leaf->hash_value();
  hash::combine(h, static_cast<size_t>(EvalOp::Adjoint));
  EvalExpr adj{EvalOp::Adjoint,     //
               ResultType::Tensor,  //
               orig.clone(),        //
               std::move(idxs),     //
               phase,               //
               h,                   //
               nullptr};
  return EvalExprNode{std::move(adj), std::move(bare_leaf),
                      std::move(sentinel)};
}
}  // namespace

EvalExprNode binarize(Tensor const& t,
                      [[maybe_unused]] const BinarizationOptions& opts) {
  // Detect adjoint-marked tensor leaves (label ending in U+207A '⁺'). These
  // arise when the user wrote an adjoint of a BraKetSymmetry::Nonsymm tensor,
  // see Tensor::adjoint() in expressions/tensor.cpp. We surface the adjoint
  // as an explicit IR op (EvalOp::Adjoint) wrapping the bare-label operand,
  // so backends can serve T† by conjugating + permuting the cached T result.
  //
  // IR shape: Adjoint(Tensor{<bare>}, Constant{1})
  // The Constant(1) right child is a sentinel — present so the FullBinaryNode
  // invariant ("every non-leaf has two children") holds; evaluate ignores it
  // for EvalOp::Adjoint dispatch.
  if (!t.label().empty() && t.label().back() == adjoint_label) {
    // The Adjoint node carries the *adjointed* tensor (so its canon_indices
    // reflect the slot order parents see).

    // Build the bare-label operand: copy and call adjoint() to toggle the
    // marker off and swap bra/ket back to natural orientation.
    Tensor bare{t};
    bare.adjoint();
    SEQUANT_ASSERT(bare.label().empty() ||
                   bare.label().back() != adjoint_label);
    return make_adjoint_over(t, EvalExprNode{EvalExpr{bare}},
                             t.indices() | ranges::to<EvalExpr::index_vector>,
                             1);
  }

  // A leaf whose canonical spelling carries the elementwise-conjugation
  // marker (a BraKetSymmetry::Conjugate tensor authored in the swapped
  // orientation) is served via an EvalOp::Adjoint wrapper over the bare
  // (unconjugated) leaf, which holds the shared cached value. Unlike the
  // '⁺' case above (an explicit Nonsymm adjoint = conjugate *and*
  // transpose), the fold already put both orientations on the same canonical
  // slot order, so the wrapper carries that *same* order as its operand: the
  // adjoint() eval degenerates to a pure elementwise conjugation
  // (result(post) = operand(pre).conj() with post == pre, no permutation).
  EvalExpr leaf{t};
  if (leaf.expr()->is<Tensor>() && leaf.expr()->as<Tensor>().conjugated()) {
    // unstar the canonical spelling: that IS the bare operand (the fold
    // already put the slots in canonical orientation)
    Tensor bare{leaf.expr()->as<Tensor>()};
    bare.conjugate();
    EvalExprNode bare_leaf{EvalExpr{bare}};
    SEQUANT_ASSERT(!bare_leaf->expr()->as<Tensor>().conjugated());
    auto idxs = bare_leaf->canon_indices();
    auto phase = bare_leaf->canon_phase();
    return make_adjoint_over(leaf.expr()->as<Tensor>(), std::move(bare_leaf),
                             std::move(idxs), phase);
  }
  return EvalExprNode{std::move(leaf)};
}

EvalExprNode binarize(Sum const& sum, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter) {
  using ranges::views::move;
  using ranges::views::transform;
  auto summands =
      sum.summands()  //
      | transform([&uncontract, &opts, &node_counter](ExprPtr const& x) {
          return impl::binarize(x, uncontract, opts, node_counter);
        })  //
      | ranges::to_vector;

  bool const all_tensors =
      ranges::all_of(summands, [](auto&& n) { return n->is_tensor(); });

  [[maybe_unused]] bool const all_scalars =
      ranges::all_of(summands, [](auto&& n) { return n->is_scalar(); });

  SEQUANT_ASSERT(all_tensors | all_scalars);

  auto hvals = summands | transform([](auto&& n) { return n->hash_value(); });

  auto make_sum = [i = 0,                    //
                   hs = imed_hashes(hvals),  //
                   all_tensors, &opts](EvalExpr const& left,
                                       EvalExpr const&) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    if (all_tensors) {
      auto const t = value_oriented(left.as_tensor());
      return {
          EvalOp::Sum,         //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(opts, bra(t.bra()), ket(t.ket()),
                                            aux(t.aux())),  //
          left.canon_indices(),                             //
          1,                                                //
          h,                                                //
          nullptr};
    } else {
      return {EvalOp::Sum,              //
              ResultType::Scalar,       //
              detail::make_variable(),  //
              {},                       //
              1,                        //
              h,                        //
              nullptr};
    }
  };

  return fold_left_to_node(summands | move, make_sum);
}

EvalExprNode binarize(Product const& prod, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter) {
  using ranges::views::filter;
  using ranges::views::move;
  using ranges::views::transform;

  if (prod.factors().empty()) {
    return binarize(Constant(prod.scalar()));
  }

  auto const ltr_uncontr_idxs = [&]() {
    auto factor_idxs = prod.factors() |
                       transform([](auto&& xpr) { return all_indices(xpr); }) |
                       ranges::to_vector;
    return left_to_right_binarization_indices<Index, IndexSet>(factor_idxs,
                                                               uncontract);
  }();

  auto factors = prod.factors()  //
                 | transform([i = 0, &ltr_uncontr_idxs, &opts,
                              &node_counter](ExprPtr const& x) mutable {
                     return impl::binarize(x, ltr_uncontr_idxs.children[i++],
                                           opts, node_counter);
                   })  //
                 | ranges::to_vector;

  auto hvals = factors | transform([](auto&& n) { return n->hash_value(); });
  auto const hs = imed_hashes(hvals) | ranges::to_vector;

  auto make_prod = [i = 0, &hs, &ltr_uncontr_idxs, &opts, &node_counter](
                       EvalExprNode const& left,
                       EvalExprNode const& right) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    auto const& uncontracted_idxs = ltr_uncontr_idxs.imed[i];
    if (left->is_scalar() && right->is_scalar()) {
      // scalar * scalar
      return {EvalOp::Product,
              ResultType::Scalar,
              detail::make_variable(),
              {},
              1,
              h,
              nullptr};
    } else if (left->is_scalar() || right->is_scalar()) {
      // scalar * tensor or tensor * scalar
      auto const& tl = left->is_tensor() ? left : right;
      auto const t = value_oriented(tl->as_tensor());
      return {
          EvalOp::Product,     //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(opts, bra(t.bra()), ket(t.ket()),
                                            aux(t.aux())),  //
          tl->canon_indices(),                              //
          tl->canon_phase(),                                //
          h,
          nullptr};
    } else {
      // tensor * tensor
      container::svector<ExprWithHash> subfacs;
      collect_tensor_factors(left, subfacs);
      collect_tensor_factors(right, subfacs);
      auto ts = subfacs | transform([](auto&& t) { return t.expr; });
      IndexGroups<IndexVec> const target_indices = [&ts, &uncontracted_idxs]() {
        // route each surviving hyperindex to its correct slot
        // (bra, ket, or aux) based on which slot it occupies in
        // the factor tensors .. if appears in multiple slots put into aux
        //
        // count on the value orientation of each factor: a folded Conjugate
        // leaf is spelled swapped+starred but its indices occupy the authored
        // slots by value; counting the folded spelling would migrate its ket
        // group into bra and merge the intermediate's partition
        auto unfolded = ts | transform([](ExprPtr const& x) -> ExprPtr {
                          if (x->is<Tensor>() && x->as<Tensor>().conjugated())
                            return ex<Tensor>(value_oriented(x->as<Tensor>()));
                          return x;
                        }) |
                        ranges::to_vector;
        auto counts = get_used_indices_with_counts(ex<Product>(unfolded));
        IndexGroups<IndexVec> result;
        for (auto&& [k, v] : counts) {
          if (v.nonproto() == 0) continue;
          if (v.total() > 1) {
            if (uncontracted_idxs.contains(k)) result.aux.emplace_back(k);
            continue;
          }
          auto& group = v.bra ? result.bra : v.ket ? result.ket : result.aux;
          group.emplace_back(k);
        }
        return result;
      }();

      auto tn = TensorNetwork(ts);
      auto named_indices = tn.ext_indices();
      for (auto&& ix : uncontracted_idxs) named_indices.emplace(ix);

      auto canon = tn.canonicalize_slots(
          TensorCanonicalizer::cardinal_tensor_labels(), &named_indices);
      hash::combine(h, canon.hash_value());
      bool const scalar_result = canon.named_indices_canonical.empty();
      EvalExpr result =
          scalar_result
              ? EvalExpr{EvalOp::Product,          //
                         ResultType::Scalar,       //
                         detail::make_variable(),  //
                         {},                       //
                         canon.phase,              //
                         h,
                         std::move(canon.graph)}
              : EvalExpr{EvalOp::Product,     //
                         ResultType::Tensor,  //
                         detail::make_tensor_wo_symmetries(
                             opts, bra(target_indices.bra),
                             ket(target_indices.ket), aux(target_indices.aux)),
                         canon.get_indices<Index::index_vector>(),  //
                         canon.phase,                               //
                         h,
                         std::move(canon.graph)};
      // This is a genuine contraction (DP) node: the optimizer's
      // node_batch_axes carries one entry per such node, in the same
      // left-first post-order (children -- built by the recursive
      // impl::binarize calls above, which all run before this lambda is
      // invoked -- fully processed before this node). Stamp it if the caller
      // supplied per-node modes; always advance node_counter regardless, so
      // the top-level SEQUANT_ASSERT(node_counter ==
      // opts.node_batch_axes.size()) in binarize(ExprPtr, ...) can catch a
      // misaligned optimizer/binarize post-order.
      if (node_counter < opts.node_batch_axes.size()) {
        auto const& ann = opts.node_batch_axes[node_counter];
        result.set_batched_here(ann.axes);
        result.set_contracted_modes(ann.contracted_modes);
        result.set_batch_order_aware(ann.order_aware);
        result.set_batch_effective_count(ann.effective_count);
      }
      ++node_counter;
      return result;
    }
  };

  if (prod.scalar() == 1) {
    return fold_left_to_node(factors | move, make_prod);
  } else {
    auto left = fold_left_to_node(factors | move, make_prod);
    auto right = binarize(Constant{prod.scalar()});

    auto expr = left->is_tensor()
                    ? detail::make_tensor(value_oriented(left->as_tensor()),
                                          false, opts)
                : left->is_constant() ? (left->expr() * right->expr())
                                      : detail::make_variable();
    auto type = left->is_tensor() ? ResultType::Tensor : ResultType::Scalar;

    auto h = left->hash_value();
    hash::combine(h, right->hash_value());
    auto result = EvalExpr{EvalOp::Product,        //
                           type,                   //
                           expr,                   //
                           left->canon_indices(),  //
                           left->canon_phase(),    //
                           h,                      //
                           nullptr};

    return EvalExprNode{std::move(result), std::move(left), std::move(right)};
  }
}

namespace impl {

EvalExprNode binarize(ExprPtr const& expr, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter) {
  if (expr->is<Constant>())  //
    return binarize(expr->as<Constant>());

  if (expr->is<Variable>())  //
    return binarize(expr->as<Variable>());

  if (expr->is<Tensor>())  //
    return binarize(expr->as<Tensor>(), opts);

  if (expr->is<Sum>())  //
    return binarize(expr->as<Sum>(), uncontract, opts, node_counter);

  if (expr->is<Product>())  //
    return binarize(expr->as<Product>(), uncontract, opts, node_counter);

  if (expr->is<Power>())  //
    return binarize(expr->as<Power>());

  throw Exception("Encountered unsupported expression in binarize.");
}

}  // namespace impl

}  // namespace sequant
