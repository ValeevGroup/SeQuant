#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/complex.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/functional/not_fn.hpp>
#include <range/v3/range/operations.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/move.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ranges>
#include <stdexcept>
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
                    ColumnSymmetry csymm, bool keep_order = false) {
  // This function is creating intermediate tensors, which don't come with
  // an externally provided "correct"/canonical order of its indices.
  // Hence, we are free to define our own canonical order, which we
  // conveniently set to the indices being sorted in each group -- EXCEPT
  // when the caller passes an order that must be kept (keep_order): the
  // placeholder of an opaque node (a Sum) is what an enclosing tensor network
  // sees, so its slots must be spelled in the order the value is laid out in
  // (see binarize(Sum)).
  if (opts.merge_indices) {
    using std::ranges::begin;
    using std::ranges::end;

    Index::index_vector indices;
    indices.insert(indices.end(), begin(b), end(b));
    indices.insert(indices.end(), begin(k), end(k));
    indices.insert(indices.end(), begin(a), end(a));

    if (!keep_order) std::ranges::sort(indices);
    return ex<Tensor>(label_tensor, bra(), ket(), aux(std::move(indices)), symm,
                      bksymm, csymm);
  } else {
    if (!keep_order) {
      std::ranges::sort(b);
      std::ranges::sort(k);
      std::ranges::sort(a);
    }

    return ex<Tensor>(label_tensor, std::move(b), std::move(k), std::move(a),
                      symm, bksymm, csymm);
  }
}

template <std::ranges::range Bra, std::ranges::range Ket,
          std::ranges::range Aux>
ExprPtr make_tensor_wo_symmetries(const BinarizationOptions& opts, bra<Bra>&& b,
                                  ket<Ket>&& k, aux<Aux>&& a,
                                  bool keep_order = false) {
  return make_tensor<Bra, Ket, Aux>(opts, b, k, a, Symmetry::Nonsymm,
                                    BraKetSymmetry::Nonsymm,
                                    ColumnSymmetry::Nonsymm, keep_order);
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

namespace {
/// maps a folded leaf's canonical indices back to the as-written flavors
/// (the Kramers flip is an involution; ordinals and order are kept)
void kramers_flip_indices_as_written(EvalExpr::index_vector& ixs) {
  if (const auto isr = get_default_context().index_space_registry())
    kramers_flip(ixs, *isr);
}

/// the eval-leaf Kramers fold is an explicit context opt-in (see
/// CanonicalizeOptions::fold_kramers_eval_leaves)
bool fold_kramers_leaf() {
  return CanonicalizeOptions::default_options().fold_kramers_eval_leaves ==
         CanonicalizeOptions::FoldKramersEvalLeaves::Yes;
}

/// Normalizes a leaf tensor's SPELLING channels into transform bits:
/// strips a '⁺' adjoint label (adjoint = conj ∘ swap) and converts the
/// elementwise-conjugation marker to a PURE {conj} bit (slots untouched;
/// orientation deltas belong to the canonicalizer fold alone). Symm markers
/// are value-redundant and dropped. Returns the accumulated transform.
CanonTransform normalize_leaf_spelling(Tensor& t) {
  CanonTransform tr{};
  if (!t.label().empty() && t.label().back() == adjoint_label) {
    t.adjoint();  // removes the label, swaps slots back
    tr = compose(tr, {.conj = true, .braket_swap = true});
  }
  if (t.conjugated()) {
    if (t.braket_symmetry() != BraKetSymmetry::Symm)
      tr = compose(tr, {.conj = true});
    t.conjugate();  // the unmarked spelling is stored
  }
  return tr;
}

struct LeafNormalization {
  CanonTransform transform;    ///< maps the stored spelling to the as-written
  bool kramers_fired = false;  ///< the Kramers fold respelled the leaf
};

/// Respells a leaf tensor IN PLACE as its canonical block form and returns
/// the retrieval transform, composing the channels in this order (their
/// inverse, in reverse order, is EvalExpr::denoted_expr):
///   1. spelling channels (normalize_leaf_spelling): adjoint label,
///      conjugation marker -> {conj, swap} / {conj} bits;
///   2. the eval-leaf Kramers fold (explicit opt-in): a down-first leaf is
///      respelled as its up-first partner; the fold's marker is a pure
///      {conj} bit and its phase multiplies in;
///   3. block canonicalization WITH the braket fold (the block
///      canonicalizer's own Kramers fold stays off): the antisymmetric
///      reorder phase multiplies in; a braket-fold marker (the fold swapped
///      a Conjugate tensor INTO its canonical orientation) becomes
///      {conj, swap}, the canonical slots are kept.
/// The stored spelling is unmarked, up-row, block-canonical: what a leaf
/// provider serves.
LeafNormalization normalize_leaf(Tensor& t) {
  LeafNormalization result;
  auto& tr = result.transform;
  tr = compose(tr, normalize_leaf_spelling(t));
  const int kramers_phase = fold_kramers_leaf() ? canonicalize_kramers(t) : 1;
  result.kramers_fired = t.conjugated();
  if (result.kramers_fired) {
    tr = compose(tr, {.conj = true});
    t.conjugate();
  }
  const auto block_byproduct = TensorBlockCanonicalizer{}.apply(t);
  tr = compose(tr, {.phase = static_cast<std::int8_t>(
                        (block_byproduct ? -1 : 1) * kramers_phase)});
  if (t.conjugated()) {
    tr = compose(tr, {.conj = true, .braket_swap = true});
    t.conjugate();
  }
  return result;
}
}  // namespace

EvalExpr::EvalExpr(Tensor const& tnsr)
    : op_type_{std::nullopt},
      result_type_{ResultType::Tensor},
      expr_{tnsr.clone()} {
  SEQUANT_ASSERT(!tnsr.indices().empty());
  // the stored spelling is canonical (unmarked, up-row, block-canonical);
  // every spelling channel becomes a CanonTransform byproduct applied on
  // retrieval, so every route to one canonical spelling lands on one slot
  auto const [transform, kramers_fired] = normalize_leaf(expr_->as<Tensor>());
  canon_transform_ = transform;
  if (is_tot(tnsr)) {
    // slot identity: the canonical labeling of the block-canonical spelling.
    // The block canonicalizer is label-blind (same-space slots keep their
    // order), so the labeling finishes the reorder: every (anti)symmetric
    // bra/ket bundle is put into the labeling's canonical slot order (the
    // order its phase is defined against) with the permutation parity as
    // the phase -- so the two spellings t{a2,a3;..} and t{a3,a2;..} store
    // ONE spelling, share one slot, and differ by the retrieval phase only.
    // (The network clones its tensors: adopt the respelled one.)
    ExprPtrList tlist{expr_};
    auto tn = TensorNetwork(tlist);
    auto md = tn.canonicalize_slots(
        {.cardinal_tensor_labels =
             TensorCanonicalizer::cardinal_tensor_labels(),
         .apply_slot_order = true});
    hash_value_ = md.hash_value();
    canon_transform_ = compose(canon_transform_, {.phase = md.phase});
    expr_ = std::dynamic_pointer_cast<Expr>(tn.tensors().front());
    SEQUANT_ASSERT(expr_ && expr_->is<Tensor>());
    // array-faithful indices in the Nested (outer;inner) convention: a ToT
    // array's outer modes are the plain slots PLUS the proto constituents,
    // deterministically ordered by NestedTensorIndices (the md list is the
    // same set in named-canonical order, which annots must not depend on)
    auto const slot_ixs =
        expr_->as<Tensor>().const_indices() | ranges::to<index_vector>;
    auto const nti = tot_indices<index_vector>(slot_ixs);
    canon_indices_ =
        ranges::views::concat(nti.outer, nti.inner) | ranges::to<index_vector>;
    connectivity_ = std::move(md.graph);
  } else {
    auto const& t = expr_->as<Tensor>();
    hash_value_ = hash_terminal_tensor(t);
    canon_indices_ = t.const_indices() | ranges::to<index_vector>;
  }
  // T19 layer 2 contract: expr() keeps the FOLDED (up-row) spelling -- what
  // a leaf provider fetches -- while canon_indices() (the parent's
  // contraction labels; TA matches annotations, not spellings) carries the
  // as-written flavors in the same canonical order, so the served up block
  // + {conj, phase} denotes the as-written value
  if (kramers_fired) kramers_flip_indices_as_written(canon_indices_);
  kramers_folded_ = kramers_fired;
}

EvalExpr::EvalExpr(Constant const& c)
    : op_type_{std::nullopt},
      result_type_{ResultType::Scalar},
      expr_{c.clone()},
      hash_value_{hash::value(c)} {}

EvalExpr::EvalExpr(Variable const& v)
    : op_type_{std::nullopt},
      result_type_{ResultType::Scalar},
      expr_{v.clone()} {
  // a conjugation marker rides the transform; the unmarked spelling is
  // stored and hashed (one cache slot for x and x^*)
  auto& vv = expr_->as<Variable>();
  if (vv.conjugated()) {
    canon_transform_.conj = true;
    vv.conjugate();
  }
  hash_value_ = hash::value(vv);
}

EvalExpr::EvalExpr(Power const& p)
    : op_type_{std::nullopt},
      result_type_{ResultType::Scalar},
      expr_{p.clone()} {
  // conj(b^n) = conj(b)^n for integer n: the marker rides the transform
  auto& pp = expr_->as<Power>();
  if (pp.conjugated()) {
    canon_transform_.conj = true;
    pp.conjugate();
  }
  hash_value_ = hash::value(pp);
}

EvalExpr::EvalExpr(EvalOp op, ResultType res, ExprPtr const& ex,
                   index_vector ixs, CanonTransform transform, size_t h,
                   std::shared_ptr<bliss::Graph> connectivity)
    : op_type_{op},
      result_type_{res},
      expr_{ex.clone()},
      canon_indices_{std::move(ixs)},
      canon_transform_{transform},
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

size_t EvalExpr::hash_value() const noexcept { return hash_value_; }

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
  } else if (op_type() == EvalOp::RealPart) {
    return "Re";
  } else if (op_type() == EvalOp::ImagPart) {
    return "Im";
  } else if (is_variable()) {
    return toUtf8(as_variable().label());
  } else {
    SEQUANT_ABORT("EvalExpr::label: unhandled expression type");
  }
}

std::int8_t EvalExpr::canon_phase() const noexcept {
  return canon_transform_.phase;
}

ExprPtr EvalExpr::denoted_expr() const {
  SEQUANT_ASSERT(is_tensor());
  auto t = expr_->as<Tensor>();
  auto const tr = canon_transform_;
  // Re-materialize the transform syntactically: the slot swap, then (for a
  // Kramers-folded leaf) the flavor flip, whose own marker toggle is already
  // inside tr.conj and is undone by toggling once more. The conj bit itself
  // is spelled as the marker WHETHER OR NOT it came with the swap: a
  // Hermitian leaf written in its non-canonical orientation denotes as
  // C^*{swapped} -- value-equivalent to the as-written spelling only up to
  // the Hermiticity the parent network does not see -- because the marker
  // COLORS the parent's graph, which is what keeps mixed products such as
  // C.C^* and C.C identity-distinct and their canonical layouts right.
  // (Spelling the swapped Hermitian leaf unmarked, i.e. taking the marker
  // out for every channel that came with a conj, was tried on 2026-09-03:
  // it broke the HSeOH PNS-CCD residual in iteration 2 -- a cached
  // intermediate served in the wrong layout -- while the unit suites stayed
  // green.)
  if (tr.braket_swap) static_cast<AbstractTensor&>(t)._swap_bra_ket();
  if (kramers_folded_) kramers_flip_slots(static_cast<AbstractTensor&>(t));
  if (tr.conj != kramers_folded_) t.conjugate();
  return ex<Tensor>(std::move(t));
}

CanonTransform EvalExpr::canon_transform() const noexcept {
  return canon_transform_;
}

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
  // the conjugation marker NEVER enters the slot hash: slot identity is the
  // canonical spelling; a value-distinctive marker (Nonsymm) rides in the
  // leaf's CanonTransform and salts the PARENT's structural hash instead
  return h;
}
}  // namespace

/// Prefix slot hashes of a range of summand hashes: element i is the hash of
/// summands 0..i IN ORDER, i.e. the slot of the sum node that combines them
/// (the left fold's i-th intermediate). Order-sensitive on purpose: a sum
/// hands up its FIRST summand's layout (the others are permuted into it at
/// evaluation), so A + B and B + A are different slots -- while relabeled
/// copies of the same ordered sum (isomorphic summands, hence isomorphic
/// leading layouts) share one. Computed eagerly and sequentially: a lazily
/// sliced prefix view evaluated by random access saw only the summands
/// BEFORE the last one, so A + B and A + C shared a slot (fixed 2026-09-03).
template <typename Rng>
container::svector<size_t> imed_hashes(Rng const& rng) {
  container::svector<size_t> result;
  container::svector<size_t> prefix;
  for (auto&& h : rng) {
    prefix.push_back(h);
    result.push_back(hash::range(prefix.begin(), prefix.end()));
  }
  return result;
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
/// a child's structural identity: slot hash + conj/swap salt (0 for a
/// trivial transform, keeping marker-free hashes byte-stable)
inline size_t salted_hash(EvalExprNode const& n) {
  auto h = n->hash_value();
  if (auto salt = n->canon_transform().structural_salt(); salt != 0)
    hash::combine(h, salt);
  return h;
}

template <typename Rng>
void collect_tensor_factors(EvalExprNode const& node,  //
                            Rng& collect) {
  static_assert(std::is_same_v<ranges::range_value_t<Rng>, ExprWithHash>);

  if (auto op = node->op_type();
      node->is_tensor() && (!op || *op == EvalOp::Sum)) {
    // Leaf tensors enter in their DENOTED spelling (transform re-materialized
    // syntactically); a Sum-rooted subtree contributes its result tensor.
    auto e = (!op && node->expr()->is<Tensor>()) ? node->denoted_expr()
                                                 : node->expr();
    collect.emplace_back(ExprWithHash{.expr = std::move(e),  //
                                      .hash = salted_hash(node)});
  } else if (node->op_type() == EvalOp::Product && !node.leaf()) {
    collect_tensor_factors(node.left(), collect);
    collect_tensor_factors(node.right(), collect);
  }
}

namespace {
EvalExprNode binarize_re_im(ExprPtr const& inner, EvalOp op,
                            IndexSet const& uncontract,
                            const BinarizationOptions& opts,
                            std::size_t& node_counter, bool shared_counter);
}  // namespace

EvalExprNode binarize(Constant const& c) { return EvalExprNode{EvalExpr{c}}; }

EvalExprNode binarize(Variable const& v) { return EvalExprNode{EvalExpr{v}}; }

EvalExprNode binarize(Power const& p) { return EvalExprNode{EvalExpr{p}}; }

EvalExprNode binarize(Tensor const& t) {
  // Every conjugation channel ('⁺' adjoint label, elementwise-conjugation
  // marker, Conjugate-braket orientation) is normalized by the EvalExpr leaf
  // ctor into the canonical unmarked spelling plus a CanonTransform served
  // on retrieval -- a tensor leaf is always just a leaf.
  EvalExpr leaf{t};
  // whose spelling normalization produced a non-trivial retrieval transform.
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

  // uniform-conj hoisting: a sum whose EVERY summand carries conj equals
  // conj of the unconjugated sum -- strip the conj salts so the slot hash
  // matches, and record {conj} on the sum nodes
  bool const hoist_conj =
      !ranges::empty(summands) && ranges::all_of(summands, [](auto&& n) {
        return n->canon_transform().conj;
      });
  // a summand's phase (its hand-up is the DENOTED value) hoists the same
  // way: uniform -> a whole-node phase; mixed -> no transform of the
  // unsigned sum (A - B vs A + B), so the negated summands salt the slot
  auto const negated = [](auto&& n) { return n->canon_phase() == -1; };
  bool const hoist_phase =
      !ranges::empty(summands) && ranges::all_of(summands, negated);
  bool const mixed_phase = !hoist_phase && ranges::any_of(summands, negated);
  auto hvals = summands | transform([hoist_conj, mixed_phase](auto&& n) {
                 auto h = n->hash_value();
                 auto tr = n->canon_transform();
                 if (hoist_conj) tr.conj = false;
                 if (auto salt = tr.structural_salt(); salt != 0)
                   hash::combine(h, salt);
                 if (mixed_phase && n->canon_phase() == -1)
                   hash::combine(h, CanonTransform::phase_salt);
                 return h;
               });
  CanonTransform const sum_transform{
      .phase = static_cast<std::int8_t>(hoist_phase ? -1 : 1),
      .conj = hoist_conj};
  auto make_sum =
      [i = 0, sum_transform,     //
       hs = imed_hashes(hvals),  //
       align = std::size_t{0},   //
       all_tensors,
       &opts](EvalExpr const& left,
              [[maybe_unused]] EvalExpr const& right) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    if (all_tensors) {
      // partition from the DENOTED orientation (stored canonical slots,
      // re-swapped per the child transform)
      auto const t = left.denoted_expr()->as<Tensor>();
      return {
          EvalOp::Sum,         //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(opts, bra(t.bra()), ket(t.ket()),
                                            aux(t.aux())),  //
          left.canon_indices(),                             //
          sum_transform,                                    //
          h,                                                //
          nullptr};
    } else {
      return {EvalOp::Sum,              //
              ResultType::Scalar,       //
              detail::make_variable(),  //
              {},                       //
              sum_transform,            //
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

  // a Re/Im wrapper factor shares the node counter iff every factor is a
  // scalar (wrappers, constants, variables): then the product has no
  // contraction nodes of its own and the wrappers' inner nodes are the
  // summand's DP nodes (mirrors optimize_impl's re-keying)
  const bool wrapper_shares_counter = ranges::all_of(
      prod.factors(), [](ExprPtr const& f) { return f->is_scalar(); });
  auto factors =
      prod.factors()  //
      | transform([i = 0, &ltr_uncontr_idxs, &opts, &node_counter,
                   wrapper_shares_counter](ExprPtr const& x) mutable {
          auto const& uncontr = ltr_uncontr_idxs.children[i++];
          if (x->is<RealPart>())
            return binarize_re_im(x->as<RealPart>().inner(), EvalOp::RealPart,
                                  uncontr, opts, node_counter,
                                  wrapper_shares_counter);
          if (x->is<ImagPart>())
            return binarize_re_im(x->as<ImagPart>().inner(), EvalOp::ImagPart,
                                  uncontr, opts, node_counter,
                                  wrapper_shares_counter);
          if (x->is<Sum>()) {
            // A Sum factor is opaque to the single-term optimizer
            // (opt_mixed_product stands a placeholder tensor in
            // for it, and puts the Sum back untouched), so the
            // contraction nodes INSIDE it are not DP nodes and have
            // no entry in opts.node_batch_axes: binarize them with
            // a private counter and no axes. Consuming the shared
            // counter here shifted every outer node's annotation
            // onto the wrong (inner) node -- e.g. a contracted-axis
            // batch mark onto a node whose result still carries
            // that axis.
            BinarizationOptions inner_opts = opts;
            inner_opts.node_batch_axes.clear();
            std::size_t inner_counter = 0;
            return impl::binarize(x, uncontr, inner_opts, inner_counter);
          }
          return impl::binarize(x, uncontr, opts, node_counter);
        })  //
      | ranges::to_vector;

  // PREFIX-uniform conj hoisting (design spec): the left-fold combines
  // factor prefixes, and a prefix that is uniformly conjugated equals the
  // conj of its unconjugated counterpart -- its node hoists {conj} and its
  // factors' conj salts are stripped, so e.g. the (A^*·B^*) intermediate of
  // A^*·B^*·C is a cache hit on the A·B slot. A broken prefix keeps the
  // salts (mixed marks stay identity-distinct).
  std::vector<char> prefix_conj;
  prefix_conj.reserve(ranges::size(factors) + 1);
  prefix_conj.push_back(false);  // 0-factor prefix
  {
    bool run = !ranges::empty(factors);
    for (auto const& n : factors) {
      run = run && n->canon_transform().conj;
      prefix_conj.push_back(run);
    }
  }
  // prefix hashes with PER-PREFIX conj-salt stripping: factor salts are
  // stripped only inside a prefix that is uniformly conjugated (where the
  // conj hoists onto that prefix's node); in a broken prefix every factor
  // contributes its full salt. Per-prefix (not per-factor) stripping keeps
  // the identity order-insensitive: C·C^* and C^*·C still agree.
  std::vector<size_t> hs;
  {
    auto const n = ranges::size(factors);
    hs.reserve(n);
    std::vector<size_t> buf;
    buf.reserve(n);
    for (std::size_t j = 1; j <= n; ++j) {
      buf.clear();
      bool const strip = static_cast<bool>(prefix_conj[j]);
      std::size_t k = 0;
      for (auto const& fac : factors) {
        if (++k > j) break;
        auto h = fac->hash_value();
        auto tr = fac->canon_transform();
        if (strip) tr.conj = false;
        if (auto salt = tr.structural_salt(); salt != 0) hash::combine(h, salt);
        buf.push_back(h);
      }
      hs.push_back(hash::range_unordered(buf.begin(), buf.end()));
    }
  }

  auto make_prod = [i = 0, &hs, &ltr_uncontr_idxs, &opts, &prefix_conj,
                    &node_counter](
                       EvalExprNode const& left,
                       EvalExprNode const& right) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    // combining the (i+1)-factor prefix: hoist iff that prefix is uniform
    bool const hoist_conj = static_cast<bool>(prefix_conj[i + 1]);
    auto const& uncontracted_idxs = ltr_uncontr_idxs.imed[i];
    if (left->is_scalar() && right->is_scalar()) {
      // scalar * scalar
      return {EvalOp::Product,
              ResultType::Scalar,
              detail::make_variable(),
              {},
              CanonTransform{},
              h,
              nullptr};
    } else if (left->is_scalar() || right->is_scalar()) {
      // scalar * tensor or tensor * scalar
      auto const& tl = left->is_tensor() ? left : right;
      auto const t = tl->denoted_expr()->as<Tensor>();  // denoted orientation
      return {
          EvalOp::Product,     //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(opts, bra(t.bra()), ket(t.ket()),
                                            aux(t.aux())),  //
          tl->canon_indices(),                              //
          tl->canon_transform(),                            //
          h,
          nullptr};
    } else {
      // tensor * tensor
      container::svector<ExprWithHash> subfacs;
      collect_tensor_factors(left, subfacs);
      collect_tensor_factors(right, subfacs);
      // Uniform-conj hoisting (design spec): when EVERY tensor factor's
      // denoted spelling is conjugated, the conjugation is a whole-node
      // transform -- strip the markers (the TN then hashes onto the
      // unconjugated product's slot) and record {conj} on this node. Mixed
      // marks stay in the TN, where the marker coloring keeps e.g. C·C^*
      // identity-distinct from C·C. (Sum-level hoisting: T7 follow-up.)
      if (hoist_conj)
        for (auto& f : subfacs)
          if (f.expr->is<Tensor>() && f.expr->as<Tensor>().conjugated())
            f.expr->as<Tensor>().conjugate();
      auto ts = subfacs | transform([](auto&& t) { return t.expr; });
      IndexGroups<IndexVec> const target_indices = [&ts, &uncontracted_idxs]() {
        // route each surviving hyperindex to its correct slot
        // (bra, ket, or aux) based on which slot it occupies in
        // the factor tensors .. if appears in multiple slots put into aux
        //
        // count on the denoted spellings (collect_tensor_factors already
        // re-materialized each leaf's authored orientation; the conjugation
        // marker does not affect slot occupancy)
        auto counts = get_used_indices_with_counts(ex<Product>(ts));
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
          {.cardinal_tensor_labels =
               TensorCanonicalizer::cardinal_tensor_labels(),
           .named_indices = &named_indices});
      hash::combine(h, canon.hash_value());
      bool const scalar_result = canon.named_indices_canonical.empty();
      // the children hand up their DENOTED values (transforms applied) while
      // the network above is spelled with their canonical slots: their
      // phases hoist multiplicatively onto this node (A * (-B) = -(A * B)),
      // so the slot -- phase-blind, shared by the spellings that differ by
      // a child's antisymmetric reorder -- holds one canonical value
      CanonTransform const transform{
          .phase = static_cast<std::int8_t>(canon.phase * left->canon_phase() *
                                            right->canon_phase()),
          .conj = hoist_conj};
      EvalExpr result =
          scalar_result
              ? EvalExpr{EvalOp::Product,          //
                         ResultType::Scalar,       //
                         detail::make_variable(),  //
                         {},                       //
                         transform,                //
                         h,
                         std::move(canon.graph)}
              : EvalExpr{EvalOp::Product,     //
                         ResultType::Tensor,  //
                         detail::make_tensor_wo_symmetries(
                             opts, bra(target_indices.bra),
                             ket(target_indices.ket), aux(target_indices.aux)),
                         canon.get_indices<Index::index_vector>(),  //
                         transform,                                 //
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
                    ? detail::make_tensor(left->denoted_expr()->as<Tensor>(),
                                          false, opts)
                : left->is_constant() ? (left->expr() * right->expr())
                                      : detail::make_variable();
    auto type = left->is_tensor() ? ResultType::Tensor : ResultType::Scalar;

    // a REAL scalar commutes with conj, so a conj-hoisted subtree hoists
    // through the wrap too (\mathcal{T}-partner terms carry real prefactors)
    bool const wrap_hoist = left->canon_transform().conj &&
                            right->is_constant() &&
                            right->as_constant().value().imag() == 0;
    auto h = left->hash_value();
    {
      auto tr = left->canon_transform();
      if (wrap_hoist) tr.conj = false;
      if (auto salt = tr.structural_salt(); salt != 0) hash::combine(h, salt);
    }
    hash::combine(h, salted_hash(right));
    auto result = EvalExpr{EvalOp::Product,          //
                           type,                     //
                           expr,                     //
                           left->canon_indices(),    //
                           left->canon_transform(),  //
                           h,                        //
                           nullptr};

    return EvalExprNode{std::move(result), std::move(left), std::move(right)};
  }
}

namespace {

// Unary Re/Im wrapper over the shared inner subtree: EvalOp::RealPart or
// ImagPart, ResultType::Scalar, Constant{1} sentinel right child (the
// FullBinaryNode "every non-leaf has two children" invariant; evaluate
// ignores it for these ops). Node hash = the inner child's salted hash
// combined with the op, so Re(s), Im(s) and bare s occupy distinct slots
// while the inner subtree itself stays on its own shared slot.
EvalExprNode binarize_re_im(ExprPtr const& inner, EvalOp op,
                            IndexSet const& uncontract,
                            const BinarizationOptions& opts,
                            std::size_t& node_counter, bool shared_counter) {
  // A wrapper at the summand root, or whose product siblings are all
  // scalars, has its inner contraction nodes as the summand's DP nodes
  // (the optimizer optimizes the inner and records its batch axes under the
  // summand): share the node counter. Otherwise the inner is opaque to the
  // single-term optimizer (like a Sum factor) -- private counter, no axes.
  EvalExprNode inner_node = [&]() {
    if (shared_counter)
      return impl::binarize(inner, uncontract, opts, node_counter);
    BinarizationOptions inner_opts = opts;
    inner_opts.node_batch_axes.clear();
    std::size_t inner_counter = 0;
    return impl::binarize(inner, uncontract, inner_opts, inner_counter);
  }();
  auto h = inner_node->hash_value();
  if (auto salt = inner_node->canon_transform().structural_salt(); salt != 0)
    hash::combine(h, salt);
  hash::combine(h, static_cast<size_t>(op));
  // the inner phase hoists through the projection (Re(-x) = -Re(x))
  EvalExpr wrap{op,
                ResultType::Scalar,
                detail::make_variable(),
                {},
                CanonTransform{.phase = inner_node->canon_phase()},
                h,
                nullptr};
  EvalExprNode sentinel{EvalExpr{Constant{1}}};
  return EvalExprNode{std::move(wrap), std::move(inner_node),
                      std::move(sentinel)};
}

}  // namespace

namespace impl {

EvalExprNode binarize(ExprPtr const& expr, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter) {
  if (expr->is<RealPart>())
    return binarize_re_im(expr->as<RealPart>().inner(), EvalOp::RealPart,
                          uncontract, opts, node_counter,
                          /*shared_counter=*/true);

  if (expr->is<ImagPart>())
    return binarize_re_im(expr->as<ImagPart>().inner(), EvalOp::ImagPart,
                          uncontract, opts, node_counter,
                          /*shared_counter=*/true);

  if (expr->is<Constant>())  //
    return binarize(expr->as<Constant>());

  if (expr->is<Variable>())  //
    return binarize(expr->as<Variable>());

  if (expr->is<Tensor>())  //
    return binarize(expr->as<Tensor>());

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
