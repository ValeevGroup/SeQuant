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
#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/equal.hpp>
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
#include <cstdlib>
#include <iostream>
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

std::size_t EvalExpr::layout_fingerprint() const noexcept {
  if (layout_fingerprint_) return *layout_fingerprint_;
  // A Kramers-folded leaf's slot holds the FOLDED (up-row) spelling, i.e.
  // expr(), while canon_indices() keeps the as-written flavors (see the ctor):
  // fingerprint the layout the slot carries, so the down-first leaf shares
  // its up-first partner's slot (the flavor difference rides the transform).
  if (kramers_folded_) {
    index_vector folded = identity_indices();
    if (const auto isr = get_default_context().index_space_registry())
      kramers_flip(folded, *isr);
    layout_fingerprint_ = layout_fingerprint_of(folded);
  } else {
    layout_fingerprint_ = layout_fingerprint_of(identity_indices());
  }
  return *layout_fingerprint_;
}

std::size_t EvalExpr::layout_fingerprint_of(
    index_vector const& modes) noexcept {
  container::map<Index, std::size_t> ids;
  auto id_of = [&ids](Index const& ix) {
    return ids.try_emplace(ix, ids.size()).first->second;
  };
  std::size_t fp = 0;
  // Walk the modes in the order the RESULT carries them, which is the order
  // indices_annot() builds: the proto-free (outer) indices in canonical order,
  // then the proto-carrying (inner) ones. Using the raw canon_indices() order
  // instead would separate nodes whose two groups interleave differently while
  // laying their modes out identically, and needlessly cost cache sharing.
  auto walk = [&](bool proto) {
    for (auto const& ix : modes) {
      if (ix.has_proto_indices() != proto) continue;
      hash::combine(fp, static_cast<std::int64_t>(ix.space().attr()));
      hash::combine(fp, id_of(ix));
      hash::combine(fp, ix.proto_indices().size());
      for (auto const& p : ix.proto_indices()) hash::combine(fp, id_of(p));
    }
  };
  walk(false);
  walk(true);
  return fp;
}

EvalExpr::index_vector const& EvalExpr::canon_indices() const noexcept {
  return canon_indices_;
}

void EvalExpr::set_identity_indices(index_vector ixs) {
  identity_indices_ = std::move(ixs);
  layout_fingerprint_.reset();
  refresh_identity_tensor();
}

void EvalExpr::refresh_identity_tensor() {
  identity_tensor_.reset();
  if (identity_indices_.empty() || !expr_ || !expr_->is<Tensor>()) return;
  SEQUANT_ASSERT(identity_indices_.size() == canon_indices_.size());
  // the erased spelling of the result tensor (graph-less nodes -- Sum roots,
  // scalar * tensor -- are compared by block on it, see
  // TreeNodeEqualityComparator)
  container::map<Index, Index> repl;
  for (std::size_t k = 0; k < canon_indices_.size(); ++k)
    if (canon_indices_[k] != identity_indices_[k])
      repl.emplace(canon_indices_[k], identity_indices_[k]);
  if (repl.empty()) return;
  Tensor t = expr_->as<Tensor>();
  t.transform_indices(repl);
  t.reset_tags();
  identity_tensor_ = std::move(t);
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
  // a leaf with a Kramers-union slot has no elementwise time-reversal image
  // (see kramers_union_index): served as written, never folded
  const int kramers_phase = (fold_kramers_leaf() && !has_kramers_union_slot(t))
                                ? canonicalize_kramers(t)
                                : 1;
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

EvalExpr::EvalExpr(Tensor const& tnsr, eval::KramersBlindness const* blindness)
    : op_type_{std::nullopt},
      result_type_{ResultType::Tensor},
      expr_{tnsr.clone()} {
  SEQUANT_ASSERT(!tnsr.indices().empty());
  // the stored spelling is canonical (unmarked, up-row, block-canonical);
  // every spelling channel becomes a CanonTransform byproduct applied on
  // retrieval, so every route to one canonical spelling lands on one slot
  auto const [transform, kramers_fired] = normalize_leaf(expr_->as<Tensor>());
  canon_transform_ = transform;
  // Kramers-blind identity (kramers_blind.hpp): the indices of the normalized
  // spelling whose flavour the served value does not depend on. Empty (the
  // default, and whenever the hook is inactive) => identity exactly as below.
  eval::ErasureMap erasure;
  if (blindness && blindness->active())
    erasure = eval::erasure_map(std::array{expr_}, *blindness);
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
    // array-faithful indices in the Nested (outer;inner) convention, i.e. the
    // layout a leaf provider serves for the STORED spelling: the outer modes
    // are the pure proto indices (those not occupying a slot of their own,
    // e.g. the pair labels of C{mu;a<ij>}, laid out i,j,mu) followed by the
    // plain slots IN SLOT ORDER, the inner modes the proto-carrying slots in
    // slot order. A plain slot that is also a proto index (the occupied kets
    // of t{a<ij>,b<ij>;i,j}) takes its SLOT position: the array's outer mode
    // k pairs with inner mode k as a column, so t{a<ij>,b<ij>;j,i} must be
    // read j,i;a,b -- the proto-first order (tot_indices) would read the same
    // array i,j;a,b and serve t^{ab}_{ij} for t^{ab}_{ji}. (The md list is the
    // same set in named-canonical order, which annots must not depend on.)
    auto const slot_ixs =
        expr_->as<Tensor>().const_indices() | ranges::to<index_vector>;
    canon_indices_.clear();
    auto const is_plain_slot = [&slot_ixs](Index const& p) {
      return ranges::any_of(slot_ixs, [&p](Index const& s) {
        return !s.has_proto_indices() && s == p;
      });
    };
    for (auto const& ix : slot_ixs)
      for (auto const& p : ix.proto_indices())
        if (!is_plain_slot(p) && !ranges::contains(canon_indices_, p))
          canon_indices_.emplace_back(p);
    for (auto const& ix : slot_ixs)
      if (!ix.has_proto_indices()) canon_indices_.emplace_back(ix);
    for (auto const& ix : slot_ixs)
      if (ix.has_proto_indices()) canon_indices_.emplace_back(ix);
    connectivity_ = std::move(md.graph);
    if (!erasure.empty()) {
      // hash and graph from the erased spelling of the STORED (canonical)
      // tensor; expr_ and canon_indices_ keep the as-written flavours. The
      // erased network must canonicalize to the same slot order as the real
      // one (the slot holds the real layout): if it does not, keep the
      // unerased identity -- a missed fold, never a transposed one.
      auto erased =
          ex<Tensor>(eval::erase_indices(expr_->as<Tensor>(), erasure));
      ExprPtrList elist{erased};
      auto etn = TensorNetwork(elist);
      auto emd = etn.canonicalize_slots(
          {.cardinal_tensor_labels =
               TensorCanonicalizer::cardinal_tensor_labels(),
           .apply_slot_order = true});
      auto const ecanon_e =
          std::dynamic_pointer_cast<Expr>(etn.tensors().front());
      SEQUANT_ASSERT(ecanon_e && ecanon_e->is<Tensor>());
      auto const& ecanon = ecanon_e->as<Tensor>();
      bool const same_order = ranges::equal(
          ecanon.const_slots(), erased->as<Tensor>().const_slots(),
          [](Index const& a, Index const& b) {
            return a.full_label() == b.full_label();
          });
      if (same_order) {
        hash_value_ = emd.hash_value();
        connectivity_ = std::move(emd.graph);
        identity_tensor_ = erased->as<Tensor>();
        identity_indices_ = canon_indices_;
        for (auto& ix : identity_indices_) ix = eval::erase_index(ix, erasure);
      }
    }
  } else {
    auto const& t = expr_->as<Tensor>();
    canon_indices_ = t.const_indices() | ranges::to<index_vector>;
    if (!erasure.empty()) {
      identity_tensor_ = eval::erase_indices(t, erasure);
      hash_value_ = hash_terminal_tensor(*identity_tensor_);
      identity_indices_ =
          identity_tensor_->const_indices() | ranges::to<index_vector>;
    } else {
      hash_value_ = hash_terminal_tensor(t);
    }
  }
  // T19 layer 2 contract: expr() keeps the FOLDED (up-row) spelling -- what
  // a leaf provider fetches -- while canon_indices() (the parent's
  // contraction labels; TA matches annotations, not spellings) carries the
  // as-written flavors in the same canonical order, so the served up block
  // + {conj, phase} denotes the as-written value
  if (kramers_fired) kramers_flip_indices_as_written(canon_indices_);
  kramers_folded_ = kramers_fired;
  fold_layout_into_hash();
}

void EvalExpr::fold_layout_into_hash() noexcept {
  // The slot identity (hash) is label-blind and orientation-blind by design,
  // so two nodes that lay their result modes out differently -- same-space
  // externals of an isomorphic network ordered either way (bliss breaks an
  // automorphic orbit by input vertex order), or a Sum handing up another
  // summand's layout -- would share it, and a cached buffer served under the
  // other layout is a transposed value. Fold the (renaming-invariant) layout
  // fingerprint in, so hash equality implies layout equality: the hash-keyed
  // value maps of the ordered (DAG) executor and the equality comparator then
  // agree on what is one value.
  if (result_type_ == ResultType::Tensor)
    hash::combine(hash_value_, layout_fingerprint());
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

bool EvalExpr::is_unary_op() const noexcept {
  auto const op = op_type();
  return op == EvalOp::RealPart || op == EvalOp::ImagPart ||
         op == EvalOp::KramersFlip;
}

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
  // prefix hashes over the summands as a MULTISET (order-independent): a Sum
  // of the same summands in another order is the same value. The layout a
  // Sum hands up (its first summand's, see binarize(Sum)) is separated by the
  // layout fingerprint every tensor-valued node folds into its hash, so two
  // orders that lay their result out differently still get distinct slots.
  container::svector<size_t> result;
  container::svector<size_t> prefix;
  for (auto&& h : rng) {
    prefix.push_back(h);
    result.push_back(hash::range_unordered(prefix.begin(), prefix.end()));
  }
  return result;
}

struct ExprWithHash {
  ExprPtr expr;
  size_t hash;
  /// the factor's own canonical phase: the part of its orientation that its
  /// spelling in the network cannot carry (an index reorder), see
  /// collect_tensor_factors
  std::int8_t phase = 1;
  /// whether the factor is a leaf (false: a Sum-rooted intermediate); see
  /// eval::erasable_indices
  bool leaf = true;
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
      node->is_tensor() &&
      (!op || *op == EvalOp::Sum || *op == EvalOp::KramersFlip)) {
    // Leaf tensors enter in their DENOTED spelling (transform re-materialized
    // syntactically); a Sum-rooted subtree contributes its result tensor, and
    // so does a KramersFlip wrapper (its expr() is the flipped-flavour
    // spelling the parent contracts through).
    auto e = (!op && node->expr()->is<Tensor>()) ? node->denoted_expr()
                                                 : node->expr();
    // The spelling carries the factor's conj / bra-ket swap but not its
    // reorder phase (a sign is not a spelling), and a Sum root enters in its
    // slot spelling outright: the phase rides along for the product fold.
    collect.emplace_back(ExprWithHash{.expr = std::move(e),  //
                                      .hash = salted_hash(node),
                                      .phase = node->canon_phase(),
                                      .leaf = !op});
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

EvalExprNode binarize(Tensor const& t, const BinarizationOptions& opts) {
  // Every conjugation channel ('⁺' adjoint label, elementwise-conjugation
  // marker, Conjugate-braket orientation) is normalized by the EvalExpr leaf
  // ctor into the canonical unmarked spelling plus a CanonTransform served
  // on retrieval -- a tensor leaf is always just a leaf.
  EvalExpr leaf{t, &opts.kramers_blindness};
  // whose spelling normalization produced a non-trivial retrieval transform.
  return EvalExprNode{std::move(leaf)};
}

EvalExprNode binarize(Sum const& sum, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter) {
  using ranges::views::move;
  using ranges::views::transform;
  // SEQUANT_BATCH_AXES_DEBUG=1: one line per summand with the contraction
  // nodes it consumed from opts.node_batch_axes, to align with the
  // optimizer's per-summand entry counts (optimize.cpp rekey_onto)
  static const bool debug = std::getenv("SEQUANT_BATCH_AXES_DEBUG");
  std::size_t smand_idx = 0;
  // a summand is added elementwise into the Sum's layout, so it must keep
  // the Sum's labels: never a KramersFlip as a whole (the Sum folds as one)
  BinarizationOptions sopts = opts;
  sopts.kramers_fold_this = false;
  auto summands =
      sum.summands()  //
      | transform([&uncontract, &opts, &sopts, &node_counter,
                   &smand_idx](ExprPtr const& x) {
          std::size_t const before = node_counter;
          auto node = impl::binarize(x, uncontract, sopts, node_counter);
          if (debug && !opts.node_batch_axes.empty())
            std::cerr << "[batch-axes] binarize summand " << smand_idx << ": "
                      << (node_counter - before) << " nodes"
                      << " type="
                      << (x->is<Product>() ? "Product"
                          : x->is<Sum>()   ? "Sum"
                                           : "other")
                      << " | " << toUtf8(x->to_latex()).substr(0, 160) << "\n";
          ++smand_idx;
          return node;
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
  // Every binary Sum produced by fold_left_to_node below folds the running
  // accumulator (the chain seed, or a prior chain Sum) in as the left
  // operand (see fold_left_to_node in binary_node.hpp: the accumulator is
  // always `l`), so every chain Sum node accumulates its left operand in
  // place -- see EvalExpr::accumulate_in_place.
  auto make_sum = [i = 0, sum_transform,                         //
                   hs = imed_hashes(hvals) | ranges::to_vector,  //
                   all_tensors, &opts](EvalExpr const& left,
                                       EvalExpr const&) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    if (all_tensors) {
      // partition from the DENOTED orientation (stored canonical slots,
      // re-swapped per the child transform)
      auto const t = left.denoted_expr()->as<Tensor>();
      // The placeholder is what an enclosing tensor network sees for this
      // (opaque) node, so spell its slots -- within each bra/ket/aux group --
      // in the node's canonical index order, i.e. the order the value is laid
      // out in. Sorted by label instead, two relabeled spellings of one sum
      // (values: transposes of each other, e.g. the bra slots of a
      // column-symmetric g projected by C's onto a_1,a_2 vs a_2,a_1) spelled
      // one identical placeholder, and an enclosing product got one hash AND
      // one canonical layout for both: the cache served one value for the
      // other untransposed (HSeOH PNS-CCD, (vv|vv) kept 4-center; test
      // sum_placeholder_is_spelled_in_its_layout).
      auto const& frame = left.canon_indices();
      auto in_canon_order = [&frame](auto const& group) {
        Index::index_vector ordered;
        for (auto const& ix : frame)
          if (std::find(group.begin(), group.end(), ix) != group.end())
            ordered.emplace_back(ix);
        SEQUANT_ASSERT(ordered.size() == ranges::size(group));
        return ordered;
      };
      // the layout is part of the value's identity (see layout_fingerprint):
      // the same summands led by a differently laid-out summand are another
      // slot, or a cached array would be served in the wrong mode order
      // the sum hands up its first summand's layout, identity layout included
      hash::combine(h,
                    EvalExpr::layout_fingerprint_of(left.identity_indices()));
      EvalExpr result{
          EvalOp::Sum,         //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(
              opts, bra(in_canon_order(t.bra())), ket(in_canon_order(t.ket())),
              aux(in_canon_order(t.aux())), /*keep_order=*/true),  //
          left.canon_indices(),                                    //
          sum_transform,                                           //
          h,                                                       //
          nullptr};
      if (left.has_identity_erasure())
        result.set_identity_indices(left.identity_indices());
      result.set_accumulate_in_place(true);
      return result;
    } else {
      EvalExpr result{EvalOp::Sum,              //
                      ResultType::Scalar,       //
                      detail::make_variable(),  //
                      {},                       //
                      sum_transform,            //
                      h,                        //
                      nullptr};
      result.set_accumulate_in_place(true);
      return result;
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
  // a factor is contracted through its own denoted labels, so it may fold
  // as a whole even inside a summand that may not
  BinarizationOptions fopts = opts;
  fopts.kramers_fold_this = true;
  auto factors =
      prod.factors()  //
      | transform([i = 0, &ltr_uncontr_idxs, &opts = fopts, &node_counter,
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
      // the identity layout follows the tensor operand's (Kramers-blind
      // erased where it was)
      hash::combine(h, EvalExpr::layout_fingerprint_of(tl->identity_indices()));
      EvalExpr result{
          EvalOp::Product,     //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(opts, bra(t.bra()), ket(t.ket()),
                                            aux(t.aux())),  //
          tl->canon_indices(),                              //
          tl->canon_transform(),                            //
          h,
          nullptr};
      if (tl->has_identity_erasure())
        result.set_identity_indices(tl->identity_indices());
      return result;
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

      // Kramers-blind identity (kramers_blind.hpp): the hash, the
      // connectivity graph, the canonical result layout and the phase all
      // come from the ERASED flattened network; only the result labels are
      // mapped back to the as-written indices (the erasure is a renaming, so
      // the map is a bijection). target_indices above keep the real slots.
      eval::ErasureMap erasure;
      if (opts.kramers_blindness.active()) {
        container::svector<bool> leaf_flags;
        for (auto const& f : subfacs) leaf_flags.push_back(f.leaf);
        erasure = eval::erasure_map(ts, opts.kramers_blindness, leaf_flags);
      }
      container::svector<ExprPtr> ts_id;
      for (ExprPtr const& e : ts)
        ts_id.push_back(erasure.empty() ? e
                                        : ex<Tensor>(eval::erase_indices(
                                              e->as<Tensor>(), erasure)));
      auto tn = TensorNetwork(ts_id);
      auto named_indices = tn.ext_indices();
      for (auto&& ix : uncontracted_idxs)
        named_indices.emplace(eval::erase_index(ix, erasure));

      auto canon = tn.canonicalize_slots(
          {.cardinal_tensor_labels =
               TensorCanonicalizer::cardinal_tensor_labels(),
           .named_indices = &named_indices});
      hash::combine(h, canon.hash_value());
      bool const scalar_result = canon.named_indices_canonical.empty();
      // The network above is the FLATTENED one: every leaf (and Sum root)
      // under this product enters in its slot spelling
      // (collect_tensor_factors), so canon.phase already relates the value this
      // node computes -- the contraction of those spellings, phases aside -- to
      // the canonical network's. The one thing a spelling cannot carry is a
      // factor's own reorder phase, so those hoist multiplicatively onto this
      // node (A * (-B) = -(A * B)). A Product child's own phase is NOT
      // re-applied: its sub-network is part of the flattened one, so its
      // reorder sign is inside canon.phase already; multiplying it again made
      // two spellings of one slot disagree by a sign whenever the child's
      // canonical phase was -1 (cache-only symptom: the slot served -R to one
      // of them).
      std::int8_t factor_phase = 1;
      for (auto const& f : subfacs)
        factor_phase = static_cast<std::int8_t>(factor_phase * f.phase);
      CanonTransform const transform{
          .phase = static_cast<std::int8_t>(canon.phase * factor_phase),
          .conj = hoist_conj};
      auto result_indices = canon.get_indices<Index::index_vector>();
      // the result layout is part of a tensor-valued node's identity (see
      // layout_fingerprint); under erasure the identity layout is the erased
      // one and the node's own labels are the as-written ones
      Index::index_vector identity_indices;
      if (!erasure.empty()) {
        identity_indices = result_indices;
        eval::ErasureMap unerase;
        for (auto const& [real, placeholder] : erasure)
          unerase.emplace(placeholder, real);
        for (auto& ix : result_indices) ix = eval::erase_index(ix, unerase);
      }
      if (!scalar_result)
        hash::combine(h,
                      EvalExpr::layout_fingerprint_of(
                          erasure.empty() ? result_indices : identity_indices));
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
                         std::move(result_indices),  //
                         transform,                  //
                         h,
                         std::move(canon.graph)};
      if (!scalar_result && !erasure.empty())
        result.set_identity_indices(std::move(identity_indices));
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
        result.set_node_slice_mask(ann.axes);
        result.set_batch_loops_opened_here(ann.opened_here);
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

/// slot occurrences (bra / ket / aux, not protos) of every index across the
/// leaves of \p e; a Sum contributes its first summand (every summand has the
/// same externals), a Re/Im wrapper nothing (scalar-valued: its indices are
/// all contracted inside it)
void count_slot_occurrences(ExprPtr const& e, container::map<Index, int>& cnt) {
  if (e->is<Tensor>()) {
    auto const& t = e->as<Tensor>();
    for (auto const& ix : t.bra()) ++cnt[ix];
    for (auto const& ix : t.ket()) ++cnt[ix];
    for (auto const& ix : t.aux()) ++cnt[ix];
  } else if (e->is<Sum>()) {
    auto const& s = e->as<Sum>().summands();
    if (!s.empty()) count_slot_occurrences(s.front(), cnt);
  } else if (e->is<Product>()) {
    for (auto const& f : e->as<Product>().factors())
      count_slot_occurrences(f, cnt);
  }
}

/// the flavour-blind key of an index: its full label (protos included) with
/// the Kramers flavour marks removed, shared by the two partners of a pair
std::wstring flavour_blind_label(Index const& ix) {
  std::wstring s(ix.full_label());
  std::erase_if(s, [](wchar_t c) { return c == L'↑' || c == L'↓'; });
  return s;
}

/// Phase 2a (mpqc doc/dev/specs/2026-09-18-union-axis-time-reversal-fold.md):
/// a Product / Sum whose flavoured externals (after the Kramers-blind erasure
/// of the pair labels) are down-majority -- a tie resolved by the flavour
/// string in the flavour-blind canonical order of the externals, so exactly
/// one of two partners folds -- and whose leaves are all time-reversal
/// symmetric is the time-reversal image of its flipped spelling: the flipped
/// expression (the canonical partner, shared with the partner family) is
/// binarized and wrapped in a KramersFlip over the union free legs with phase
/// (-1)^{n_down}. Null when the fold does not apply.
std::optional<EvalExprNode> maybe_kramers_fold(ExprPtr const& expr,
                                               IndexSet const& uncontract,
                                               BinarizationOptions const& opts,
                                               std::size_t& node_counter) {
  auto const isr = get_default_context().index_space_registry();
  if (!isr) return std::nullopt;
  container::svector<ExprPtr> leaves;
  bool all_tr = true;
  expr->visit(
      [&](ExprPtr const& x) {
        if (!x->is<Tensor>()) return;
        leaves.push_back(x);
        all_tr = all_tr && x->as<Tensor>().kramers_symmetry() ==
                               KramersSymmetry::TimeReversal;
      },
      /*atoms_only=*/true);
  if (!all_tr || leaves.empty()) return std::nullopt;

  container::map<Index, int> cnt;
  count_slot_occurrences(expr, cnt);
  // the blind pair labels (kramers_blind.hpp): the node's value does not
  // depend on their flavour, so they are not flavoured externals here
  auto const blind = opts.kramers_blindness.active()
                         ? eval::blind_indices(leaves, opts.kramers_blindness)
                         : container::set<Index>{};
  auto const down = [&isr](Index const& ix) {
    return !isr->kramers_canonical(ix.space());
  };
  // the flavoured externals: used once across the leaves or kept
  // uncontracted; the union legs, the spin-free indices and the blind pair
  // labels do not count
  container::svector<Index> flav;
  for (auto const& [ix, n] : cnt) {
    if (n != 1 && !uncontract.contains(ix)) continue;
    // a nested (proto-carrying) union axis cannot be flipped by
    // Result::kramers_flip, which acts on the outer modes
    if (ix.has_proto_indices() && kramers_union_index(ix, *isr))
      return std::nullopt;
    if (blind.contains(ix) || !isr->kramers_partner(ix.space())) continue;
    flav.push_back(ix);
  }
  if (flav.empty()) return std::nullopt;
  auto const n_down = std::count_if(flav.begin(), flav.end(), down);
  auto const n_up = static_cast<std::ptrdiff_t>(flav.size()) - n_down;
  bool noncanonical = n_down > n_up;
  if (n_down == n_up) {
    std::sort(flav.begin(), flav.end(), [](Index const& x, Index const& y) {
      return flavour_blind_label(x) < flavour_blind_label(y);
    });
    std::wstring own, flipped;
    for (auto const& ix : flav) {
      own += down(ix) ? L'b' : L'a';
      flipped += down(ix) ? L'a' : L'b';
    }
    noncanonical = flipped < own;
  }
  if (!noncanonical) return std::nullopt;

  // the canonical partner: every flavoured index occurrence flipped, the
  // pair labels inside unflavoured (union) composites included
  auto flipped_expr = expr->clone();
  flipped_expr->visit(
      [](ExprPtr& x) {
        if (!x->is<Tensor>()) return;
        auto& t = x->as<Tensor>();
        kramers_flip_slots_deep(t);
        t.reset_tags();
      },
      /*atoms_only=*/true);
  IndexSet flipped_uncontract;
  for (auto const& ix : uncontract)
    flipped_uncontract.emplace(kramers_flipped_deep(ix, *isr));
  auto inner =
      impl::binarize(flipped_expr, flipped_uncontract, opts, node_counter);
  SEQUANT_ASSERT(inner->is_tensor());
  // the union free legs as OUTER mode positions (the proto-free indices in
  // canonical order, see EvalExpr::indices_annot)
  container::svector<std::size_t> modes;
  std::size_t k = 0;
  for (auto const& ix : inner->canon_indices()) {
    if (ix.has_proto_indices()) continue;
    if (kramers_union_index(ix, *isr)) modes.push_back(k);
    ++k;
  }
  auto const phase = static_cast<std::int8_t>((n_down % 2) ? -1 : 1);
  Tensor denoted = inner->as_tensor();
  kramers_flip_slots_deep(denoted);
  denoted.reset_tags();
  return make_kramers_flip_node(std::move(inner), std::move(modes), phase,
                                std::move(denoted));
}

}  // namespace

EvalExprNode make_kramers_flip_node(EvalExprNode inner,
                                    container::svector<std::size_t> modes,
                                    std::int8_t phase, Tensor denoted) {
  SEQUANT_ASSERT(inner->result_type() == ResultType::Tensor);
  auto const isr = get_default_context().index_space_registry();
  auto h = inner->hash_value();
  if (auto salt = inner->canon_transform().structural_salt(); salt != 0)
    hash::combine(h, salt);
  hash::combine(h, static_cast<size_t>(EvalOp::KramersFlip));
  for (auto m : modes) hash::combine(h, m);
  hash::combine(h, static_cast<std::int64_t>(phase));
  // the wrapper's labels: the child's layout with every flavoured index
  // flipped (the denoted spelling)
  auto ixs = inner->canon_indices();
  if (isr) kramers_flip_deep(ixs, *isr);
  hash::combine(h, EvalExpr::layout_fingerprint_of(ixs));
  // the child's phase / conj hoist through F (linear, commutes with conj), as
  // the inner phase hoists through Re/Im in binarize_re_im
  auto const& itr = inner->canon_transform();
  EvalExpr wrap{EvalOp::KramersFlip,
                ResultType::Tensor,
                ex<Tensor>(std::move(denoted)),
                std::move(ixs),
                CanonTransform{.phase = itr.phase, .conj = itr.conj},
                h,
                nullptr};
  wrap.set_kramers_flip(std::move(modes), phase);
  if (inner->has_identity_erasure()) {
    auto id = inner->identity_indices();
    if (isr) kramers_flip_deep(id, *isr);
    wrap.set_identity_indices(std::move(id));
  }
  EvalExprNode sentinel{EvalExpr{Constant{1}}};
  return EvalExprNode{std::move(wrap), std::move(inner), std::move(sentinel)};
}

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

  if (opts.kramers_fold_intermediates && opts.kramers_fold_this &&
      (expr->is<Sum>() || expr->is<Product>()))
    if (auto folded = maybe_kramers_fold(expr, uncontract, opts, node_counter))
      return std::move(*folded);

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
