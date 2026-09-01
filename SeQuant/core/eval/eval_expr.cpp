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

namespace {
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
}  // namespace

EvalExpr::EvalExpr(Tensor const& tnsr)
    : op_type_{std::nullopt},
      result_type_{ResultType::Tensor},
      expr_{tnsr.clone()} {
  SEQUANT_ASSERT(!tnsr.indices().empty());
  if (is_tot(tnsr)) {
    // ToT leaf: normalize the spelling channels first, then let the
    // slot-canonicalization (fold ON, the default) report the orientation
    // fold via conjugated_tensors; respell to the canonical orientation so
    // the stored spelling is canonical for flat and ToT leaves alike.
    auto& t0 = expr_->as<Tensor>();
    canon_transform_ = compose(canon_transform_, normalize_leaf_spelling(t0));
    ExprPtrList tlist{expr_};
    auto tn = TensorNetwork(tlist);
    auto md = tn.canonicalize_slots(
        {.cardinal_tensor_labels =
             TensorCanonicalizer::cardinal_tensor_labels()});
    hash_value_ = md.hash_value();
    canon_transform_.phase = md.phase;
    if (!md.conjugated_tensors.empty()) {
      // single-tensor network: the canonical labeling spells this leaf in
      // the swapped orientation -- the fold map is the delta
      canon_transform_ =
          compose(canon_transform_, {.conj = true, .braket_swap = true});
      static_cast<AbstractTensor&>(t0)._swap_bra_ket();
    }
    canon_indices_ = md.get_indices<index_vector>();
    connectivity_ = std::move(md.graph);
  } else {
    // Single (protoindex-free) tensor leaf: normalize to the canonical
    // unmarked spelling; every conjugation channel becomes a CanonTransform
    // byproduct applied on retrieval. Transform bits compose syntactically:
    // marker => conj, slot swap relative to canonical => braket_swap (for a
    // Conjugate tensor the two compose to adjoint, the identity on Hermitian
    // values -- every spelling route lands on one slot + a correct map).
    auto& t = expr_->as<Tensor>();
    canon_transform_ = compose(canon_transform_, normalize_leaf_spelling(t));
    // 3. block-canonicalize WITH the fold (the eval-boundary exception is
    //    gone); a fold performed here toggles the marker, which converts to
    //    transform bits the same way
    auto phase = TensorBlockCanonicalizer{}.apply(t);
    canon_transform_.phase = phase ? -1 : 1;
    if (t.conjugated()) {  // fold byproduct: canonicalize_braket swapped the
      // slots INTO the canonical orientation and marked; convert the marker
      // to transform bits and keep the canonical slots
      canon_transform_ =
          compose(canon_transform_, {.conj = true, .braket_swap = true});
      t.conjugate();
    }
    hash_value_ = hash_terminal_tensor(t);
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

std::int8_t EvalExpr::canon_phase() const noexcept {
  return canon_transform_.phase;
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
/// the leaf's DENOTED symbolic spelling: the stored canonical spelling with
/// the transform re-materialized syntactically -- bra<->ket swapped back when
/// braket_swap is set, the conjugation marker restored when conj is set. TN
/// building and slot counting must see the orientation and conjugation the
/// surrounding expression wrote (the marker also colors the TN graph, which
/// keeps mixed-conj products like C·C^* identity-distinct).
inline ExprPtr denoted_spelling(EvalExpr const& ee) {
  auto t = ee.expr()->as<Tensor>();
  auto const tr = ee.canon_transform();
  if (tr.braket_swap) static_cast<AbstractTensor&>(t)._swap_bra_ket();
  if (tr.conj) t.conjugate();
  return ex<Tensor>(std::move(t));
}

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
    auto e = (!op && node->expr()->is<Tensor>()) ? denoted_spelling(*node)
                                                 : node->expr();
    collect.emplace_back(ExprWithHash{.expr = std::move(e),  //
                                      .hash = salted_hash(node)});
  } else if (node->op_type() == EvalOp::Product && !node.leaf()) {
    collect_tensor_factors(node.left(), collect);
    collect_tensor_factors(node.right(), collect);
  }
}

EvalExprNode binarize(Constant const& c) { return EvalExprNode{EvalExpr{c}}; }

EvalExprNode binarize(Variable const& v) { return EvalExprNode{EvalExpr{v}}; }

EvalExprNode binarize(Power const& p) { return EvalExprNode{EvalExpr{p}}; }

EvalExprNode binarize(Tensor const& t) {
  // Every conjugation channel ('⁺' adjoint label, elementwise-conjugation
  // marker, Conjugate-braket orientation) is normalized by the EvalExpr leaf
  // ctor into the canonical unmarked spelling plus a CanonTransform served
  // on retrieval -- a tensor leaf is always just a leaf.
  return EvalExprNode{EvalExpr{t}};
}

EvalExprNode binarize(Sum const& sum, IndexSet const& uncontract,
                      const BinarizationOptions& opts) {
  using ranges::views::move;
  using ranges::views::transform;
  auto summands = sum.summands()  //
                  | transform([&uncontract, &opts](ExprPtr const& x) {
                      return impl::binarize(x, uncontract, opts);
                    })  //
                  | ranges::to_vector;

  bool const all_tensors =
      ranges::all_of(summands, [](auto&& n) { return n->is_tensor(); });

  [[maybe_unused]] bool const all_scalars =
      ranges::all_of(summands, [](auto&& n) { return n->is_scalar(); });

  SEQUANT_ASSERT(all_tensors | all_scalars);

  auto hvals = summands | transform([](auto&& n) { return salted_hash(n); });

  auto make_sum = [i = 0,                    //
                   hs = imed_hashes(hvals),  //
                   all_tensors, &opts](EvalExpr const& left,
                                       EvalExpr const&) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    if (all_tensors) {
      // partition from the DENOTED orientation (stored canonical slots,
      // re-swapped per the child transform)
      auto const t = denoted_spelling(left)->as<Tensor>();
      return {
          EvalOp::Sum,         //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(opts, bra(t.bra()), ket(t.ket()),
                                            aux(t.aux())),  //
          left.canon_indices(),                             //
          CanonTransform{},                                 //
          h,                                                //
          nullptr};
    } else {
      return {EvalOp::Sum,              //
              ResultType::Scalar,       //
              detail::make_variable(),  //
              {},                       //
              CanonTransform{},         //
              h,                        //
              nullptr};
    }
  };

  return fold_left_to_node(summands | move, make_sum);
}

EvalExprNode binarize(Product const& prod, IndexSet const& uncontract,
                      const BinarizationOptions& opts) {
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

  auto factors =
      prod.factors()  //
      | transform([i = 0, &ltr_uncontr_idxs, &opts](ExprPtr const& x) mutable {
          return impl::binarize(x, ltr_uncontr_idxs.children[i++], opts);
        })  //
      | ranges::to_vector;

  auto hvals = factors | transform([](auto&& n) { return salted_hash(n); });
  auto const hs = imed_hashes(hvals) | ranges::to_vector;

  auto make_prod = [i = 0, &hs, &ltr_uncontr_idxs, &opts](
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
              CanonTransform{},
              h,
              nullptr};
    } else if (left->is_scalar() || right->is_scalar()) {
      // scalar * tensor or tensor * scalar
      auto const& tl = left->is_tensor() ? left : right;
      auto const t =
          denoted_spelling(*tl)->as<Tensor>();  // denoted orientation
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
      bool const hoist_conj =
          !subfacs.empty() &&
          ranges::all_of(subfacs, [](ExprWithHash const& f) {
            return f.expr->is<Tensor>() && f.expr->as<Tensor>().conjugated();
          });
      if (hoist_conj)
        for (auto& f : subfacs) f.expr->as<Tensor>().conjugate();
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
      if (scalar_result) {
        return {EvalOp::Product,                                           //
                ResultType::Scalar,                                        //
                detail::make_variable(),                                   //
                {},                                                        //
                CanonTransform{.phase = canon.phase, .conj = hoist_conj},  //
                h,
                std::move(canon.graph)};
      } else {
        return {EvalOp::Product,     //
                ResultType::Tensor,  //
                detail::make_tensor_wo_symmetries(opts, bra(target_indices.bra),
                                                  ket(target_indices.ket),
                                                  aux(target_indices.aux)),
                canon.get_indices<Index::index_vector>(),  //
                CanonTransform{.phase = canon.phase},      //
                h,
                std::move(canon.graph)};
      }
    }
  };

  if (prod.scalar() == 1) {
    return fold_left_to_node(factors | move, make_prod);
  } else {
    auto left = fold_left_to_node(factors | move, make_prod);
    auto right = binarize(Constant{prod.scalar()});

    auto expr = left->is_tensor()
                    ? detail::make_tensor(denoted_spelling(*left)->as<Tensor>(),
                                          false, opts)
                : left->is_constant() ? (left->expr() * right->expr())
                                      : detail::make_variable();
    auto type = left->is_tensor() ? ResultType::Tensor : ResultType::Scalar;

    auto h = salted_hash(left);
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

namespace impl {

EvalExprNode binarize(ExprPtr const& expr, IndexSet const& uncontract,
                      const BinarizationOptions& opts) {
  if (expr->is<Constant>())  //
    return binarize(expr->as<Constant>());

  if (expr->is<Variable>())  //
    return binarize(expr->as<Variable>());

  if (expr->is<Tensor>())  //
    return binarize(expr->as<Tensor>());

  if (expr->is<Sum>())  //
    return binarize(expr->as<Sum>(), uncontract, opts);

  if (expr->is<Product>())  //
    return binarize(expr->as<Product>(), uncontract, opts);

  if (expr->is<Power>())  //
    return binarize(expr->as<Power>());

  throw Exception("Encountered unsupported expression in binarize.");
}

}  // namespace impl

}  // namespace sequant
