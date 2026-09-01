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
#include <SeQuant/core/utility/exception.hpp>
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
#include <string>
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
                    aux<Aux> a, const TensorSymmetries& syms) {
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
    return ex<Tensor>(label_tensor, bra(), ket(), aux(std::move(indices)),
                      syms);
  } else {
    std::ranges::sort(b);
    std::ranges::sort(k);
    std::ranges::sort(a);

    return ex<Tensor>(label_tensor, std::move(b), std::move(k), std::move(a),
                      syms);
  }
}

template <std::ranges::range Bra, std::ranges::range Ket,
          std::ranges::range Aux>
ExprPtr make_tensor_wo_symmetries(const BinarizationOptions& opts, bra<Bra>&& b,
                                  ket<Ket>&& k, aux<Aux>&& a) {
  return make_tensor<Bra, Ket, Aux>(
      opts, b, k, a,
      TensorSymmetries{.perm = Symmetry::Nonsymm,
                       .braket = BraKetSymmetry::Nonsymm,
                       .column = ColumnSymmetry::Nonsymm});
}

ExprPtr make_tensor(Tensor const& t, bool with_symm,
                    const BinarizationOptions& opts) {
  // carrying the symmetries means carrying the traits: the intermediate's
  // field-dependent symmetries are derived from them and its own slots
  const TensorSymmetries syms =
      with_symm ? t.symmetries()
                : TensorSymmetries{.perm = Symmetry::Nonsymm,
                                   .braket = BraKetSymmetry::Nonsymm,
                                   .column = ColumnSymmetry::Nonsymm};

  return make_tensor(opts, bra(t.bra()), ket(t.ket()), aux(t.aux()), syms);
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
    auto md = tn.canonicalize_slots(
        {.cardinal_tensor_labels =
             TensorCanonicalizer::cardinal_tensor_labels()});
    hash_value_ = md.hash_value();
    canon_transform_.phase = md.phase;
    canon_indices_ = md.get_indices<index_vector>();
    connectivity_ = std::move(md.graph);
  } else {
    // Single (protoindex-free) tensor: block-canonicalize it in place. This is
    // a lightweight per-tensor canonicalization (no deep tensor-network
    // canonicalization is needed for a tensor that is not itself a network),
    // and it normalizes bra<->ket orientation for braket-symmetric tensors so
    // that equivalent half-tensor forms (e.g. X{a;;x} and X{;a;x}) fold.
    auto& t = expr_->as<Tensor>();
    // A leaf keeps its as-written orientation: only a Symm bra<->ket
    // exchange, a free respelling, folds; a signed exchange (AntiSymm) does
    // not, and a Conjugate tensor's bundles are never exchanged. The yielder
    // is asked for the leaf's canonical spelling and the engine takes its
    // array as the leaf's own value, so a leaf's phase is a
    // cache-orientation round trip that never reaches the value, which is
    // why a respelling that costs a sign is not folded here.
    auto phase =
        TensorBlockCanonicalizer{/*fold_signed_braket=*/false}.apply(t);
    canon_transform_.phase = phase ? -1 : 1;
    // The leaf hash (hash_terminal_tensor) keys the array by label, slot
    // layout, and states: a K-conjugated leaf over a complex basis is its
    // own array and stays a leaf, so the state must separate it from its
    // unmarked twin. An Adjoint IR node carries the marked spelling as its
    // expr() but hashes the bare leaf salted by the op (make_adjoint_node).
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
  // The leaf hash keys an array by its label, slot layout, and states; the
  // symmetry traits are not part of it. Both states are part of the leaf's
  // value identity: a marked spelling is a different array from its unmarked
  // twin with the same slots, so the two must not share a cache slot. The
  // default state adds nothing, so unmarked leaves keep their hash.
  if (tnsr.adjointed() || tnsr.kconjugated())
    hash::combine(h, static_cast<std::uint8_t>((tnsr.adjointed() ? 1 : 0) |
                                               (tnsr.kconjugated() ? 2 : 0)));
  // Over a real basis an Odd-parity array is imaginary while the Even one is
  // real, so the two must not share a cache slot. Even and None add nothing,
  // so the default keys stay.
  if (tnsr.conjugation_symmetry() == ConjugationSymmetry::AntiSymm)
    hash::combine(h, static_cast<std::uint8_t>(tnsr.conjugation_symmetry()));
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

/// Assembles the Adjoint IR node binarize(Tensor) serves a state with:
/// Adjoint(@p adjointed as the expr, Constant{1} sentinel) over @p bare_leaf,
/// laid out by @p canon_ix (the backends' adjoint kernel permutes the operand
/// from its own annotation to this one, then conjugates elementwise), with
/// the node hash = the bare-leaf hash salted by EvalOp::Adjoint so cache
/// lookups don't collide.
EvalExprNode make_adjoint_node(EvalExprNode bare_leaf, ExprPtr adjointed,
                               EvalExpr::index_vector canon_ix,
                               CanonTransform transform) {
  EvalExprNode sentinel{EvalExpr{Constant{1}}};
  auto h = bare_leaf->hash_value();
  hash::combine(h, static_cast<size_t>(EvalOp::Adjoint));
  EvalExpr adj{EvalOp::Adjoint,
               ResultType::Tensor,
               std::move(adjointed),
               std::move(canon_ix),
               transform,
               h,
               nullptr};
  return EvalExprNode{std::move(adj), std::move(bare_leaf),
                      std::move(sentinel)};
}

EvalExprNode binarize(Constant const& c) { return EvalExprNode{EvalExpr{c}}; }

EvalExprNode binarize(Variable const& v) { return EvalExprNode{EvalExpr{v}}; }

EvalExprNode binarize(Power const& p) { return EvalExprNode{EvalExpr{p}}; }

EvalExprNode binarize(Tensor const& t, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter) {
  // A leaf keeps its as-written orientation; the states are served as IR
  // ops over the bare leaf. A marked tensor never carries a sign (signs live
  // in scalars), so clearing a state consumes none.
  if (!t.adjointed() && !t.kconjugated()) return EvalExprNode{EvalExpr{t}};
  if (t.adjointed()) {
    // T⁺: Adjoint (bra/ket permute, then conjugate) over the bare array,
    // which is the adjoint of the adjointed spelling: Tensor::adjoint()
    // exchanges the bundles back and clears the state, and the normalization
    // it runs re-examines only the K state, already normalized, so the sign
    // is 1. That K state, if any, stays on the operand leaf, its own array.
    // The '⁺' is served the same way in every basis: an index-mutating API
    // can leave a '⁺' over a real basis, and the Adjoint node is
    // value-correct there too. IR shape: Adjoint(Tensor{<bare>},
    // Constant{1}); the Constant(1) right child is a sentinel so the
    // FullBinaryNode invariant ("every non-leaf has two children") holds.
    Tensor bare{t};
    [[maybe_unused]] const auto sign = bare.adjoint();
    SEQUANT_ASSERT(sign == 1);
    SEQUANT_ASSERT(!bare.adjointed());
    SEQUANT_ASSERT(bare.kconjugated() == t.kconjugated());
    return make_adjoint_node(
        binarize(bare, uncontract, opts, node_counter), ex<Tensor>(t),
        t.const_indices() | ranges::to<EvalExpr::index_vector>,
        CanonTransform{});
  }
  // T꙳ over a real basis is the elementwise conjugate of the bare array: the
  // Adjoint kernel with an identity layout (permute-then-conjugate with the
  // operand's own annotation) is that conjugation. Over a complex basis the
  // K-conjugated operator's matrix is an array of its own, a plain leaf
  // whose hash carries the state.
  if (t.base_field() == Field::Real) {
    Tensor bare{t};
    [[maybe_unused]] const auto sign = bare.set_states(false, false);
    SEQUANT_ASSERT(sign == 1);
    EvalExprNode leaf{EvalExpr{bare}};
    auto canon_ix = leaf->canon_indices();
    return make_adjoint_node(std::move(leaf), ex<Tensor>(t),
                             std::move(canon_ix), CanonTransform{});
  }
  return EvalExprNode{EvalExpr{t}};
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

  // Every binary Sum produced by fold_left_to_node below folds the running
  // accumulator (the chain seed, or a prior chain Sum) in as the left
  // operand (see fold_left_to_node in binary_node.hpp: the accumulator is
  // always `l`), so every chain Sum node accumulates its left operand in
  // place -- see EvalExpr::accumulate_in_place.
  auto make_sum = [i = 0,                                        //
                   hs = imed_hashes(hvals) | ranges::to_vector,  //
                   all_tensors, &opts](EvalExpr const& left,
                                       EvalExpr const&) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    if (all_tensors) {
      // This node takes its slot layout from the left operand, and nothing
      // else: a marked operand carries no sign (an Adjoint node, or a
      // K-conjugated leaf that is its own array, has phase 1; a sign a state
      // consumed lives in a scalar), so every operand here is spelled by its
      // value.
      auto const& t = left.as_tensor();
      EvalExpr result{
          EvalOp::Sum,         //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(opts, bra(t.bra()), ket(t.ket()),
                                            aux(t.aux())),  //
          left.canon_indices(),                             //
          CanonTransform{},                                 //
          h,                                                //
          nullptr};
      result.set_accumulate_in_place(true);
      return result;
    } else {
      EvalExpr result{EvalOp::Sum,              //
                      ResultType::Scalar,       //
                      detail::make_variable(),  //
                      {},                       //
                      CanonTransform{},         //
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
              CanonTransform{},
              h,
              nullptr};
    } else if (left->is_scalar() || right->is_scalar()) {
      // scalar * tensor or tensor * scalar
      auto const& tl = left->is_tensor() ? left : right;
      // the tensor operand supplies this node's slot layout, and nothing
      // else: a marked operand carries no sign (an Adjoint node, or a
      // K-conjugated leaf that is its own array, has phase 1; a sign a state
      // consumed lives in a scalar such as the one beside it here)
      auto const& t = tl->as_tensor();
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
      auto ts = subfacs | transform([](auto&& t) { return t.expr; });
      IndexGroups<IndexVec> const target_indices = [&ts, &uncontracted_idxs]() {
        // route each surviving hyperindex to its correct slot
        // (bra, ket, or aux) based on which slot it occupies in
        // the factor tensors .. if appears in multiple slots put into aux
        auto counts =
            get_used_indices_with_counts(ex<Product>(ts | ranges::to_vector));
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
      EvalExpr result =
          scalar_result
              ? EvalExpr{EvalOp::Product,          //
                         ResultType::Scalar,       //
                         detail::make_variable(),  //
                         {},                       //
                         CanonTransform{.phase = canon.phase},  //
                         h,
                         std::move(canon.graph)}
              : EvalExpr{EvalOp::Product,     //
                         ResultType::Tensor,  //
                         detail::make_tensor_wo_symmetries(
                             opts, bra(target_indices.bra),
                             ket(target_indices.ket), aux(target_indices.aux)),
                         canon.get_indices<Index::index_vector>(),  //
                         CanonTransform{.phase = canon.phase},      //
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

    ExprPtr expr;
    if (left->is_tensor()) {
      // the operand supplies the layout, and nothing else: a marked factor
      // carries no sign of its own, so the product scalar this node applies
      // is the only sign here
      expr = detail::make_tensor(left->as_tensor(), false, opts);
    } else if (left->is_constant()) {
      expr = left->expr() * right->expr();
    } else {
      expr = detail::make_variable();
    }
    auto type = left->is_tensor() ? ResultType::Tensor : ResultType::Scalar;

    auto h = left->hash_value();
    hash::combine(h, right->hash_value());
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
                      const BinarizationOptions& opts,
                      std::size_t& node_counter) {
  if (expr->is<Constant>())  //
    return binarize(expr->as<Constant>());

  if (expr->is<Variable>())  //
    return binarize(expr->as<Variable>());

  if (expr->is<Tensor>())  //
    return binarize(expr->as<Tensor>(), uncontract, opts, node_counter);

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
