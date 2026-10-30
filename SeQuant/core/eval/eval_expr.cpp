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
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>
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
#include <cstddef>
#include <optional>
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
        get_default_context_snapshot().cardinal_tensor_labels());
    hash_value_ = md.hash_value();
    canon_phase_ = md.phase;
    canon_indices_ = md.get_indices<index_vector>();
    connectivity_ = std::move(md.graph);
  } else {
    // Single (protoindex-free) tensor: block-canonicalize it in place. This is
    // a lightweight per-tensor canonicalization (no deep tensor-network
    // canonicalization is needed for a tensor that is not itself a network),
    // and it normalizes bra<->ket orientation for braket-symmetric tensors so
    // that equivalent half-tensor forms (e.g. X{a;;x} and X{;a;x}) fold.
    auto& t = expr_->as<Tensor>();
    auto phase = TensorBlockCanonicalizer{}.apply(t);
    canon_phase_ = phase ? -1 : 1;
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

bool EvalExpr::is_elementwise_power() const noexcept {
  return op_type() == EvalOp::Power;
}

bool is_unary(EvalOp op) noexcept {
  return op == EvalOp::Adjoint || op == EvalOp::Power;
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
/// @note An Index object's IndexSpace type and quantum numbers, and its basis
///       instance and name (if any), contribute to the hash.
///
template <typename T>
size_t hash_indices(T const& indices) noexcept {
  size_t h = 0;
  for (auto const& idx : indices) {
    hash::combine(h, hash::value(idx.space().type().to_int32()));
    hash::combine(h, hash::value(idx.space().qns().to_int32()));
    if (idx.basis().has_basis_instance()) {
      hash::combine(h, hash::value(*idx.basis().basis_instance()));
      if (idx.basis().has_name())
        hash::combine(h, hash::value(idx.basis().name()));
    }
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
    for (auto&& ix : expr->as<Tensor>().const_indices()) {
      result.emplace(ix);
      // proto-only indices are still live outer modes of the ToT array
      for (auto&& p : ix.proto_indices()) result.emplace(p);
    }
  else if (expr->is<TensorPower>())
    all_indices(result, expr->as<TensorPower>().base());
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
      node->is_tensor() && (!op || *op == EvalOp::Sum || is_unary(*op)))
    // Treat a unary node (Adjoint, Power) the same as Sum here: it produces a
    // tensor result that enters a parent Product as a single factor — the
    // parent shouldn't recurse past its boundary, just collect the tensor
    // held in node->expr() and move on.
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

/// Assembles a unary node (see is_unary) of type \p op over \p operand, with
/// a Constant \p sentinel as the right child. The hash is the operand's
/// combined with \p op and, for a Power, the exponent.
EvalExprNode make_unary_node(EvalOp op, EvalExprNode operand,
                             Constant const& sentinel, ExprPtr expr,
                             EvalExpr::index_vector canon_indices,
                             std::int8_t phase) {
  SEQUANT_ASSERT(is_unary(op));
  auto h = operand->hash_value();
  hash::combine(h, static_cast<size_t>(op));
  if (op == EvalOp::Power) {
    SEQUANT_ASSERT(sentinel.value().imag() == 0);
    hash::combine(h, hash::value(sentinel.value().real()));
  }
  EvalExpr node{op,
                ResultType::Tensor,
                std::move(expr),
                std::move(canon_indices),
                phase,
                h,
                nullptr};
  return EvalExprNode{std::move(node), std::move(operand),
                      EvalExprNode{EvalExpr{sentinel}}};
}

}  // namespace

EvalExprNode binarize(Tensor const& t) {
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
  if (is_adjoint_label(t.label())) {
    // The Adjoint node carries the *adjointed* tensor (so its canon_indices
    // reflect the slot order parents see).

    // Build the bare-label operand: copy and call adjoint() to toggle the
    // marker off and swap bra/ket back to natural orientation.
    Tensor bare{t};
    bare.adjoint();
    SEQUANT_ASSERT(!is_adjoint_label(bare.label()));

    // Sentinel right child is Constant(1); the hash differs from the bare
    // leaf's so cache lookups don't collide.
    return make_unary_node(EvalOp::Adjoint, EvalExprNode{EvalExpr{bare}},
                           Constant{1}, t.clone(),
                           t.indices() | ranges::to<EvalExpr::index_vector>, 1);
  }
  return EvalExprNode{EvalExpr{t}};
}

EvalExprNode binarize(TensorPower const& power, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter,
                      detail::InlineContext* inline_ctx) {
  auto operand =
      impl::binarize(power.base(), uncontract, opts, node_counter, inline_ctx);
  auto const& t = operand->as_tensor();
  auto expr = detail::make_tensor_wo_symmetries(opts, bra(t.bra()),
                                                ket(t.ket()), aux(t.aux()));
  auto canon_indices = operand->canon_indices();
  // base = phase * canonical operand, so base^e = phase^e * operand^e; a
  // non-integer exponent never meets phase -1 (pow() rejects antisymmetric
  // bases there)
  auto const& e = power.exponent();
  bool const even = denominator(e) == 1 && numerator(e) % 2 == 0;
  std::int8_t const phase = even ? 1 : operand->canon_phase();
  return make_unary_node(EvalOp::Power, std::move(operand),
                         Constant{power.exponent()}, std::move(expr),
                         std::move(canon_indices), phase);
}

EvalExprNode binarize(Sum const& sum, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter,
                      detail::InlineContext* inline_ctx) {
  using ranges::views::move;
  using ranges::views::transform;
  auto summands =
      sum.summands()  //
      | transform([&uncontract, &opts, &node_counter,
                   inline_ctx](ExprPtr const& x) {
          return impl::binarize(x, uncontract, opts, node_counter, inline_ctx);
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
      auto const& t = left.as_tensor();
      EvalExpr result{
          EvalOp::Sum,         //
          ResultType::Tensor,  //
          detail::make_tensor_wo_symmetries(opts, bra(t.bra()), ket(t.ket()),
                                            aux(t.aux())),  //
          left.canon_indices(),                             //
          1,                                                //
          h,                                                //
          nullptr};
      result.set_accumulate_in_place(true);
      return result;
    } else {
      EvalExpr result{EvalOp::Sum,              //
                      ResultType::Scalar,       //
                      detail::make_variable(),  //
                      {},                       //
                      1,                        //
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
                      std::size_t& node_counter,
                      detail::InlineContext* inline_ctx) {
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
                 | transform([i = 0, &ltr_uncontr_idxs, &opts, &node_counter,
                              inline_ctx](ExprPtr const& x) mutable {
                     return impl::binarize(x, ltr_uncontr_idxs.children[i++],
                                           opts, node_counter, inline_ctx);
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
      auto const& t = tl->as_tensor();
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
      IndexGroups<IndexVec> const target_indices = [prod = ex<Product>(ts),
                                                    &uncontracted_idxs]() {
        // route each surviving hyperindex to its correct slot
        // (bra, ket, or aux) based on which slot it occupies in
        // the factor tensors .. if appears in multiple slots put into aux
        auto counts = get_used_indices_with_counts(prod);
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
          get_default_context_snapshot().cardinal_tensor_labels(),
          &named_indices);
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
                    ? detail::make_tensor(left->as_tensor(), false, opts)
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

namespace detail {

namespace {

/// \return true if \p x and \p y have the same bra, ket and aux ranks and the
///         same index space in every slot
template <typename X, typename Y>
bool same_slot_spaces(X const& x, Y const& y) {
  auto same = [](auto const& a, auto const& b) {
    return std::ranges::equal(a, b, {}, &Index::space, &Index::space);
  };
  return same(x.bra(), y.bra()) && same(x.ket(), y.ket()) &&
         same(x.aux(), y.aux());
}

/// \return true if two heads cannot be told apart by the leaves naming them
bool same_shape(ResultExpr const& x, ResultExpr const& y) {
  return x.produces_tensor() == y.produces_tensor() && same_slot_spaces(x, y);
}

}  // namespace

/// The definitions of a system of equations and the state of inlining them;
/// see binarize(std::span<const ResultExpr>, const BinarizationOptions&).
struct InlineContext {
  /// \param system the equations
  /// \param consumed receives true for every equation inlined somewhere
  InlineContext(std::span<const ResultExpr> system,
                container::svector<bool>& consumed)
      : system_(system), consumed_(consumed), fresh_([this](Index const& ix) {
          return !used_labels_.contains(std::wstring(ix.label()));
        }) {
    for (std::size_t k = 0; k != system_.size(); ++k) {
      ResultExpr const& res = system_[k];
      if (!res.has_label())
        throw Exception("inline: every equation needs a labelled head");
      auto& same_label = by_label_[res.label()];
      for (std::size_t j : same_label)
        if (same_shape(system_[j], res))
          throw Exception("inline: duplicate definition of " +
                          toUtf8(res.label()));
      same_label.push_back(k);
      heads_.push_back(res.produces_tensor() ? res.result_as_tensor()
                                             : Tensor{});

      auto use = [this](Index const& ix) {
        used_labels_.emplace(ix.label());
        for (auto const& p : ix.proto_indices())
          used_labels_.emplace(p.label());
      };
      for (auto const& ix : res.indices()) use(ix);
      for (auto const& ix : get_used_indices(res.expression())) use(ix);
    }
  }

  InlineContext(InlineContext const&) = delete;
  InlineContext& operator=(InlineContext const&) = delete;

  /// \return the equation whose tensor head \p leaf names, if any
  /// \throw Exception if \p leaf matches a head but not its symmetries
  std::optional<std::size_t> match(Tensor const& leaf) const {
    auto const it = by_label_.find(std::wstring(leaf.label()));
    if (it == by_label_.end()) return std::nullopt;
    for (std::size_t k : it->second) {
      if (!system_[k].produces_tensor()) continue;
      Tensor const& head = heads_[k];
      if (!same_slot_spaces(leaf, head)) continue;
      if (leaf.symmetry() != head.symmetry() ||
          leaf.braket_symmetry() != head.braket_symmetry() ||
          leaf.column_symmetry() != head.column_symmetry())
        throw Exception(
            "inline: leaf " + toUtf8(io::serialization::to_string(leaf)) +
            " names head " + toUtf8(io::serialization::to_string(head)) +
            " but declares different symmetries; use result_as_tensor()");
      return k;
    }
    return std::nullopt;
  }

  /// \return the equation whose scalar head \p leaf names, if any
  /// \throw Exception if \p leaf is conjugated and names a scalar head
  std::optional<std::size_t> match(Variable const& leaf) const {
    auto const it = by_label_.find(std::wstring(leaf.label()));
    if (it == by_label_.end()) return std::nullopt;
    for (std::size_t k : it->second) {
      if (system_[k].produces_tensor()) continue;
      if (leaf.conjugated())
        throw Exception("inline: conjugated use of head " +
                        toUtf8(leaf.label()) +
                        " is not supported; define its conjugate as its own "
                        "equation");
      return k;
    }
    return std::nullopt;
  }

  /// \return the tree of equation \p k, with every definition inlined
  EvalExprNode binarize_equation(std::size_t k, BinarizationOptions const& opts,
                                 std::size_t& node_counter) {
    return binarize_guarded(k, system_[k], opts, node_counter);
  }

  /// \return the tree of equation \p k as a factor in place of \p leaf (null
  ///         for a Variable leaf): its head's indices renamed to \p leaf's and
  ///         every other index of its rhs to a fresh one; not stamped with
  ///         \p opts.node_batch_axes
  EvalExprNode inline_definition(std::size_t k, Tensor const* leaf,
                                 BinarizationOptions const& opts) {
    consumed_[k] = true;
    ResultExpr const& def = system_[k];

    container::map<Index, Index> renaming;
    Tensor head;
    if (leaf) {
      head = heads_[k];
      auto map_slots = [&renaming](auto const& from, auto const& to) {
        SEQUANT_ASSERT(std::ranges::size(from) == std::ranges::size(to));
        auto t = std::ranges::begin(to);
        for (auto const& ix : from) {
          if (!ix.nonnull()) {
            ++t;
            continue;
          }
          auto const [it, inserted] = renaming.emplace(ix, *t);
          if (!inserted && it->second != *t)
            throw Exception("inline: a leaf names head index " +
                            toUtf8(ix.full_label()) + " by two indices");
          ++t;
        }
      };
      map_slots(head.bra(), leaf->bra());
      map_slots(head.ket(), leaf->ket());
      map_slots(head.aux(), leaf->aux());
      transform_indices(head, renaming);
      reset_tags(head);
    }
    // every other index of the rhs (a dummy, or every index under a scalar
    // head) gets a fresh one, so the subtree cannot collide with the
    // consumer's or another inlined subtree's indices
    container::svector<Index> dummies;
    for (auto const& ix : get_used_indices(def.expression()))
      if (!renaming.contains(ix)) dummies.push_back(ix);
    for (auto const& ix : dummies) {
      Index f = fresh_.make(ix.drop_proto_indices());
      used_labels_.emplace(f.label());
      renaming.emplace(ix, std::move(f));
    }
    // the proto indices of a renamed dummy follow the renaming
    for (auto const& ix : dummies) {
      if (!ix.has_proto_indices()) continue;
      container::vector<Index> protos;
      for (auto const& p : ix.proto_indices()) {
        auto const it = renaming.find(p);
        protos.push_back(it == renaming.end() ? p : it->second);
      }
      Index& f = renaming.at(ix);
      f = Index(Index(f), std::move(protos), ix.symmetric_proto_indices());
    }

    ExprPtr rhs = transform_expr(def.expression(), renaming);
    ResultExpr const renamed =
        leaf ? ResultExpr(head, std::move(rhs))
             : ResultExpr(def.result_as_variable(), std::move(rhs));
    BinarizationOptions sub_opts = opts;
    sub_opts.node_batch_axes.clear();
    std::size_t node_counter = 0;
    return binarize_guarded(k, renamed, sub_opts, node_counter);
  }

 private:
  /// binarizes \p res, the (renamed) equation \p k, with \p k marked active
  /// \throw Exception if \p k is active already
  EvalExprNode binarize_guarded(std::size_t k, ResultExpr const& res,
                                BinarizationOptions const& opts,
                                std::size_t& node_counter) {
    if (!active_.emplace(k).second)
      throw Exception("inline: cyclic definition of " +
                      toUtf8(system_[k].label()));
    auto tree = impl::binarize(res, opts, node_counter, this);
    active_.erase(k);
    return tree;
  }

  std::span<const ResultExpr> system_;
  container::svector<bool>& consumed_;
  /// head label -> the equations with that label
  container::map<std::wstring, container::svector<std::size_t>> by_label_;
  /// the head of every equation as a Tensor (empty for a scalar head)
  container::svector<Tensor> heads_;
  /// the equations whose trees are being built
  container::set<std::size_t> active_;
  /// the labels of every index of the system and of every fresh index
  container::set<std::wstring> used_labels_;
  IndexFactory fresh_;
};

}  // namespace detail

namespace impl {

EvalExprNode binarize(ExprPtr const& expr, IndexSet const& uncontract,
                      const BinarizationOptions& opts,
                      std::size_t& node_counter,
                      detail::InlineContext* inline_ctx) {
  if (expr->is<Constant>())  //
    return binarize(expr->as<Constant>());

  if (expr->is<Variable>()) {
    auto const& v = expr->as<Variable>();
    if (inline_ctx)
      if (auto const k = inline_ctx->match(v))
        return inline_ctx->inline_definition(*k, nullptr, opts);
    return binarize(v);
  }

  if (expr->is<Tensor>()) {
    auto const& t = expr->as<Tensor>();
    if (inline_ctx)
      if (auto const k = inline_ctx->match(t))
        return inline_ctx->inline_definition(*k, &t, opts);
    return binarize(t);
  }

  if (expr->is<Sum>())  //
    return binarize(expr->as<Sum>(), uncontract, opts, node_counter,
                    inline_ctx);

  if (expr->is<Product>())  //
    return binarize(expr->as<Product>(), uncontract, opts, node_counter,
                    inline_ctx);

  if (expr->is<Power>())  //
    return binarize(expr->as<Power>());

  if (expr->is<TensorPower>())  //
    return binarize(expr->as<TensorPower>(), uncontract, opts, node_counter,
                    inline_ctx);

  throw Exception("Encountered unsupported expression in binarize.");
}

EvalExprNode binarize(ResultExpr const& res, const BinarizationOptions& opts,
                      std::size_t& node_counter,
                      detail::InlineContext* inline_ctx) {
  ResultExpr const lowered = lower_broadcast_sum(res);

  // collect result indices so they are not contracted in the expression tree
  IndexSet uncontract;
  for (auto&& ix : lowered.indices()) uncontract.emplace(ix);

  // the root's bra/ket would be positional (see the deprecated
  // binarize(ExprPtr)), but it is overwritten below with
  // lowered.result_as_tensor(), so the layout is caller-determined
  EvalExprNode tree = binarize(lowered.expression(), uncontract, opts,
                               node_counter, inline_ctx);

  const bool is_scalar =
      lowered.bra().empty() && lowered.ket().empty() && lowered.aux().empty();

  if (tree.size() < 2) {
    // We want to have a result node with the result from the ResultExpr.
    // In order for that to work, we need a dedicated result node in the first
    // place. Hence, we adapt the represented expression for terminals to be
    // that terminal multiplied by 1.
    EvalExpr result = is_scalar
                          ? *sequant::binarize(lowered.result_as_variable())
                          : *sequant::binarize(lowered.result_as_tensor());
    EvalOpSetter{}.set(result, EvalOp::Product);

    tree = EvalExprNode(std::move(result), std::move(tree),
                        sequant::binarize(Constant(1)));
  }
  SEQUANT_ASSERT(tree.size() > 1);

  if (is_scalar) {
    if (lowered.has_label()) {
      tree->expr()->as<Variable>().set_label(lowered.label());
    }
  } else {
    tree->expr()->as<Tensor>() = lowered.result_as_tensor();
  }

  return tree;
}

container::svector<EvalExprNode> binarize_system(
    std::span<const ResultExpr> system, const BinarizationOptions& opts,
    container::svector<bool>& consumed) {
  consumed.assign(system.size(), false);
  detail::InlineContext ctx(system, consumed);
  container::svector<EvalExprNode> trees;
  for (std::size_t k = 0; k != system.size(); ++k) {
    std::size_t node_counter = 0;
    trees.emplace_back(ctx.binarize_equation(k, opts, node_counter));
    // see binarize(ExprPtr const&, ...)
    if (!opts.node_batch_axes.empty()) {
      SEQUANT_ASSERT(node_counter == opts.node_batch_axes.size());
    }
  }
  return trees;
}

}  // namespace impl

ResultExpr lower_broadcast_sum(ResultExpr const& res) {
  if (!res.expression()->is<Sum>()) return res;

  const auto& summands = res.expression()->as<Sum>().summands();
  // created at the first summand that needs broadcasting
  std::optional<Sum::summands_type> lowered;
  for (std::size_t s = 0; s != summands.size(); ++s) {
    const ExprPtr& summand = summands[s];
    const container::set<Index> supplied = detail::supplied_result_indices(
        res, *detail::without_symmetrizers(summand));

    std::optional<Product> broadcast;
    auto append_ones = [&](const auto& head_indices, SlotType slot) {
      for (const Index& idx : head_indices) {
        if (!idx.nonnull() || supplied.contains(idx)) continue;
        if (slot == SlotType::Aux) {
          throw Exception(
              "lower_broadcast_sum: cannot broadcast over aux index " +
              toUtf8(idx.full_label()));
        }
        if (!broadcast) {
          broadcast.emplace();
          broadcast->append(1, summand, Product::Flatten::Once);
        }
        broadcast->append(1, make_ones(idx, slot), Product::Flatten::No);
      }
    };
    append_ones(res.bra(), SlotType::Bra);
    append_ones(res.ket(), SlotType::Ket);
    append_ones(res.aux(), SlotType::Aux);

    if (broadcast && !lowered) {
      lowered.emplace(summands.begin(), summands.begin() + s);
    }
    if (lowered) {
      lowered->push_back(broadcast ? ex<Product>(std::move(*broadcast))
                                   : summand);
    }
  }
  if (!lowered) return res;

  ResultExpr result = res;
  result = ex<Sum>(std::move(*lowered), Sum::move_only_tag{});
  return result;
}

}  // namespace sequant
