#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/wstring.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <range/v3/action.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/functional.hpp>
#include <range/v3/iterator.hpp>
#include <range/v3/view.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <memory>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace sequant {

namespace {

size_t hash_terminal_tensor(Tensor const&) noexcept;

bool is_tot(Tensor const& t) noexcept {
  return ranges::any_of(t.const_indices(), &Index::has_proto_indices);
}

}  // namespace

std::string to_label_annotation(const Index& idx) {
  using namespace ranges::views;
  using ranges::to;

  return sequant::to_string(idx.label()) +
         (idx.proto_indices() | transform(&Index::label) |
          transform([](auto&& str) { return sequant::to_string(str); }) |
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
  if (is_tot(tnsr)) {
    ExprPtrList tlist{expr_};
    auto tn = TensorNetworkV2(tlist);
    auto md =
        tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels());
    hash_value_ = md.hash_value();
    canon_phase_ = md.phase;
    canon_indices_ = md.get_indices<index_vector>();
  } else {
    hash_value_ = hash_terminal_tensor(tnsr);
    canon_phase_ = 1;
    canon_indices_ = tnsr.indices() | ranges::to<index_vector>;
  }
}

EvalExpr::EvalExpr(Constant const& c)
    : op_type_{std::nullopt},
      result_type_{ResultType::Scalar},
      hash_value_{hash::value(c)},
      expr_{c.clone()} {}

EvalExpr::EvalExpr(Variable const& v)
    : op_type_{std::nullopt},
      result_type_{ResultType::Scalar},
      hash_value_{hash::value(v)},
      expr_{v.clone()} {}

EvalExpr::EvalExpr(EvalOp op, ResultType res, ExprPtr const& ex,
                   index_vector ixs, std::int8_t p, size_t h)
    : op_type_{op},
      result_type_{res},
      expr_{ex.clone()},
      canon_indices_{std::move(ixs)},
      canon_phase_{p},
      hash_value_{h} {}

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

bool EvalExpr::is_primary() const noexcept { return !op_type(); }

bool EvalExpr::is_sum() const noexcept { return op_type() == EvalOp::Sum; }

bool EvalExpr::is_product() const noexcept {
  return op_type() == EvalOp::Product;
}

Tensor const& EvalExpr::as_tensor() const { return expr().as<Tensor>(); }

Constant const& EvalExpr::as_constant() const { return expr().as<Constant>(); }

Variable const& EvalExpr::as_variable() const { return expr().as<Variable>(); }

std::string EvalExpr::label() const noexcept {
  if (is_tensor())
    return to_string(as_tensor().label()) + "(" + indices_annot() + ")";
  else if (is_constant()) {
    return sequant::to_string(sequant::deparse(as_constant()));
  } else {
    assert(is_variable());
    return to_string(as_variable().label());
  }
}

std::int8_t EvalExpr::canon_phase() const noexcept { return canon_phase_; }

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
  hash::combine(h, hash_indices(tnsr.const_indices()));
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
  using ranges::views::transform;
  return inits(rng) | transform([](auto&& v) {
           return hash::range_unordered(ranges::begin(v), ranges::end(v));
         });
}

struct ExprWithHash {
  ExprPtr expr;
  size_t hash;
};

namespace dummy {
inline constexpr std::wstring_view label_tensor{L"I"};
inline constexpr std::wstring_view label_scalar{L"Z"};

template <typename... Args>
ExprPtr make_tensor(Args&&... args) {
  return ex<Tensor>(label_tensor, std::forward<Args>(args)...);
}

ExprPtr make_tensor(Tensor const& t, bool with_symm) {
  return with_symm ? ex<Tensor>(label_tensor,           //
                                bra(t.bra()),           //
                                ket(t.ket()),           //
                                aux(t.aux()),           //
                                t.symmetry(),           //
                                t.braket_symmetry(),    //
                                t.particle_symmetry())  //
                   : ex<Tensor>(label_tensor,           //
                                bra(t.bra()),           //
                                ket(t.ket()),           //
                                aux(t.aux()));
}

ExprPtr make_variable() { return ex<Variable>(label_scalar); }

}  // namespace dummy

using EvalExprNode = FullBinaryNode<EvalExpr>;

///
/// \brief Collect tensors appearing as a factor at the leaf node of a product
///        sub-tree, or, at the root node of a sum sub-tree.
///
template <typename Rng>
void collect_tensor_factors(EvalExprNode const& node,  //
                            Rng& collect) {
  static_assert(std::is_same_v<ranges::range_value_t<Rng>, ExprWithHash>);

  if (auto op = node->op_type();
      node->is_tensor() && (!op || *op == EvalOp::Sum))
    collect.emplace_back(ExprWithHash{node->expr(), node->hash_value()});
  else if (node->op_type() == EvalOp::Product && !node.leaf()) {
    collect_tensor_factors(node.left(), collect);
    collect_tensor_factors(node.right(), collect);
  }
}

EvalExprNode binarize(Constant const& c) { return EvalExprNode{EvalExpr{c}}; }

EvalExprNode binarize(Variable const& v) { return EvalExprNode{EvalExpr{v}}; }

EvalExprNode binarize(Tensor const& t) { return EvalExprNode{EvalExpr{t}}; }

EvalExprNode binarize(Sum const& sum) {
  using ranges::views::move;
  using ranges::views::transform;
  auto summands = sum.summands()                                             //
                  | transform([](ExprPtr const& x) { return binarize(x); })  //
                  | ranges::to_vector;

  bool const all_tensors =
      ranges::all_of(summands, [](auto&& n) { return n->is_tensor(); });

  bool const all_scalars =
      ranges::all_of(summands, [](auto&& n) { return n->is_scalar(); });

  assert(all_tensors | all_scalars);

  auto hvals = summands | transform([](auto&& n) { return n->hash_value(); });

  auto make_sum = [i = 0,                    //
                   hs = imed_hashes(hvals),  //
                   all_tensors](EvalExpr const& left,
                                EvalExpr const& right) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    if (all_tensors) {
      auto const& t = left.as_tensor();
      return {EvalOp::Sum,                                                   //
              ResultType::Tensor,                                            //
              dummy::make_tensor(bra(t.bra()), ket(t.ket()), aux(t.aux())),  //
              left.canon_indices(),                                          //
              1,                                                             //
              h};
    } else {
      return {EvalOp::Sum,             //
              ResultType::Scalar,      //
              dummy::make_variable(),  //
              {},                      //
              1,                       //
              h};
    }
  };

  return fold_left_to_node(summands | move, make_sum);
}

EvalExprNode binarize(Product const& prod) {
  using ranges::views::move;
  using ranges::views::transform;
  auto factors = prod.factors()                                             //
                 | transform([](ExprPtr const& x) { return binarize(x); })  //
                 | ranges::to_vector;

  auto hvals = factors | transform([](auto&& n) { return n->hash_value(); });
  auto const hs = imed_hashes(hvals) | ranges::to_vector;

  auto make_prod = [i = 0, &hs](EvalExprNode const& left,
                                EvalExprNode const& right) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    if (left->is_scalar() && right->is_scalar()) {
      // scalar * scalar
      return {EvalOp::Product,
              ResultType::Scalar,
              dummy::make_variable(),
              {},
              1,
              h};
    } else if (left->is_scalar() || right->is_scalar()) {
      // scalar * tensor or tensor * scalar
      auto const& tl = left->is_tensor() ? left : right;
      auto const& t = tl->as_tensor();
      return {EvalOp::Product,                                               //
              ResultType::Tensor,                                            //
              dummy::make_tensor(bra(t.bra()), ket(t.ket()), aux(t.aux())),  //
              tl->canon_indices(),                                           //
              1,                                                             //
              h};
    } else {
      // tensor * tensor
      container::svector<ExprWithHash> subfacs;
      collect_tensor_factors(left, subfacs);
      collect_tensor_factors(right, subfacs);
      auto ts = subfacs | transform([](auto&& t) { return t.expr; });
      auto tn = TensorNetworkV2(ts);
      auto canon =
          tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels());
      hash::combine(h, canon.hash_value());
      bool const scalar_result = canon.named_indices_canonical.empty();
      if (scalar_result) {
        return {EvalOp::Product,         //
                ResultType::Scalar,      //
                dummy::make_variable(),  //
                {},                      //
                canon.phase,             //
                h};
      } else {
        auto idxs = get_unique_indices(Product(ts));
        return {EvalOp::Product,     //
                ResultType::Tensor,  //
                dummy::make_tensor(bra(idxs.bra), ket(idxs.ket), aux(idxs.aux)),
                canon.get_indices<Index::index_vector>(),  //
                canon.phase,                               //
                h};
      }
    }
  };

  if (prod.scalar() == 1) {
    return fold_left_to_node(factors | move, make_prod);
  } else {
    auto left = fold_left_to_node(factors | move, make_prod);
    auto right = binarize(Constant{prod.scalar()});

    auto expr = left->is_tensor()     ? dummy::make_tensor(left->as_tensor(),
                                                           /*with_symm = */ false)
                : left->is_constant() ? (left->expr() * right->expr())
                                      : dummy::make_variable();
    auto type = left->is_tensor() ? ResultType::Tensor : ResultType::Scalar;

    auto h = left->hash_value();
    hash::combine(h, right->hash_value());
    auto result = EvalExpr{EvalOp::Product,        //
                           type,                   //
                           expr,                   //
                           left->canon_indices(),  //
                           1,                      //
                           h};
    return EvalExprNode{std::move(result), std::move(left), std::move(right)};
  }
}

namespace impl {

EvalExprNode binarize(ExprPtr const& expr) {
  if (expr->is<Constant>())  //
    return binarize(expr->as<Constant>());

  if (expr->is<Variable>())  //
    return binarize(expr->as<Variable>());

  if (expr->is<Tensor>())  //
    return binarize(expr->as<Tensor>());

  if (expr->is<Sum>())  //
    return binarize(expr->as<Sum>());

  if (expr->is<Product>())  //
    return binarize(expr->as<Product>());

  throw std::logic_error("Encountered unsupported expression in binarize.");
}

}  // namespace impl

}  // namespace sequant
