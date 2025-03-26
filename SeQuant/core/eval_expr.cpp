#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
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
#include <atomic>
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

std::size_t detail::next_eval_expr_id() {
  static std::atomic<std::size_t> next_id = 0;

  std::size_t id = next_id.load();

  // This ensures that we are updating next_id with the next id while
  // also ensuring that no other thread is currently producing the same
  // id that this one is doing.
  while (!next_id.compare_exchange_weak(id, id + 1)) {
  }

  return id;
}

namespace {

size_t hash_terminal_tensor(Tensor const&) noexcept;

size_t hash_imed(EvalExpr const&, EvalExpr const&, EvalOp) noexcept;

ExprPtr make_imed(EvalExpr const&, EvalExpr const&, EvalOp) noexcept;

bool is_tot(Tensor const& t) noexcept {
  return ranges::any_of(t.const_indices(), &Index::has_proto_indices);
}

std::wstring_view const var_label = L"Z";

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
    : op_type_{EvalOp::Id},
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
    : op_type_{EvalOp::Id},
      result_type_{ResultType::Scalar},
      hash_value_{hash::value(c)},
      expr_{c.clone()} {}

EvalExpr::EvalExpr(Variable const& v)
    : op_type_{EvalOp::Id},
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

EvalOp EvalExpr::op_type() const noexcept { return op_type_; }

ResultType EvalExpr::result_type() const noexcept { return result_type_; }

size_t EvalExpr::hash_value() const noexcept { return hash_value_; }

std::size_t EvalExpr::id() const noexcept { return id_; }

ExprPtr EvalExpr::expr() const noexcept { return expr_; }

void EvalExpr::set_expr(ExprPtr expr) { expr_ = std::move(expr); }

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

Tensor const& EvalExpr::as_tensor() const noexcept {
  return expr().as<Tensor>();
}

Constant const& EvalExpr::as_constant() const noexcept {
  return expr().as<Constant>();
}

Variable const& EvalExpr::as_variable() const noexcept {
  return expr().as<Variable>();
}

std::string EvalExpr::label() const noexcept {
  if (is_tensor())
    return to_string(as_tensor().label()) + "(" + indices_annot() + ")";
  else if (is_constant()) {
    return sequant::to_string(sequant::to_latex(as_constant()));
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

///
/// \return hash value to identify the connectivity between a pair of tensors.
///
/// @note Let [(i,j)] be the list of ordered pair of index positions that are
///       connected. i is the position in the indices of the first tensor (T1)
///       and j is that of the second tensor (T2). Then this function combines
///       the hash values of the elements of this list.
///
/// @warning O(N^2) algorithm
///
size_t hash_tensor_pair_topology(Tensor const& t1, Tensor const& t2) noexcept {
  using ranges::views::enumerate;
  size_t h = 0;
  for (auto&& [pos1, idx1] : t1.const_indices() | enumerate)
    for (auto&& [pos2, idx2] : t2.const_indices() | enumerate)
      if (idx1.label() == idx2.label())
        hash::combine(h, hash::value(std::pair(pos1, pos2)));
  return h;
}

size_t hash_terminal_tensor(Tensor const& tnsr) noexcept {
  size_t h = 0;
  hash::combine(h, hash::value(tnsr.label()));
  hash::combine(h, hash_indices(tnsr.const_indices()));
  return h;
}

size_t hash_imed(EvalExpr const& left, EvalExpr const& right,
                 EvalOp op) noexcept {
  size_t h = 0;
  hash::combine(h, hash::value(op));

  auto lh = hash::value(left);
  auto rh = hash::value(right);

  hash::combine(h, lh < rh ? lh : rh);
  hash::combine(h, lh < rh ? rh : lh);

  if (left.result_type() == ResultType::Tensor &&
      right.result_type() == ResultType::Tensor)
    hash::combine(h, hash_tensor_pair_topology(left.expr()->as<Tensor>(),
                                               right.expr()->as<Tensor>()));
  return h;
}

Symmetry tensor_symmetry_sum(EvalExpr const& left,
                             EvalExpr const& right) noexcept {
  auto const& t1 = left.expr()->as<Tensor>();
  auto const& t2 = right.expr()->as<Tensor>();

  auto sym1 = t1.symmetry();
  auto sym2 = t2.symmetry();
  if (sym1 == sym2)
    return sym1;  // sum of symm/symm or antisymm/antisymm tensors

  // sum of one symmetric and one antisymmetric tensor
  if (sym1 != sym2 && sym1 != Symmetry::nonsymm && sym2 != Symmetry::nonsymm)
    return Symmetry::symm;

  return Symmetry::nonsymm;
}

Symmetry tensor_symmetry_prod(EvalExpr const& left,
                              EvalExpr const& right) noexcept {
  using index_set_t = container::set<Index, Index::LabelCompare>;

  // HELPER LAMBDA
  // check if all the indices in cont1 are in cont2 AND vice versa
  auto all_common_indices = [](const auto& cont1, const auto& cont2) -> bool {
    return (cont1.size() == cont2.size()) &&
           (cont1 | ranges::to<index_set_t>) ==
               (cont2 | ranges::to<index_set_t>);
  };
  // //////

  auto const& tnsr1 = left.expr()->as<Tensor>();
  auto const& tnsr2 = right.expr()->as<Tensor>();

  if (hash::value(left) == hash::value(right)) {
    // potential outer product of the same tensor
    auto const uniq_idxs =
        ranges::views::concat(tnsr1.const_indices(), tnsr2.const_indices()) |
        ranges::to<index_set_t>;

    if (static_cast<std::size_t>(ranges::distance(uniq_idxs)) ==
        tnsr1.const_indices().size() + tnsr2.const_indices().size()) {
      // outer product confirmed
      return Symmetry::antisymm;
    }
  }

  // not an outer product of same tensor confirmed
  auto imed_sym = Symmetry::invalid;
  bool whole_bk_contracted = (all_common_indices(tnsr1.bra(), tnsr2.ket()) ||
                              all_common_indices(tnsr1.ket(), tnsr2.bra()));
  auto sym1 = tnsr1.symmetry();
  auto sym2 = tnsr2.symmetry();

  assert(sym1 != Symmetry::invalid);
  assert(sym2 != Symmetry::invalid);

  if (whole_bk_contracted &&
      !(sym1 == Symmetry::nonsymm || sym2 == Symmetry::nonsymm)) {
    imed_sym = sym1 == sym2 ? sym1 : Symmetry::symm;

  } else {
    imed_sym = Symmetry::nonsymm;
  }

  assert(imed_sym != Symmetry::invalid);
  return imed_sym;
}

ParticleSymmetry particle_symmetry(Symmetry s) noexcept {
  return (s == Symmetry::symm || s == Symmetry::antisymm)
             ? ParticleSymmetry::symm
             : ParticleSymmetry::nonsymm;
}

ExprPtr make_sum(EvalExpr const& left, EvalExpr const& right) noexcept {
  assert(left.is_tensor() && right.is_tensor());

  auto const& t1 = left.as_tensor();
  auto const& t2 = right.as_tensor();

  assert(t1.bra_rank() + t1.ket_rank()         //
             == t2.bra_rank() + t2.ket_rank()  //
         && "differing ranks for summed tensors");

  auto ts = tensor_symmetry_sum(left, right);
  auto ps = particle_symmetry(ts);
  auto bks = get_default_context().braket_symmetry();
  return ex<Tensor>(L"I", t1.bra(), t1.ket(), t1.aux(), ts, bks, ps);
}

ExprPtr make_prod(EvalExpr const& left, EvalExpr const& right) noexcept {
  assert(left.is_tensor() && right.is_tensor());

  auto const& t1 = left.as_tensor();
  auto const& t2 = right.as_tensor();

  auto [b, k, a] = get_uncontracted_indices(t1, t2);
  if (b.empty() && k.empty() && a.empty()) {
    // dot product
    return ex<Variable>(var_label);
  } else {
    // regular tensor product
    auto ts = tensor_symmetry_prod(left, right);
    auto ps = particle_symmetry(ts);
    auto bks = get_default_context().braket_symmetry();
    return ex<Tensor>(L"I", bra(std::move(b)), ket(std::move(k)),
                      aux(std::move(a)), ts, bks, ps);
  }
}

ExprPtr make_imed(EvalExpr const& left, EvalExpr const& right,
                  EvalOp op) noexcept {
  assert(op != EvalOp::Id);

  auto lres = left.result_type();
  auto rres = right.result_type();

  if (lres == ResultType::Scalar && rres == ResultType::Scalar) {
    // scalar (+|*) scalar

    return ex<Variable>(var_label);

  } else if (lres == ResultType::Scalar && rres == ResultType::Tensor) {
    // scalar (*) tensor

    assert(op == EvalOp::Prod && "scalar + tensor not supported");
    auto const& t = right.expr()->as<Tensor>();
    return ex<Tensor>(Tensor{L"I", t.bra(), t.ket(), t.aux(), t.symmetry(),
                             t.braket_symmetry(), t.particle_symmetry()});

  } else if (lres == ResultType::Tensor && rres == ResultType::Scalar) {
    // tensor (*) scalar

    return make_imed(right, left, op);

  } else {
    // tensor (+|*) tensor

    auto lh = hash::value(left);
    auto rh = hash::value(right);
    auto const& left_ = lh <= rh ? left : right;
    auto const& right_ = lh <= rh ? right : left;

    if (op == EvalOp::Sum) {
      // tensor (+) tensor
      return make_sum(left_, right_);
    } else {
      // tensor (*) tensor
      return make_prod(left_, right_);
    }
  }
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
      node->is_tensor() && op == EvalOp::Id || op == EvalOp::Sum)
    collect.emplace_back(ExprWithHash{node->expr(), node->hash_value()});
  else if (node->op_type() == EvalOp::Prod && !node.leaf()) {
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

  auto make_prod = [i = 0, hs = imed_hashes(hvals)](
                       EvalExprNode const& left,
                       EvalExprNode const& right) mutable -> EvalExpr {
    auto h = ranges::at(hs, ++i);
    if (left->is_scalar() && right->is_scalar()) {
      // scalar * scalar
      return {
          EvalOp::Prod, ResultType::Scalar, dummy::make_variable(), {}, 1, h};
    } else if (left->is_scalar() || right->is_scalar()) {
      // scalar * tensor or tensor * scalar
      auto const& tl = left->is_tensor() ? left : right;
      auto const& t = tl->as_tensor();
      return {EvalOp::Prod,                                                  //
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
        return {EvalOp::Prod,            //
                ResultType::Scalar,      //
                dummy::make_variable(),  //
                {},                      //
                canon.phase,             //
                h};
      } else {
        auto idxs = get_unique_indices(Product(ts));
        return {EvalOp::Prod,        //
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

    auto h = left->hash_value();
    hash::combine(h, right->hash_value());
    auto result = EvalExpr{EvalOp::Prod,           //
                           left->result_type(),    //
                           left->expr(),           //
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
