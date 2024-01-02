#include "eval_expr.hpp"
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/wstring.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/functional.hpp>
#include <range/v3/view.hpp>

#include <cmath>

namespace sequant {

namespace {

size_t hash_terminal_tensor(Tensor const&) noexcept;

size_t hash_imed(EvalExpr const&, EvalExpr const&, EvalOp) noexcept;

ExprPtr make_imed(EvalExpr const&, EvalExpr const&, EvalOp) noexcept;

bool is_tot(Tensor const& t) noexcept {
  return ranges::any_of(t.const_braket(), &Index::has_proto_indices);
}

std::wstring_view const var_label = L"Z";

}  // namespace

NestedTensorIndices::NestedTensorIndices(const sequant::Tensor& tnsr) {
  auto push_ix = [this](Index const& ix) {
    if (ix.has_proto_indices())
      inner.push_back(ix);
    else
      outer.push_back(ix);
  };

  for (auto const& ix : tnsr.const_braket()) {
    push_ix(ix);
    for (auto const& ix_proto : ix.proto_indices()) push_ix(ix_proto);
  }

  if (!inner.empty()) {
    ranges::actions::stable_sort(outer, Index::LabelCompare{});
    ranges::actions::unique(outer, [](Index const& ix1, Index const& ix2) {
      return ix1.label() == ix2.label();
    });
  }
}

std::string EvalExpr::braket_annot() const noexcept {
  if (!is_tensor()) return {};

  // given an iterable of sequant::Index objects, returns a string made
  // of their full labels separted by comma
  //   eg. (a_1^{i_1,i_2},a_2^{i_2,i_3}) -> "a_1i_1i_2,a_2i_2i_3"
  //   eg. (i_1, i_2) -> "i_1,i_2"
  auto annot = [](auto&& ixs) -> std::string {
    using namespace ranges::views;

    auto full_labels = ixs                              //
                       | transform(&Index::full_label)  //
                       | transform([](auto&& fl) {      //
                           return sequant::to_string(fl);
                         });
    return full_labels                      //
           | intersperse(std::string{","})  //
           | join                           //
           | ranges::to<std::string>;
  };

  auto nested = NestedTensorIndices{as_tensor()};

  return nested.inner.empty()  //
             ? annot(nested.outer)
             : annot(nested.outer) + ";" + annot(nested.inner);
}

size_t EvalExpr::global_id_{};

EvalExpr::EvalExpr(Tensor const& tnsr)
    : op_type_{EvalOp::Id},
      result_type_{ResultType::Tensor},
      hash_value_{hash_terminal_tensor(tnsr)},
      id_{},
      expr_{tnsr.clone()},
      tot_{is_tot(tnsr)} {}

EvalExpr::EvalExpr(Constant const& c)
    : op_type_{EvalOp::Id},
      result_type_{ResultType::Scalar},
      hash_value_{hash::value(c)},
      id_{},
      expr_{c.clone()},
      tot_{false} {}

EvalExpr::EvalExpr(Variable const& v)
    : op_type_{EvalOp::Id},
      result_type_{ResultType::Scalar},
      hash_value_{hash::value(v)},
      id_{},
      expr_{v.clone()},
      tot_{false} {}

EvalExpr::EvalExpr(EvalExpr const& left, EvalExpr const& right, EvalOp op)
    : op_type_{op},
      hash_value_{hash_imed(left, right, op)},
      id_{++global_id_},
      expr_{make_imed(left, right, op)} {
  result_type_ = expr_->is<Tensor>() ? ResultType::Tensor : ResultType::Scalar;
  tot_ = expr_->is<Tensor>() && is_tot(expr_->as<Tensor>());
}

EvalOp EvalExpr::op_type() const noexcept { return op_type_; }

ResultType EvalExpr::result_type() const noexcept { return result_type_; }

size_t EvalExpr::hash_value() const noexcept { return hash_value_; }

size_t EvalExpr::id() const noexcept { return id_; }

ExprPtr EvalExpr::expr() const noexcept { return expr_; }

bool EvalExpr::tot() const noexcept { return tot_; }

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
    return to_string(as_tensor().label()) + "(" + braket_annot() + ")";
  else if (is_constant()) {
    auto const& c = as_constant();
    auto real = Constant{c.value().real()}.value<double>();
    auto imag = Constant{c.value().imag()}.value<double>();
    assert(real != 0 || imag != 0);
    std::string r = std::to_string(real);
    std::string i = std::to_string(imag);
    if (real == 0) return i;
    if (imag == 0) return r;
    return "(" + r + "," + i + ")";
  } else {
    assert(is_variable());
    return to_string(as_variable().label());
  }
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
size_t hash_braket(T const& bk) noexcept {
  size_t h = 0;
  for (auto const& idx : bk) {
    hash::combine(h, hash::value(idx.space().type().to_int32()));
    hash::combine(h, hash::value(idx.space().qns().to_int32()));
  }
  return h;
}

///
/// \return hash value to identify the connectivity between a pair of tensors.
///
/// @note Let [(i,j)] be the list of ordered pair of index positions that are
///       connected. i is the position in the braket of the first tensor (T1)
///       and j is that of the second tensor (T2). Then this function combines
///       the hash values of the elements of this list.
///
/// @warning O(N^2) algorithm
///
size_t hash_tensor_pair_topology(Tensor const& t1, Tensor const& t2) noexcept {
  using ranges::views::enumerate;
  size_t h = 0;
  for (auto&& [pos1, idx1] : t1.const_braket() | enumerate)
    for (auto&& [pos2, idx2] : t2.const_braket() | enumerate)
      if (idx1.label() == idx2.label())
        hash::combine(h, hash::value(std::pair(pos1, pos2)));
  return h;
}

size_t hash_terminal_tensor(Tensor const& tnsr) noexcept {
  size_t h = 0;
  hash::combine(h, hash::value(tnsr.label()));
  hash::combine(h, hash_braket(tnsr.const_braket()));
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

std::pair<container::svector<Index>,  // bra
          container::svector<Index>   // ket
          >
target_braket(Tensor const& t1, Tensor const& t2) noexcept {
  using ranges::views::keys;
  using ranges::views::values;
  using ranges::views::zip;
  using index_container = container::svector<Index>;

  auto remove_item = [](auto& vec, size_t pos) -> void {
    std::swap(vec[pos], vec.back());
    vec.pop_back();
  };

  auto left = zip(t1.bra(), t1.ket()) | ranges::to_vector;
  auto right = zip(t2.bra(), t2.ket()) | ranges::to_vector;

  while (!right.empty()) {
    for (std::size_t rr = 0; rr < right.size(); ++rr) {
      auto& [rb, rk] = right[rr];
      for (std::size_t ll = 0; ll < left.size(); ++ll) {
        auto& [lb, lk] = left[ll];
        if (lb == rk && rb == lk) {
          remove_item(left, ll);
          remove_item(right, rr);
          goto next_contract;
        } else if (lb == rk) {
          std::swap(lk, rk);
          remove_item(left, ll);
          goto next_contract;
        } else if (lk == rb) {
          std::swap(lb, rb);
          remove_item(left, ll);
          goto next_contract;
        }
      }  // ll
      left.emplace_back(rb, rk);
      remove_item(right, rr);
    }  // rr
  next_contract:
      /* just point to the start of the while loop */;
  }
  // the result is now in left

  return {keys(left) | ranges::to<index_container>,
          values(left) | ranges::to<index_container>};
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
        ranges::views::concat(tnsr1.const_braket(), tnsr2.const_braket()) |
        ranges::to<index_set_t>;

    if (static_cast<std::size_t>(ranges::distance(uniq_idxs)) ==
        tnsr1.const_braket().size() + tnsr2.const_braket().size()) {
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
  return ex<Tensor>(L"I", t1.bra(), t1.ket(), ts, bks, ps);
}

ExprPtr make_prod(EvalExpr const& left, EvalExpr const& right) noexcept {
  assert(left.is_tensor() && right.is_tensor());

  auto const& t1 = left.as_tensor();
  auto const& t2 = right.as_tensor();

  auto [bra, ket] = target_braket(t1, t2);
  if (bra.empty() && ket.empty()) {
    // dot product
    return ex<Variable>(var_label);
  } else {
    // regular tensor product
    auto ts = tensor_symmetry_prod(left, right);
    auto ps = particle_symmetry(ts);
    auto bks = get_default_context().braket_symmetry();
    return ex<Tensor>(L"I", bra, ket, ts, bks, ps);
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
    return ex<Tensor>(Tensor{L"I", t.bra(), t.ket(), t.symmetry(),
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

}  // namespace sequant
