#include "eval_expr.hpp"

#include <cassert>

namespace sequant {

EvalExpr::hash_t EvalExpr::hash() const { return hash_; }

EvalOp EvalExpr::op() const { return op_; }

const Tensor& EvalExpr::tensor() const { return tensor_; }

const Constant& EvalExpr::scalar() const { return scalar_; }

EvalExpr::EvalExpr(const Tensor& tnsr) : op_{EvalOp::Id}, tensor_{tnsr},
hash_{EvalExpr::hash_terminal_tensor(tnsr)}, annot_{braket_to_annot(tnsr.const_braket())}{}

EvalExpr::EvalExpr(const EvalExpr& xpr1,
                   const EvalExpr& xpr2,
                   EvalOp op) {
  assert(op != EvalOp::Id);
  auto const& expr1 = xpr1.hash() < xpr2.hash() ? xpr1 : xpr2;
  auto const& expr2 = xpr1.hash() < xpr2.hash() ? xpr2 : xpr1;

  auto bk = braket_type{};
  if (op == EvalOp::Prod)
    bk = target_braket_prod(expr1.tensor(), expr2.tensor());
  else {
    std::get<0>(bk) = expr1.tensor().bra()
                      | ranges::to<index_container_type>;
    std::get<1>(bk) = expr1.tensor().ket()
                      | ranges::to<index_container_type>;
  }

  auto const& t1 = expr1.tensor();
  auto const& t2 = expr2.tensor();
  Symmetry s = Symmetry::invalid;
  switch (op) {
    case EvalOp::Sum: s = infer_tensor_symmetry_sum(expr1, expr2); break;
    case EvalOp::Prod: s = infer_tensor_symmetry_prod(expr1, expr2); break;
    case EvalOp::Symm: s = Symmetry::symm; break;
    case EvalOp::Antisymm: s = Symmetry::antisymm; break;
    default: assert(false && "Unsupported operation for symmetry detect.");
  }
  op_ = op;

  tensor_ = Tensor{L"I",
      std::get<0>(bk),
      std::get<1>(bk),
      s, infer_braket_symmetry(), infer_particle_symmetry(s)};

  hash_ = hash_imed(expr1, expr2, op);

  annot_ = braket_to_annot(tensor_.const_braket());
}

Symmetry EvalExpr::infer_tensor_symmetry_sum(EvalExpr const& xpr1,
                                             EvalExpr const& xpr2) {
  auto sym1 = xpr1.tensor().symmetry();
  auto sym2 = xpr2.tensor().symmetry();
  if (sym1 == sym2)
    return sym1;  // sum of symm/symm or antisymm/antisymm tensors

  if (sym1 != sym2 && sym1 != Symmetry::nonsymm && sym2 != Symmetry::nonsymm)
    return Symmetry::symm;  // sum of one symmetric and one antisymmetric
  // tensor

  return Symmetry::nonsymm;
}

Symmetry EvalExpr::infer_tensor_symmetry_prod(EvalExpr const& xpr1,
                                              EvalExpr const& xpr2) {
  using index_set_t = container::set<Index, Index::LabelCompare>;
  // HELPER LAMBDA
  // check if all the indices in cont1 are in cont2 AND vice versa
  auto all_common_indices = [](const auto& cont1, const auto& cont2) -> bool {
    return (cont1.size() == cont2.size()) &&
           (cont1 | ranges::to<index_set_t>) ==
               (cont2 | ranges::to<index_set_t>);
  };
  // //////

  auto const& tnsr1 = xpr1.tensor();
  auto const& tnsr2 = xpr2.tensor();

  if (xpr1.hash() == xpr2.hash()){
    // potential outer product

    auto const uniq_idxs = ranges::views::concat(tnsr1.const_braket(),
                                                tnsr2.const_braket())
                          | ranges::to<index_set_t>;

    if (ranges::distance(uniq_idxs) == tnsr1.const_braket().size()
                                     + tnsr2.const_braket().size()) {
      return Symmetry::antisymm;
    }
  }

  bool whole_bk_contracted = (all_common_indices(tnsr1.bra(), tnsr2.ket()) ||
                              all_common_indices(tnsr1.ket(), tnsr2.bra()));

  // sym/sym or antisym/antisym with whole braket contraction
  if (whole_bk_contracted && tnsr1.symmetry() == tnsr2.symmetry())
    return tnsr1.symmetry();

  // non symmetric intermediate
  return Symmetry::nonsymm;
}

ParticleSymmetry EvalExpr::infer_particle_symmetry(Symmetry s) {
  return (s == Symmetry::symm || s == Symmetry::antisymm)
             ? ParticleSymmetry::symm
             : ParticleSymmetry::nonsymm;
}

BraKetSymmetry EvalExpr::infer_braket_symmetry() {
  return get_default_context().braket_symmetry();
}

EvalExpr::braket_type EvalExpr::target_braket_prod(const Tensor& tnsr1,
                                                   const Tensor& tnsr2) {
  using ranges::views::keys;
  using ranges::views::values;
  using ranges::views::zip;

  auto remove_item = [](auto& vec, size_t pos) -> void {
    std::swap(vec[pos], vec.back());
    vec.pop_back();
  };

  auto left = zip(tnsr1.bra(), tnsr1.ket()) | ranges::to_vector;
  auto right = zip(tnsr2.bra(), tnsr2.ket()) | ranges::to_vector;

next_contract:
  while (!right.empty()) {
    for (auto rr = 0; rr < right.size(); ++rr) {
      auto& [rb, rk] = right[rr];
      for (auto ll = 0; ll < left.size(); ++ll) {
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
  }
  // the result is now in left

  return {keys(left) | ranges::to<index_container_type>,
          values(left) | ranges::to<index_container_type>};
}

EvalExpr::hash_t EvalExpr::hash_braket(
    const decltype(std::declval<Tensor>().const_braket())& braket) {
  EvalExpr::hash_t bkHash = 0;
  for (auto const& idx: braket) {
    hash::combine(bkHash, hash::value(idx.space().type().to_int32()));
    hash::combine(bkHash, hash::value(idx.space().qns().to_int32()));
  }
  return bkHash;
}

EvalExpr::hash_t EvalExpr::hash_terminal_tensor(const Tensor& tnsr) {
  EvalExpr::hash_t tHash = 0;
  hash::combine(tHash, hash::value<std::wstring>(tnsr.label().data()));
  hash::combine(tHash, EvalExpr::hash_braket(tnsr.const_braket()));
  return tHash;
}

EvalExpr::hash_t EvalExpr::hash_tensor_pair_topology(const Tensor& tnsr1,
                                                     const Tensor& tnsr2) {
  EvalExpr::hash_t hashTopo = 0;
  for (auto&& [pos1, idx1] : tnsr1.const_braket() | ranges::views::enumerate)
    for (auto&& [pos2, idx2] : tnsr2.const_braket() | ranges::views::enumerate)
      if (idx1.label() == idx2.label())
        hash::combine(hashTopo, hash::value(std::pair(pos1, pos2)));
  return hashTopo;
}

EvalExpr::hash_t EvalExpr::hash_imed(const EvalExpr& expr1,
                                     const EvalExpr& expr2, EvalOp op) {
  EvalExpr::hash_t imedHash = 0;
  hash::combine(imedHash, expr1.hash());
  hash::combine(imedHash, expr2.hash());

  const auto& t1 = expr1.tensor();
  const auto& t2 = expr2.tensor();

  hash::combine(imedHash, EvalExpr::hash_braket(t1.const_braket()));
  hash::combine(imedHash, EvalExpr::hash_braket(t2.const_braket()));
  hash::combine(imedHash, EvalExpr::hash_tensor_pair_topology(t1, t2));

  if (op == EvalOp::Sum) {
    hash::combine(imedHash, expr1.scalar());
    hash::combine(imedHash, expr2.scalar());
  }

  return imedHash;
}

}  // namespace sequant
