#include "eval_expr.hpp"

namespace sequant::utils {

size_t eval_expr::hash() const { return hash_; }

eval_expr::eval_op eval_expr::op() const { return op_; }

const Tensor& eval_expr::tensor() const { return tensor_; }

const Constant& eval_expr::scalar() const { return scalar_; }

eval_expr::eval_expr(const Tensor& tnsr) : op_{eval_op::Id}, tensor_{tnsr} {
  if (auto canon_biprod = tensor_.canonicalize(); canon_biprod)
    scalar_ *= canon_biprod->as<Constant>();

  hash_ = eval_expr::hash_terminal_tensor(tensor_);
}

eval_expr::eval_expr(const eval_expr& expr1, const eval_expr& expr2) {
  const auto& sxpr1 = expr1.tensor();
  const auto& sxpr2 = expr2.tensor();

  // canonicalization based on hash value
  // re-order eval sequence ( -- doesn't swap values though)
  // to obtain the equivalence of the kind
  // AB == BA for a binary operation between A and B

  auto swap_on = expr2.hash() < expr1.hash();

  op_ = eval_expr::infer_eval_op(sxpr1, sxpr2);
  tensor_ = swap_on ? eval_expr::make_imed_expr(expr2, expr1, op())
                    : eval_expr::make_imed_expr(expr1, expr2, op());
  // compute phase
  if (auto canon_biprod = tensor_.canonicalize(); canon_biprod)
    scalar_ *= canon_biprod->as<Constant>();

  hash_ = swap_on ? eval_expr::hash_imed(expr2, expr1, op())
                  : eval_expr::hash_imed(expr1, expr2, op());
}

eval_expr::eval_op eval_expr::infer_eval_op(const Tensor& tnsr1,
                                            const Tensor& tnsr2) {
  // HELPER LAMBDA
  // checks if a given index exists in a container
  auto index_exists = [](const auto& container, const auto& index) -> bool {
    return ranges::any_of(container, [l = index.label()](const auto& x) {
      return x.label() == l;
    });
  };
  // HELPER LAMBDA END

  if (tnsr1.rank() == tnsr2.rank()) {
    for (const auto& idx : tnsr1.const_braket())
      if (!(index_exists(tnsr2.const_braket(), idx))) {
        // found an index that doesn't exist
        // in the other tensor
        return eval_op::Prod;
      }

    return eval_op::Sum;  // all indices in any, exist in the other
  }

  // rank not equal: must be a product
  return eval_op::Prod;
}

Symmetry eval_expr::infer_tensor_symmetry_sum(const Tensor& tnsr1,
                                              const Tensor& tnsr2) {
  auto sym1 = tnsr1.symmetry();
  auto sym2 = tnsr2.symmetry();
  if (sym1 == sym2)
    return sym1;  // sum of symm/symm or antisymm/antisymm tensors

  if (sym1 != sym2 && sym1 != Symmetry::nonsymm && sym2 != Symmetry::nonsymm)
    return Symmetry::symm;  // sum of one symmetric and one antisymmetric
  // tensor

  return Symmetry::nonsymm;
}

Symmetry eval_expr::infer_tensor_symmetry_prod(const Tensor& tnsr1,
                                               const Tensor& tnsr2) {
  // HELPER LAMBDA
  // check if all the indices in cont1 are in cont2 AND vice versa
  auto all_common_indices = [](const auto& cont1, const auto& cont2) -> bool {
    if (cont1.size() != cont2.size()) return false;

    return (cont1 | ranges::to<container::set<Index, Index::LabelCompare>>) ==
           (cont2 | ranges::to<container::set<Index, Index::LabelCompare>>);
  };
  // //////

  bool whole_bk_contracted = (all_common_indices(tnsr1.bra(), tnsr2.ket()) ||
                              all_common_indices(tnsr1.ket(), tnsr2.bra()));

  // sym/sym or antisym/antisym with whole braket contraction
  if (whole_bk_contracted && tnsr1.symmetry() == tnsr2.symmetry())
    return tnsr1.symmetry();

  // non symmetric intermediate
  return Symmetry::nonsymm;
}

ParticleSymmetry eval_expr::infer_particle_symmetry(Symmetry s) {
  return (s == Symmetry::symm || s == Symmetry::antisymm)
             ? ParticleSymmetry::symm
             : ParticleSymmetry::nonsymm;
}

BraKetSymmetry eval_expr::infer_braket_symmetry() {
  return get_default_context().braket_symmetry();
}

eval_expr::braket_type eval_expr::target_braket_sum(const Tensor& tnsr1,
                                                    const Tensor& tnsr2) {
  return eval_expr::braket_type{tnsr1.bra(), tnsr1.ket()};
}

eval_expr::braket_type eval_expr::target_braket_prod(const Tensor& tnsr1,
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

size_t eval_expr::hash_braket(
    const decltype(std::declval<Tensor>().const_braket())& braket) {
  size_t bkHash = 0;
  for (auto&& ispace : braket | ranges::views::transform([](const Index& x) {
                         return x.space().type().to_int32();
                       })) {
    hash::combine(bkHash, hash::value(ispace));
  }
  return bkHash;
}

size_t eval_expr::hash_terminal_tensor(const Tensor& tnsr) {
  size_t tHash = 0;
  hash::combine(tHash, hash::value<std::wstring>(tnsr.label().data()));
  hash::combine(tHash, eval_expr::hash_braket(tnsr.const_braket()));
  return tHash;
}

size_t eval_expr::hash_tensor_pair_topology(const Tensor& tnsr1,
                                            const Tensor& tnsr2) {
  size_t hashTopo = 0;
  for (auto&& [pos1, idx1] : tnsr1.const_braket() | ranges::views::enumerate)
    for (auto&& [pos2, idx2] : tnsr2.const_braket() | ranges::views::enumerate)
      if (idx1.label() == idx2.label())
        hash::combine(hashTopo, hash::value(std::pair(pos1, pos2)));
  return hashTopo;
}

size_t eval_expr::hash_imed(const eval_expr& expr1, const eval_expr& expr2,
                            eval_op op) {
  size_t imedHash = 0;
  hash::combine(imedHash, expr1.hash());
  hash::combine(imedHash, expr2.hash());

  const auto& t1 = expr1.tensor();
  const auto& t2 = expr2.tensor();

  hash::combine(imedHash, eval_expr::hash_braket(t1.const_braket()));
  hash::combine(imedHash, eval_expr::hash_braket(t2.const_braket()));
  hash::combine(imedHash, eval_expr::hash_tensor_pair_topology(t1, t2));

  if (op == eval_op::Sum) {
    hash::combine(imedHash, expr1.scalar());
    hash::combine(imedHash, expr2.scalar());
  }

  return imedHash;
}

Tensor eval_expr::make_imed_expr(const eval_expr& expr1, const eval_expr& expr2,
                                 eval_op op) {
  assert(op != eval_op::Id);

  const auto& t1 = expr1.tensor();
  const auto& t2 = expr2.tensor();

  auto [bra, ket] = (op == eval_op::Sum)                        //
                        ? eval_expr::target_braket_sum(t1, t2)  //
                        : eval_expr::target_braket_prod(t1, t2);

  Symmetry tensorSym;  // init only

  if ((expr1.hash() == expr2.hash()) &&
      (bra.size() + ket.size() == 2 * (t1.rank() + t2.rank()))) {
    // outer product
    tensorSym = Symmetry::symm;
  } else {
    tensorSym = (op == eval_op::Sum)
                    ? eval_expr::infer_tensor_symmetry_sum(t1, t2)
                    : eval_expr::infer_tensor_symmetry_prod(t1, t2);
  }

  return Tensor{L"I",
                bra,
                ket,
                tensorSym,
                eval_expr::infer_braket_symmetry(),
                eval_expr::infer_particle_symmetry(tensorSym)};
}

}  // namespace sequant::utils
