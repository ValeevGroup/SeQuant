#include "eval_expr.hpp"
#include "SeQuant/core/sequant.hpp"

namespace sequant::utils {

size_t eval_expr::hash() const { return hash_; }

eval_expr::eval_op eval_expr::op() const { return op_; }

const Tensor& eval_expr::tensor() const { return tensor_; }

const Constant& eval_expr::phase() const { return phase_; }

const Constant& eval_expr::scalar() const { return scalar_; }

void eval_expr::scale(const Constant& scal) { scalar_ = scal; }

eval_expr::eval_expr(Tensor tnsr) : op_{eval_op::Id}, tensor_{std::move(tnsr)} {
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  if (auto canon_biprod = tensor_.canonicalize(); canon_biprod)
    phase_ *= canon_biprod->as<Constant>();

  hash_ = eval_expr::hash_terminal_tensor(tensor_);
}

eval_expr::eval_expr(const eval_expr& expr1, const eval_expr& expr2) {
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  const auto& sxpr1 = expr1.tensor();
  const auto& sxpr2 = expr2.tensor();

  // canonicalization based on hash value
  // re-order eval sequence ( -- doesn't swap values though)
  // to obtain the equivalence of the kind
  // AB == BA for a binary operation between A and B

  auto swap_on = expr2.hash() < expr1.hash();

  op_ = eval_expr::infer_eval_op(sxpr1, sxpr2);
  tensor_ = std::move(swap_on ? eval_expr::make_imed_expr(expr2, expr1, op())
                              : eval_expr::make_imed_expr(expr1, expr2, op()));
  // compute phase
  phase_ *= expr1.phase();
  phase_ *= expr2.phase();
  if (auto canon_biprod = tensor_.canonicalize(); canon_biprod)
    phase_ *= canon_biprod->as<Constant>();

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
  //
  // TODO:
  // write a better algorithm
  //
  bool swap_on = tnsr1.rank() > tnsr2.rank();

  auto target_bra = (swap_on ? tnsr2.bra() : tnsr1.bra()) |
                    ranges::to<eval_expr::index_container_type>;

  auto target_ket = (swap_on ? tnsr2.ket() : tnsr1.ket()) |
                    ranges::to<eval_expr::index_container_type>;

  const auto& incoming_bra = swap_on ? tnsr1.bra() : tnsr2.bra();
  const auto& incoming_ket = swap_on ? tnsr1.ket() : tnsr2.ket();

  for (auto&& [ipos, iparticle] : ranges::views::enumerate(
           ranges::views::zip(incoming_bra, incoming_ket))) {
    for (auto&& [tpos, tparticle] :
         ranges::views::enumerate(ranges::views::zip(target_bra, target_ket))) {
      if (iparticle.first.label() == tparticle.second.label()) {
        // bra index of incoming particle contracted
        // with ket index of target particle
        if (iparticle.second.label() == tparticle.first.label()) {
          // AND, ket index contracted with bra index
          target_bra.erase(target_bra.begin() + tpos);
          target_ket.erase(target_ket.begin() + tpos);

        } else {
          // only bra index of incoming particle contracted
          // with the ket index of the target particle
          target_ket[tpos] = incoming_ket[ipos];
        }
        goto next_incoming_particle;
      } else if (iparticle.second.label() == tparticle.first.label()) {
        // ket index of incoming particle contracted
        // with bra index of target particle
        target_bra[tpos] = incoming_bra[ipos];
        goto next_incoming_particle;
      } else {
        continue;  // to next target particle
      }
    }
    // got an incoming paticle which is not contracted at all
    target_bra.emplace_back(incoming_bra.at(ipos));
    target_ket.emplace_back(incoming_ket.at(ipos));

  next_incoming_particle:;
  }

  return {target_bra, target_ket};
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

  auto tensorSym = Symmetry::invalid;  // init only

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
