#ifndef SEQUANT_INFERRED_IMED_DATA_HPP
#define SEQUANT_INFERRED_IMED_DATA_HPP

#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_network.hpp>

namespace sequant::factorize {

struct InferredImedData {
  container::svector<Index> bra_indices;

  container::svector<Index> ket_indices;

  Symmetry symmetry{Symmetry::nonsymm};

  BraKetSymmetry braket_symmetry = get_default_context().braket_symmetry();

  ParticleSymmetry particle_symmetry{ParticleSymmetry::nonsymm};

  ExprPtr phase = std::make_shared<Constant>(1);

  bool is_sum;

  Expr::hash_type hash_value{};

  InferredImedData(const ExprPtr& left, const ExprPtr& right) {
    auto& ltensor = left->as<Tensor>();
    auto& rtensor = right->as<Tensor>();

    // figure out if it is a sum of left and right tensors
    is_sum = ltensor.rank() == rtensor.rank();
    // not only rank but also all indices in one tensor must be
    // present in the other
    auto index_exists = [](const auto& container, const auto& index) {
      for (const auto& idx : container)
        if (idx == index) return true;
      return false;
    };
    for (const auto& idx : ltensor.const_braket())
      if (!(index_exists(rtensor.const_braket(), idx))) {
        is_sum = false;
        break;
      }  // now is_sum is set properly

    // Set bra ket indices for intermediate
    if (is_sum) {
      std::copy(ltensor.bra().begin(), ltensor.bra().end(),
                std::back_inserter(bra_indices));
      std::copy(ltensor.ket().begin(), ltensor.ket().end(),
                std::back_inserter(ket_indices));
      left_ = left->clone();
      right_ = right->clone();
    } else {
      // is_sum false implies product type operation
      // need to canonicalize
      // figure out target indices
      TensorNetwork::named_indices_t target_indices{};
      for (const auto& idx : ltensor.const_braket()) target_indices.insert(idx);
      for (const auto& idx : rtensor.const_braket()) {
        if (target_indices.contains(idx))
          target_indices.erase(idx);
        else
          target_indices.insert(idx);
      }  // target_indices figured out
      auto tnet = TensorNetwork(*(left->clone() * right->clone()));
      // canonicalize
      phase = tnet.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                                false, &target_indices);
      //
      left_ = std::dynamic_pointer_cast<Tensor>(tnet.tensors().at(0));
      right_ = std::dynamic_pointer_cast<Tensor>(tnet.tensors().at(1));
      //
      // set indices
      auto set_indices = [&target_indices, this](const auto& tnetTnsr) {
        auto tnsr = std::dynamic_pointer_cast<Tensor>(tnetTnsr);
        for (const auto& idx : tnsr->bra())
          if (target_indices.contains(idx)) bra_indices.push_back(idx);
        for (const auto& idx : tnsr->ket())
          if (target_indices.contains(idx)) ket_indices.push_back(idx);
      };
      set_indices(tnet.tensors().at(0));
      set_indices(tnet.tensors().at(1));
      // intermediate indices set
    }

    // Infer Tensor Symmetry
    if (is_sum) {
      if (ltensor.symmetry() == rtensor.symmetry()) {
        if (ltensor.symmetry() == Symmetry::antisymm)
          symmetry = Symmetry::antisymm;  // sum of two antisymm tensors
        else if (ltensor.symmetry() == Symmetry::symm)
          symmetry = Symmetry::symm;  // sum of two symmetric tensors
      } else if (auto sym_tuple =
                     std::make_tuple(ltensor.symmetry(), rtensor.symmetry());
                 sym_tuple ==
                     std::make_tuple(Symmetry::antisymm, Symmetry::symm) ||
                 sym_tuple ==
                     std::make_tuple(Symmetry::symm, Symmetry::antisymm))
        symmetry = Symmetry::symm;  // sum of one symmetric and one
                                    // antisymmetric tensor
    } else {
      // is_sum false implies product type operation
      // since we don't consider partial symmetry, let's check if all of
      // bra(ket) indices from one tensor are contracted with ket(bra) from
      // another tensor
      auto all_common_indices = [&index_exists](const auto& cont1,
                                                const auto& cont2) {
        TensorNetwork::named_indices_t commons;
        for (const auto& idx : cont1) commons.insert(idx);
        for (const auto& idx : cont2) {
          if (commons.contains(idx))
            commons.erase(idx);
          else
            commons.insert(idx);
        }
        return commons.size() == 0;
      };
      bool whole_bk_contracted =
          (all_common_indices(ltensor.bra(), rtensor.ket()) ||
           all_common_indices(ltensor.ket(), rtensor.bra()) ||
           all_common_indices(ltensor.bra(), rtensor.bra()) ||
           all_common_indices(ltensor.ket(), rtensor.ket()));
      if (whole_bk_contracted) {
        if (ltensor.symmetry() == Symmetry::antisymm ||
            rtensor.symmetry() == Symmetry::antisymm)
          symmetry = Symmetry::antisymm;
        else if (ltensor.symmetry() == rtensor.symmetry() &&
                 ltensor.symmetry() == Symmetry::symm)
          symmetry = Symmetry::symm;
      }
      if (ltensor.hash_value() == rtensor.hash_value() &&
          bra_indices.size() + ket_indices.size() ==
              ltensor.rank() + rtensor.rank()) {
        symmetry = Symmetry::symm;
        // outer product of same tensor is symmetric
      }
      // else symmetry is set nonsymm by default
    }
    if (symmetry == Symmetry::antisymm || symmetry == Symmetry::symm)
      particle_symmetry = ParticleSymmetry::symm;  // set particle symmetry
  };

  const auto& left() const { return left_; }
  const auto& right() const { return right_; }

 private:
  ExprPtr left_;
  ExprPtr right_;
};

}  // namespace sequant::factorize

#endif /* SEQUANT_INFERRED_IMED_DATA_HPP */
