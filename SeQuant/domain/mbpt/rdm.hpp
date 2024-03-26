//
// Created by Conner Masteran on 7/1/21.
//

#ifndef SEQUANT_DOMAIN_MBPT_RDM_HPP
#define SEQUANT_DOMAIN_MBPT_RDM_HPP

#include <SeQuant/domain/mbpt/antisymmetrizer.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

namespace sequant {
namespace mbpt {
namespace decompositions {

ExprPtr cumu_to_density(ExprPtr ex_) {
  assert(ex_->is<Tensor>());
  assert(ex_->as<Tensor>().rank() == 1);
  assert(ex_->as<Tensor>().label() == optype2label.at(OpType::RDMCumulant));
  auto down_0 = ex_->as<Tensor>().ket()[0];
  auto up_0 = ex_->as<Tensor>().bra()[0];

  auto density = ex<Tensor>(
      optype2label.at(OpType::RDM), std::initializer_list<Index>{up_0},
      std::initializer_list<Index>{down_0}, std::initializer_list<Index>{});
  return density;
}

ExprPtr cumu2_to_density(ExprPtr ex_) {
  assert(ex_->is<Tensor>());
  assert(ex_->as<Tensor>().rank() == 2);
  assert(ex_->as<Tensor>().label() == optype2label.at(OpType::RDMCumulant));

  auto down_0 = ex_->as<Tensor>().ket()[0];
  auto up_0 = ex_->as<Tensor>().bra()[0];
  auto down_1 = ex_->as<Tensor>().ket()[1];
  auto up_1 = ex_->as<Tensor>().bra()[1];

  const auto rdm_label = optype2label.at(OpType::RDM);
  auto density2 =
      ex<Tensor>(rdm_label, std::initializer_list<Index>{up_0, up_1},
                 std::initializer_list<Index>{down_0, down_1},
                 std::initializer_list<Index>{});
  auto density_1 = ex<Tensor>(rdm_label, std::initializer_list<Index>{up_0},
                              std::initializer_list<Index>{down_0},
                              std::initializer_list<Index>{});
  auto density_2 = ex<Tensor>(rdm_label, std::initializer_list<Index>{up_1},
                              std::initializer_list<Index>{down_1},
                              std::initializer_list<Index>{});

  auto d1_d2 = antisymmetrize(density_1 * density_2);
  return density2 + ex<Constant>(-1) * d1_d2.result;
}

ExprPtr cumu3_to_density(ExprPtr ex_) {
  assert(ex_->is<Tensor>());
  assert(ex_->as<Tensor>().rank() == 3);
  assert(ex_->as<Tensor>().label() == optype2label.at(OpType::RDMCumulant));

  auto down_0 = ex_->as<Tensor>().ket()[0];
  auto up_0 = ex_->as<Tensor>().bra()[0];
  auto down_1 = ex_->as<Tensor>().ket()[1];
  auto up_1 = ex_->as<Tensor>().bra()[1];
  auto down_2 = ex_->as<Tensor>().ket()[2];
  auto up_2 = ex_->as<Tensor>().bra()[2];

  const auto rdm_label = optype2label.at(OpType::RDM);
  auto cumulant2 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                              std::initializer_list<Index>{up_1, up_2},
                              std::initializer_list<Index>{down_1, down_2},
                              std::initializer_list<Index>{});
  auto density_1 = ex<Tensor>(rdm_label, std::initializer_list<Index>{up_0},
                              std::initializer_list<Index>{down_0},
                              std::initializer_list<Index>{});
  auto density_2 = ex<Tensor>(rdm_label, std::initializer_list<Index>{up_1},
                              std::initializer_list<Index>{down_1},
                              std::initializer_list<Index>{});
  auto density_3 = ex<Tensor>(rdm_label, std::initializer_list<Index>{up_2},
                              std::initializer_list<Index>{down_2},
                              std::initializer_list<Index>{});
  auto density3 =
      ex<Tensor>(rdm_label, std::initializer_list<Index>{up_0, up_1, up_2},
                 std::initializer_list<Index>{down_0, down_1, down_2},
                 std::initializer_list<Index>{});

  auto d1_d2 =
      antisymmetrize(density_1 * density_2 * density_3 + density_1 * cumulant2);
  auto temp_result = density3 * ex<Constant>(-1) * d1_d2.result;

  for (auto&& product : temp_result->as<Sum>().summands()) {
    for (auto&& factor : product->as<Product>().factors()) {
      if (factor->is<Tensor>() &&
          (factor->as<Tensor>().label() ==
           optype2label.at(OpType::RDMCumulant)) &&
          (factor->as<Tensor>().rank() == 2)) {
        factor = cumu2_to_density(factor);
      }
    }
  }
  for (auto&& product : temp_result->as<Sum>().summands()) {
    for (auto&& factor : product->as<Product>().factors()) {
      if (factor->is<Tensor>() &&
          factor->as<Tensor>().label() ==
              optype2label.at(OpType::RDMCumulant) &&
          factor->as<Tensor>().rank() == 1) {
        factor = cumu_to_density(factor);
      }
    }
  }
  return temp_result;
}

ExprPtr one_body_sub(
    ExprPtr ex_) {  // J. Chem. Phys. 132, 234107 (2010);
                    // https://doi.org/10.1063/1.3439395 eqn 15 for
  assert(ex_->is<FNOperator>());
  assert(ex_->as<FNOperator>().rank() == 1);
  auto down_0 = ex_->as<FNOperator>().annihilators()[0].index();
  auto up_0 = ex_->as<FNOperator>().creators()[0].index();

  const auto a = ex<FNOperator>(std::initializer_list<Index>{up_0},
                                std::initializer_list<Index>{down_0});
  const auto cumu1 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                                std::initializer_list<Index>{down_0},
                                std::initializer_list<Index>{up_0},
                                std::initializer_list<Index>{});

  auto result = a + (ex<Constant>(-1) * cumu1);
  return (result);
}

ExprPtr two_body_decomp(
    ExprPtr ex_, bool approx = false) {  // J. Chem. Phys. 132, 234107 (2010);
                                         // https://doi.org/10.1063/1.3439395
                                         // eqn 16 for \tilde{a}^{pr}_{qs}
  assert(ex_->is<FNOperator>());
  assert(ex_->as<FNOperator>().rank() == 2);

  auto down_0 = ex_->as<FNOperator>().annihilators()[0].index();
  auto down_1 = ex_->as<FNOperator>().annihilators()[1].index();

  auto up_0 = ex_->as<FNOperator>().creators()[0].index();
  auto up_1 = ex_->as<FNOperator>().creators()[1].index();

  const auto cumu1 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                                std::initializer_list<Index>{down_0},
                                std::initializer_list<Index>{up_0},
                                std::initializer_list<Index>{});
  const auto cumu2 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                                std::initializer_list<Index>{down_1},
                                std::initializer_list<Index>{up_1},
                                std::initializer_list<Index>{});
  const auto a = ex<FNOperator>(std::initializer_list<Index>{up_1},
                                std::initializer_list<Index>{down_1});
  const auto a2 = ex<FNOperator>(std::initializer_list<Index>{up_0, up_1},
                                 std::initializer_list<Index>{down_0, down_1});
  const auto double_cumu = ex<Tensor>(
      optype2label.at(OpType::RDMCumulant),
      std::initializer_list<Index>{down_0, down_1},
      std::initializer_list<Index>{up_0, up_1}, std::initializer_list<Index>{});

  auto term1 = cumu1 * a;
  auto term2 = cumu1 * cumu2;
  auto term3 = double_cumu;

  auto sum_of_terms = antisymmetrize(term1 + term2 + term3);
  sum_of_terms.result = ex<Constant>(-1) * sum_of_terms.result;
  auto result = a2 + sum_of_terms.result;
  return (result);
}

// express 3-body term as sums of 1 and 2-body term. as described in J. Chem.
// Phys. 132, 234107 (2010); https://doi.org/10.1063/1.3439395 eqn 17.
std::pair<ExprPtr, std::pair<std::vector<Index>, std::vector<Index>>>
three_body_decomp(ExprPtr ex_, bool approx = true) {
  assert(ex_->is<FNOperator>());
  assert(ex_->as<FNOperator>().rank() == 3);

  auto down_0 = ex_->as<FNOperator>().annihilators()[0].index();
  auto down_1 = ex_->as<FNOperator>().annihilators()[1].index();
  auto down_2 = ex_->as<FNOperator>().annihilators()[2].index();

  std::vector<Index> initial_lower{down_0, down_1, down_2};

  auto up_0 = ex_->as<FNOperator>().creators()[0].index();
  auto up_1 = ex_->as<FNOperator>().creators()[1].index();
  auto up_2 = ex_->as<FNOperator>().creators()[2].index();

  std::vector<Index> initial_upper{up_0, up_1, up_2};

  const auto cumulant = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                                   std::initializer_list<Index>{down_0},
                                   std::initializer_list<Index>{up_0},
                                   std::initializer_list<Index>{});
  const auto a = ex<FNOperator>(std::initializer_list<Index>{up_1, up_2},
                                std::initializer_list<Index>{down_1, down_2});
  auto a_cumulant = cumulant * a;

  auto cumulant2 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                              std::initializer_list<Index>{down_1},
                              std::initializer_list<Index>{up_1},
                              std::initializer_list<Index>{});
  auto cumulant3 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                              std::initializer_list<Index>{down_2},
                              std::initializer_list<Index>{up_2},
                              std::initializer_list<Index>{});
  auto cumulant_3x = cumulant * cumulant2 * cumulant3;

  auto a1 = ex<FNOperator>(std::initializer_list<Index>{up_0},
                           std::initializer_list<Index>{down_0});
  auto a1_cumu1_cumu2 = a1 * cumulant2 * cumulant3;

  auto two_body_cumu = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                                  std::initializer_list<Index>{down_1, down_2},
                                  std::initializer_list<Index>{up_1, up_2},
                                  std::initializer_list<Index>{});
  auto a1_cumu2 = a1 * two_body_cumu;

  auto cumu1_cumu2 = cumulant * two_body_cumu;
  auto sum_of_terms = antisymmetrize(a_cumulant + cumulant_3x + a1_cumu1_cumu2 +
                                     a1_cumu2 + cumu1_cumu2);

  if (!approx) {
    auto cumu3 =
        ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                   std::initializer_list<Index>{down_0, down_1, down_2},
                   std::initializer_list<Index>{up_0, up_1, up_2},
                   std::initializer_list<Index>{});

    sum_of_terms.result = cumu3 + sum_of_terms.result;
  }

  auto temp_result = sum_of_terms.result;
  simplify(temp_result);
  // std::wcout << "result before substitiutions: " <<
  // to_latex_align(temp_result) << std::endl;

  for (auto&& product :
       temp_result->as<Sum>().summands()) {  // replace all the two body terms
                                             // with one body terms.
    if (product->is<Product>()) {
      for (auto&& factor : product->as<Product>().factors()) {
        if (factor->is<FNOperator>() && factor->as<FNOperator>().rank() == 2) {
          factor = two_body_decomp(factor);
        }
      }
    } else {
    }
  }
  simplify(temp_result);
  for (auto&& product :
       temp_result->as<Sum>().summands()) {  // replace the one body terms with
                                             // the substituted expression
    if (product->is<Product>()) {
      for (auto&& factor : product->as<Product>().factors()) {
        if (factor->is<FNOperator>() && factor->as<FNOperator>().rank() == 1) {
          factor = one_body_sub(factor);
        }
      }
    }
  }
  std::pair<std::vector<Index>, std::vector<Index>> initial_pairing(
      initial_lower, initial_upper);
  std::pair<ExprPtr, std::pair<std::vector<Index>, std::vector<Index>>> result(
      temp_result, initial_pairing);
  // simplify(temp_result);
  // std::wcout << "result before substitiutions: " <<
  // to_latex_align(temp_result,20,7) << std::endl;
  return result;
}

std::pair<ExprPtr, std::pair<std::vector<Index>, std::vector<Index>>>
three_body_decomposition(ExprPtr ex_, int rank, bool fast = false) {
  std::pair<std::vector<Index>, std::vector<Index>> initial_pairing;
  if (rank == 3) {
    auto ex_pair = three_body_decomp(ex_);
    ex_ = ex_pair.first;
    initial_pairing = ex_pair.second;
    simplify(ex_);
    for (auto&& product : ex_->as<Sum>().summands()) {
      if (product->is<Product>()) {
        for (auto&& factor : product->as<Product>().factors()) {
          if (factor->is<Tensor>()) {
            if (factor->as<Tensor>().label() ==
                    optype2label.at(OpType::RDMCumulant) &&
                factor->as<Tensor>().rank() == 3) {
              factor = cumu3_to_density(factor);
            } else if (factor->as<Tensor>().label() ==
                           optype2label.at(OpType::RDMCumulant) &&
                       factor->as<Tensor>().rank() == 2) {
              factor = cumu2_to_density(factor);
            } else if (factor->as<Tensor>().label() ==
                           optype2label.at(OpType::RDMCumulant) &&
                       factor->as<Tensor>().rank() == 1) {
              factor = cumu_to_density(factor);
            } else {
              assert(factor->as<Tensor>().label() !=
                     optype2label.at(OpType::RDMCumulant));
            }
          }
        }
      }
    }
    simplify(ex_);

    }
    else if (rank == 2) {
      if (fast) {
        assert(ex_->is<FNOperator>());
        //FNOp does not store a list of indices so I have to do this
        auto down_0 = ex_->as<FNOperator>().annihilators()[0].index();
        auto down_1 = ex_->as<FNOperator>().annihilators()[1].index();
        auto down_2 = ex_->as<FNOperator>().annihilators()[2].index();

        std::vector<Index> initial_lower{down_0, down_1, down_2};

        auto up_0 = ex_->as<FNOperator>().creators()[0].index();
        auto up_1 = ex_->as<FNOperator>().creators()[1].index();
        auto up_2 = ex_->as<FNOperator>().creators()[2].index();

        std::vector<Index> initial_upper{up_0, up_1, up_2};
        initial_pairing.first = initial_lower; initial_pairing.second = initial_upper;
        // make tensors which can be decomposed into the constituent pieces later in the procedure.
        auto DE2 = ex<Tensor>(L"DE2", IndexList{down_0,down_1,down_2}, IndexList{up_0,up_1,up_2}, IndexList{});
        auto DDE = ex<Tensor>(L"DDE", IndexList{down_0,down_1,down_2}, IndexList{up_0,up_1,up_2}, IndexList{});
        auto D2E =  ex<Tensor>(L"D2E", IndexList{down_0,down_1,down_2}, IndexList{up_0,up_1,up_2}, IndexList{});
        auto result = DE2 + D2E - ex<Constant>(2) * DDE;
        return {result,initial_pairing};
      }
      auto ex_pair = three_body_decomp(ex_, true);
      ex_ = ex_pair.first;
      initial_pairing = ex_pair.second;
      simplify(ex_);
      for (auto&& product : ex_->as<Sum>().summands()) {
        if (product->is<Product>()) {
          for (auto&& factor : product->as<Product>().factors()) {
            if (factor->is<Tensor>()) {
              if (factor->as<Tensor>().label() ==
                      optype2label.at(OpType::RDMCumulant) &&
                  factor->as<Tensor>().rank() > 2) {
                factor = ex<Constant>(0);
              } else if (factor->as<Tensor>().label() ==
                             optype2label.at(OpType::RDMCumulant) &&
                         factor->as<Tensor>().rank() == 2) {
                factor = cumu2_to_density(factor);
              } else if (factor->as<Tensor>().label() ==
                         optype2label.at(OpType::RDMCumulant)) {
                factor = cumu_to_density(factor);
              } else {
                assert(factor->as<Tensor>().label() !=
                       optype2label.at(OpType::RDMCumulant));
              }
            }
          }
        }
      }
      simplify(ex_);
      // std::wcout << " cumulant replacment: " << to_latex_align(_ex,20, 7) <<
      // std::endl;
    }
    else if (rank == 1) {
      auto ex_pair = three_body_decomp(ex_, true);
      ex_ = ex_pair.first;
      initial_pairing = ex_pair.second;
      simplify(ex_);
      for (auto&& product : ex_->as<Sum>().summands()) {
        if (product->is<Product>()) {
          for (auto&& factor : product->as<Product>().factors()) {
            if (factor->is<Tensor>()) {
              if (factor->as<Tensor>().label() ==
                      optype2label.at(OpType::RDMCumulant) &&
                  factor->as<Tensor>().rank() > 1) {
                factor = ex<Constant>(0);
              } else if (factor->as<Tensor>().label() ==
                         optype2label.at(OpType::RDMCumulant)) {
                factor = cumu_to_density(factor);
              } else {
                assert(factor->as<Tensor>().label() !=
                       optype2label.at(OpType::RDMCumulant));
              }
            }
          }
        }
      }
      simplify(ex_);
    }
    else {
      throw "rank not supported!";
    }
    return {ex_, initial_pairing};

}

// in general a three body substitution can be approximated with 1, 2, or 3 body
// terms(3 body has no approximation). this is achieved by replacing densities
// with with particle number > rank by the each successive cumulant
// approximation followed by neglect of the particle rank sized term.
// TODO this implementation is ambitious and currently we only support rank 2
// decompositions.
//
// fast implementation represent non-constant solution interms of like terms and
// permutation operators.
ExprPtr three_body_substitution(ExprPtr& input, int rank, bool fast = false) {
  // just return back if the input is zero.
  if (input == ex<Constant>(0)) {
    return input;
  }
  if(fast){
    assert(rank == 2);
    if (input->is<Sum>()) {
      for (auto&& product : input->as<Sum>().summands()) {
        if (product->is<Product>()) {
          for (auto&& factor : product->as<Product>().factors()) {
            if (factor->is<FNOperator>() && (factor->as<FNOperator>().rank() ==3)) {  // find the 3-body terms
              auto fac_pair = decompositions::three_body_decomposition(factor,rank,fast);  // decompose that term and replace the existing term.
              factor = fac_pair.first;
            }
          }
        }
      }
      simplify(input);
      return input;
    }
    else if (input->is<Product>()) {
      for (auto&& factor : input->as<Product>().factors()) {
        if (factor->is<FNOperator>() &&
            (factor->as<FNOperator>().rank() == 3)) {  // find the 3-body terms
          auto fac_pair = decompositions::three_body_decomposition(
              factor,
              rank,fast);  // decompose that term and replace the existing term.
          factor = fac_pair.first;
        }
      }
      simplify(input);
      return input;
    }
    else if (input->is<FNOperator>()) {
      auto fac_pair = decompositions::three_body_decomposition(
          input, rank,fast);  // decompose that term and replace the existing term.
      input = fac_pair.first;
      simplify(input);
      return input;
      }
  }
  std::pair<std::vector<Index>, std::vector<Index>> initial_pairing;
  if (input->is<Sum>()) {
    for (auto&& product : input->as<Sum>().summands()) {
      if (product->is<Product>()) {
        for (auto&& factor : product->as<Product>().factors()) {
          if (factor->is<FNOperator>() && (factor->as<FNOperator>().rank() ==
                                           3)) {  // find the 3-body terms
            auto fac_pair = decompositions::three_body_decomposition(
                factor,
                rank);  // decompose that term and replace the existing term.
            factor = fac_pair.first;
            initial_pairing = fac_pair.second;
            if (get_default_context().spbasis() == SPBasis::spinfree) {
              factor = antisymm::spin_sum(initial_pairing.second,
                                          initial_pairing.first, factor);
              non_canon_simplify(factor);
            } else {
              throw " wrong spin basis";
            }
          }
        }
      }
    }
  } else if (input->is<Product>()) {
    for (auto&& factor : input->as<Product>().factors()) {
      if (factor->is<FNOperator>() &&
          (factor->as<FNOperator>().rank() == 3)) {  // find the 3-body terms
        auto fac_pair = decompositions::three_body_decomposition(
            factor,
            rank);  // decompose that term and replace the existing term.
        factor = fac_pair.first;
        initial_pairing =
            fac_pair
                .second;  // decompose that term and replace the existing term.
        if (get_default_context().spbasis() == SPBasis::spinfree) {
          factor = antisymm::spin_sum(initial_pairing.second,
                                      initial_pairing.first, factor);
          non_canon_simplify(factor);
        }
      }
    }
  } else if (input->is<FNOperator>()) {
    auto fac_pair = decompositions::three_body_decomposition(
        input, rank);  // decompose that term and replace the existing term.
    input = fac_pair.first;
    initial_pairing = fac_pair.second;
    if (get_default_context().spbasis() == SPBasis::spinfree) {
      // std::wcout << to_latex_align(input,20) << std::endl;
      input = antisymm::spin_sum(initial_pairing.second, initial_pairing.first,
                                 input);
      non_canon_simplify(input);
    }
  } else {
    throw "cannot handle this type";
  }

  return input;
};
}  // namespace decompositions
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_RDM_HPP
