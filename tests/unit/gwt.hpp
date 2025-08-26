//
// Created by Eduard Valeyev on 8/7/25.
//

#ifndef SEQUANT_TESTS_UNIT_GWT_HPP
#define SEQUANT_TESTS_UNIT_GWT_HPP

#include <SeQuant/core/container.hpp>

#include <range/v3/algorithm/equal_range.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

#include <vector>

namespace sequant {

/// generalized Wick theorem for a product of particle-conserving operators in
/// general space normal ordered wrt to single-product vacuum
class GWT {
 public:
  /// @param nop_ranks ranks of normal-ordered operators (`nop`s)
  GWT(std::vector<std::size_t> nop_nbody_ranks,
      bool full_contractions_only = true, bool use_topology = true)
      : full_contractions_only_(full_contractions_only),
        use_topology_(use_topology),
        nop_ranks_(nop_nbody_ranks |
                   ranges::views::transform(
                       [](const auto nbody_rank) { return 2 * nbody_rank; }) |
                   ranges::to_vector),
        rank_(ranges::accumulate(nop_ranks_, 0ul)) {
    if (!full_contractions_only_) use_topology_ = false;

    // compute map from op ordinal to nop ordinal
    for (std::size_t nop_ord = 0; nop_ord != nop_ranks_.size(); ++nop_ord) {
      for (std::size_t op = 0; op != nop_ranks_[nop_ord]; ++op) {
        op_ord_to_nop_ord_.push_back(nop_ord);
      }
    }

    // compute all single contractions
    for (std::size_t op_i = 0; op_i != rank_; ++op_i) {
      for (std::size_t op_j = op_i + 1; op_j < rank_;
           op_j += 2) {  // cre-ann or ann-cre contractions only
        if (op_ord_to_nop_ord_[op_i] == op_ord_to_nop_ord_[op_j]) continue;
        pcontrs_.emplace_back(op_i, op_j);
      }
    }

    // initialize result_
    result_.emplace(ranges::views::iota(0ul, rank_) | ranges::to_vector,
                    contrs_t{});

    // recursively apply WT
    {
      result_t current_wf_terms = result_;  // initial wavefront = input
      bool done = false;
      while (!done) {
        result_t next_wf_terms;
        for (const auto& [oper, contrs] : current_wf_terms) {
          // return true if op_it is first op of its type (cre or ann) in its
          // nop
          auto first_in_nop = [this, &oper](const auto& op_it) {
            if (op_it == oper.begin()) return true;
            const auto op_nop_ord = this->op_ord_to_nop_ord_[*op_it];
            const auto op_ann = *op_it % 2 == 0;
            bool same_nop = true;
            auto prev_op_it = op_it - 1;
            while (same_nop) {
              auto prev_op_ann = *prev_op_it % 2 == 0;
              auto prev_op_nop_ord = this->op_ord_to_nop_ord_[*(prev_op_it)];
              same_nop = op_nop_ord == prev_op_nop_ord;
              if (same_nop && op_ann == prev_op_ann) return false;
              if (prev_op_it == oper.begin()) return true;
              --prev_op_it;
            }
            return true;
          };

          for (const auto& pcontr : pcontrs_) {
            const auto& [op1, op2] = pcontr;

            auto [op1_it, op1_it2] = ranges::equal_range(oper, op1);
            if (op1_it == op1_it2) continue;
            if (use_topology_ && !first_in_nop(op1_it)) continue;
            const auto op1_ord = op1_it - oper.begin();

            auto [op2_it, op2_it2] = ranges::equal_range(oper, op2);
            if (op2_it != op2_it2) {
              // if (use_topology_ && !first_in_nop(op2_it)) continue;
              const auto op2_ord = op2_it - oper.begin();
              auto new_oper = oper;

              // n.b. erase in reverse order to keep ordinals stable
              new_oper.erase(new_oper.begin() + op2_ord);
              new_oper.erase(new_oper.begin() + op1_ord);

              contrs_t new_contrs(contrs);
              new_contrs.emplace(pcontr);
              next_wf_terms.emplace(new_oper, new_contrs);
            }
          }
        }
        done = next_wf_terms.empty();
        // if want full contraction and have new terms generated, result_ only
        // contains partial contractions, purge them
        if (full_contractions_only_) {
          if (!done) result_ = next_wf_terms;
        } else {
          for (const auto& new_term : next_wf_terms) {
            result_.emplace(new_term);
          }
        }
        current_wf_terms = std::move(next_wf_terms);
      }
    }

    // if full contractions only, filter out partial contractions from result_
    if (full_contractions_only_) {
      result_ = result_ | ranges::views::filter([](const auto& ops_contrs) {
                  return ops_contrs.oper.empty();
                }) |
                ranges::to<result_t>;
    }
  }

  const auto& result() const { return result_; }

 private:
  const bool full_contractions_only_;
  bool use_topology_;
  const std::vector<std::size_t>
      nop_ranks_;           // op ranks, i.e. 2 body operator has op rank of 4
  const std::size_t rank_;  // total particle rank
  std::vector<std::size_t>
      op_ord_to_nop_ord_;  // maps op ordinal to nop ordinal; odd op ordinals
                           // are cre, even ann

  using pcontr_t =
      std::pair<std::size_t, std::size_t>;  // primitive contraction
  std::vector<pcontr_t>
      pcontrs_;  // list of all possible primitive contractions

  using contrs_t = container::set<pcontr_t>;  // set of primitive contractions
  using oper_t =
      std::vector<std::size_t>;  // general operator is a set of ops (=list of
                                 // op ordinals in increasing order)
  using opers_t = std::vector<oper_t>;  // set of operators
  // contracted operator
  struct ContractedOperator {
    oper_t oper;
    contrs_t contrs;
    ContractedOperator(oper_t o, contrs_t c)
        : oper(std::move(o)), contrs(std::move(c)) {}

    friend bool operator<(const ContractedOperator& a,
                          const ContractedOperator& b) {
      if (a.oper != b.oper) {
        return a.oper < b.oper;
      } else {
        return a.contrs < b.contrs;
      }
    }
  };
  using coper_t = ContractedOperator;
  using result_t =
      container::set<coper_t>;  // result = set of contracted operators

  result_t result_;
};

}  // namespace sequant

#endif  // SEQUANT_TESTS_UNIT_GWT_HPP
