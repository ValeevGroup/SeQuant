//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_WICK_HPP
#define SEQUANT2_WICK_HPP

#include <mutex>
#include <utility>

#include "op.hpp"
#include "tensor.hpp"
#include "ranges.hpp"

namespace sequant2 {

/// Applies Wick's theorem to a sequence of normal-ordered operators.
///
/// @tparam S particle statistics
template<Statistics S>
class WickTheorem {
 public:
  static constexpr const Statistics statistics = S;
  static_assert(S == Statistics::FermiDirac, "WickTheorem not yet implemented for Bose-Einstein");

  explicit WickTheorem(const NormalOperatorSequence<S> &input) : input_(input) {
    assert(input.empty() || input.vacuum() != Vacuum::Invalid);
    assert(input.empty() || input.vacuum() != Vacuum::Invalid);
  }

  /// Controls whether next call to compute() will full contractions only or all (including partial) contractions.
  /// By default compute() generates all contractions.
  /// @param sf if true, will complete full contractions only.
  /// @return reference to @c *this , for daisy-chaining
  WickTheorem &full_contractions(bool fc) {
    full_contractions_ = fc;
    return *this;
  }
  /// Controls whether next call to compute() will assume spin-free or spin-orbital normal-ordered operators
  /// By default compute() assumes spin-orbital operators.
  /// @param sf if true, will complete full contractions only.
  WickTheorem &spinfree(bool sf) {
    spinfree_ = sf;
    return *this;
  }
  /// Controls whether next call to compute() will reduce the result
  /// By default compute() will not perform reduction.
  /// @param r if true, compute() will reduce the result.
  WickTheorem &reduce(bool r) {
    reduce_ = r;
    return *this;
  }

  /// Specifies the external indices; by default assume all indices are summed over
  /// @param ext_inds external (nonsummed) indices
  WickTheorem &set_external_indices(IndexList external_indices) {
    external_indices_ = external_indices;
    return *this;
  }

  /// Computes and returns the result
  /// @param count_only if true, will return a vector of default-initialized values, useful if only interested in the total count
  /// @return the result of applying Wick's theorem, i.e. a sum of {prefactor, normal operator} pairs
  ExprPtr
  compute(const bool count_only = false) const {
    if (!full_contractions_)
      throw std::logic_error("WickTheorem::compute: full_contractions=false not yet supported");
    if (spinfree_)
      throw std::logic_error("WickTheorem::compute: spinfree=true not yet supported");
    auto result = compute_nontensor_wick(count_only);
    if (reduce_ && !count_only) {
      reduce(result);
      canonicalize(result);
    }
    return std::move(result);
  }

 private:
  const NormalOperatorSequence<S> &input_;
  bool full_contractions_ = false;
  bool spinfree_ = false;
  bool reduce_ = false;
  container::vector<Index> external_indices_;

  /// carries state down the stack of recursive calls
  struct NontensorWickState {
    NontensorWickState(const NormalOperatorSequence<S>& opseq) : opseq(opseq), level(0), count_only(false) {
      compute_size();
    }
    NormalOperatorSequence<S> opseq;  //!< current state of operator sequence
    std::size_t opseq_size;  //!< current size of opseq
    Product sp;  //!< current prefactor
    int level;  //!< level in recursive wick call stack
    bool count_only;  //!< if true, only update result size

    void compute_size() {
      opseq_size = 0;
      for(const auto& op: opseq)
        opseq_size += op.size();
    }
    void reset(const NormalOperatorSequence<S>& o) {
      sp = Product{};
      opseq = o;
      compute_size();
    }
  };
  /// Applies most naive version of Wick's theorem, where sign rule involves counting Ops
  ExprPtr
  compute_nontensor_wick(const bool count_only) const {
    std::vector<std::pair<Product, NormalOperator<S>>> result;  //!< current value of the result
    std::mutex mtx;  // used in critical sections updating the result
    auto result_plus_mutex = std::make_pair(&result, &mtx);
    NontensorWickState state(input_);
    state.count_only = count_only;

    recursive_nontensor_wick(result_plus_mutex, state);

    // convert result to an Expr
    // if result.size() == 0, return null ptr
    // TODO revise if decide to use Constant(0)
    ExprPtr result_expr;
    if (result.size() == 1) {  // if result.size() == 1, return Product
      result_expr = make<Product>(std::move(result[0].first));
    }
    else if (result.size() > 1) {
      auto sum = std::make_shared<Sum>();
      for(auto& term: result) {
        sum->append(make<Product>(std::move(term.first)));
      }
      result_expr = sum;
    }
    return result_expr;
  }

  void recursive_nontensor_wick(std::pair<std::vector<std::pair<Product, NormalOperator<S>>>*, std::mutex*>& result,
                                NontensorWickState& state) const {
    // if full contractions needed, make contractions involving first index with another index, else contract any index i with index j (i<j)
    if (full_contractions_) {
      using opseq_view_type = flattened_rangenest<NormalOperatorSequence<S>>;
      auto opseq_view = opseq_view_type(&state.opseq);
      using std::begin;
      using std::end;
      auto opseq_view_begin = begin(opseq_view);

      // optimization: can't contract fully if first op is not a qp annihilator
      if (!is_qpannihilator(*opseq_view_begin, input_.vacuum()))
        return;

      auto op_iter = opseq_view_begin;
      ++op_iter;
      for(; op_iter != end(opseq_view);) {
        if (op_iter != opseq_view_begin && ranges::get_cursor(op_iter).range_iter() != ranges::get_cursor(opseq_view_begin).range_iter()) {
          if (can_contract(*opseq_view_begin, *op_iter, input_.vacuum())) {
//            std::wcout << "level " << state.level << ": contracting " << to_latex(*opseq_view_begin) << " with " << to_latex(*op_iter) << std::endl;
//            std::wcout << "  current opseq = " << to_latex(state.opseq) << std::endl;

            // update the phase, if needed
            double phase = 1;
            if (statistics == Statistics::FermiDirac) {
              const auto
                  distance = ranges::get_cursor(op_iter).ordinal() - ranges::get_cursor(opseq_view_begin).ordinal() - 1;
              if (distance % 2) {
                phase *= -1;
              }
            }

            // update the prefactor and opseq
            Product sp_copy = state.sp;
            state.sp.append(phase, contract(*opseq_view_begin, *op_iter, input_.vacuum()));
            // remove from back to front
            Op<S> right = *op_iter;
            ranges::get_cursor(op_iter).erase();
            --state.opseq_size;
            Op<S> left = *opseq_view_begin;
            ranges::get_cursor(opseq_view_begin).erase();
            --state.opseq_size;

//            std::wcout << "  opseq after contraction = " << to_latex(state.opseq) << std::endl;

            // update the result if nothing left to contract and have a nonzero result
            if (state.opseq_size == 0 && !state.sp.empty()) {
              result.second->lock();
//              std::wcout << "got " << to_latex(state.sp) << std::endl;
              if (!state.count_only)
                result.first->push_back(std::make_pair(std::move(state.sp), NormalOperator<S>{}));
              else
                result.first->resize(result.first->size() + 1);
//              std::wcout << "now up to " << result.first->size() << " terms" << std::endl;
              result.second->unlock();
            }

            if (state.opseq_size != 0) {
              ++state.level;
              recursive_nontensor_wick(result, state);
              --state.level;
            }

            // restore the prefactor and opseq
            state.sp = std::move(sp_copy);
            // restore from front to back
            ranges::get_cursor(opseq_view_begin).insert(std::move(left));
            ++state.opseq_size;
            ranges::get_cursor(op_iter).insert(std::move(right));
            ++state.opseq_size;
//            std::wcout << "  restored opseq = " << to_latex(state.opseq) << std::endl;

          }
          ++op_iter;
        }
        else
          ++op_iter;
      }
    }
    else
      assert(false);   // full_contraction_=false not implemented yet, should result in error earlier
  }

 public:
  static bool can_contract(const Op<S>& left, const Op<S>& right, Vacuum vacuum = get_default_context().vacuum()) {
    if (is_qpannihilator<S>(left, vacuum) && is_qpcreator<S>(right, vacuum)) {
      const auto qpspace_left = qpannihilator_space<S>(left, vacuum);
      const auto qpspace_right = qpcreator_space<S>(right, vacuum);
      const auto qpspace_common = intersection(qpspace_left, qpspace_right);
      if (qpspace_common != IndexSpace::null_instance())
        return true;
    }
    return false;
  }

  static std::shared_ptr<Expr> contract(const Op<S>& left, const Op<S>& right, Vacuum vacuum = get_default_context().vacuum()) {
    assert(can_contract(left, right, vacuum));
    assert(!left.index().has_proto_indices() && !right.index().has_proto_indices());  // I don't think the logic is correct for dependent indices
    if (is_pure_qpannihilator<S>(left, vacuum) && is_pure_qpcreator<S>(right, vacuum))
      return overlap(left.index(), right.index());
    else {
      const auto qpspace_left = qpannihilator_space<S>(left, vacuum);
      const auto qpspace_right = qpcreator_space<S>(right, vacuum);
      const auto qpspace_common = intersection(qpspace_left, qpspace_right);
      const auto index_common = Index::make_tmp_index(qpspace_common);
      if (qpspace_common != left.index().space() && qpspace_common != right.index().space()) {  // may need 2 overlaps if neither space is pure qp creator/annihilator
        auto result = std::make_shared<Product>();
        result->append(1, overlap(left.index(), index_common));
        result->append(1, overlap(index_common, right.index()));
        return result;
      }
      else {
        return overlap(left.index(), right.index());
      }
    }
  }

 public:  // TODO make these members private once WickTheorem can work on full expressions (not on sequences of normal operators) directly
  /// @param[in,out] on input, Wick theorem result, on output the result of reducing the overlaps
  void reduce(ExprPtr& expr) const;

};

using BWickTheorem = WickTheorem<Statistics::BoseEinstein>;
using FWickTheorem = WickTheorem<Statistics::FermiDirac>;

}  // namespace sequant2

#include "wick.impl.hpp"

#endif //SEQUANT2_WICK_HPP
