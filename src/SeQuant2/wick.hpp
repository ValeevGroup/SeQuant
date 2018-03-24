//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_WICK_HPP
#define SEQUANT2_WICK_HPP

#include <utility>

#include "op.hpp"
#include "tensor.hpp"

namespace sequant2 {

/// Applies Wick's theorem to a sequence of normal-ordered operators.
///
/// @tparam S particle statistics
template<Statistics S>
class WickTheorem {
 public:
  static_assert(S == Statistics::FermiDirac, "WickTheorem not yet implemented for Bose-Einstein");

  explicit WickTheorem(const NormalOperatorSequence<S> &input) : input_(input) {
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

  std::vector<std::pair<ScaledProduct, NormalOperator<S>>>
  compute() const {
    if (full_contractions_)
      throw std::logic_error("WickTheorem::compute: full_contractions=true not yet supported");
    if (spinfree_)
      throw std::logic_error("WickTheorem::compute: spinfree=true not yet supported");
    return compute_nontensor_wick();
  }

 private:
  const NormalOperatorSequence<S> &input_;
  bool full_contractions_ = false;
  bool spinfree_ = false;

  /// Applies most naive version of Wick's theorem, where sign rule involves counting Ops
  std::vector<std::pair<ScaledProduct, NormalOperator<S>>>
  compute_nontensor_wick() const {
    return std::vector<std::pair<ScaledProduct, NormalOperator<S>>>{};
  };
};

using BWickTheorem = WickTheorem<Statistics::BoseEinstein>;
using FWickTheorem = WickTheorem<Statistics::FermiDirac>;

}

#endif //SEQUANT2_WICK_HPP
