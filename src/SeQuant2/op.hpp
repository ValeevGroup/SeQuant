//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT2_OP_H
#define SEQUANT2_OP_H

#include <numeric>

#include <range/v3/all.hpp>

#include "index.hpp"
#include "sequant.hpp"
#include "vacuum.hpp"
#include "iterator.hpp"

namespace sequant2 {

enum class Statistics { BoseEinstein, FermiDirac };

enum class Action { create, annihilate };

/// applies (Hermitian) adjoint to @c action
inline Action adjoint(Action action) { return action == Action::create ? Action::annihilate : Action::create; }

/// @brief Op is a creator/annihilator operator
///
/// Op = Index + Action
/// @tparam S specifies the particle statistics
template<Statistics S = Statistics::FermiDirac>
class Op {
 public:
  static constexpr Statistics statistics = S;

  Op() = default;
  Op(Index index, Action action) noexcept : index_(index), action_(action) {}
  Op(std::wstring_view index_label, Action action) noexcept : index_(Index{index_label}), action_(action) {}

  Index index() const { return index_; }
  Action action() const { return action_; }
  /// applies (Hermitian) adjoint to this operator
  Op &adjoint() {
    action_ = sequant2::adjoint(action_);
    return *this;
  }

 private:
  Index index_;
  Action action_;
};

template<Statistics S>
bool operator==(const Op<S> &op1, const Op<S> &op2) {
  return op1.index() == op2.index() && op1.action() == op2.action();
}

template<Statistics S>
bool operator!=(const Op<S> &op1, const Op<S> &op2) {
  return !(op1 == op2);
}

/// @brief Operator is a sequence of Op objects
///
/// @tparam S specifies the particle statistics
template<Statistics S = Statistics::FermiDirac>
class Operator : public std::vector<Op<S>> {
 public:
  static constexpr Statistics statistics = S;

  Operator() = default;
  explicit Operator(std::initializer_list<Op<S>> ops)
      : std::vector<Op<S>>(ops) {}
  explicit Operator(std::vector<Op<S>> &&ops)
      : std::vector<Op<S>>(std::move(ops)) {}
  Operator(Action action, std::initializer_list<Index> indices)
      : std::vector<Op<S>>(make_ops(action, indices)) {}
  Operator(Action action, std::initializer_list<std::wstring_view> index_labels)
      : std::vector<Op<S>>(make_ops(action, index_labels)) {}

  operator const std::vector<Op<S>> &() const & { return *this; }
  operator std::vector<Op<S>> &&() && { return *this; }

  /// applies (Hermitian) adjoint operation to this
  /// @return reference to @c *this , for daisy-chaining
  Operator &adjoint() {
    std::reverse(this->begin(), this->end());
    std::for_each(this->begin(), this->end(), [](Op<S> &op) { op.adjoint(); });
    return *this;
  }

 private:
  std::vector<Op<S>> make_ops(Action action,
                              std::initializer_list<Index> indices) {
    std::vector<Op<S>> result;
    result.reserve(indices.size());
    for (const auto &idx : indices)
      result.emplace_back(idx, action);
    return result;
  }

  std::vector<Op<S>>
  make_ops(Action action,
           std::initializer_list<std::wstring_view> index_labels) {
    std::vector<Op<S>> result;
    result.reserve(index_labels.size());
    for (const auto &idx_label : index_labels)
      result.emplace_back(Index{idx_label}, action);
    return result;
  }
};

/// @brief NormalOperator is an Operator normal-ordered with respect to vacuum.

/// @note Normal ordering means all creators are to the left of all annihilators. It is natural to express
/// at least number-conserving normal operators (i.e. those with equal number of creators and annihilators)
/// as tensors with creators as superscripts and annihilators as subscripts. Operator
/// cre(p1) cre(p2) ... cre(pN) ann(qN) ... ann(q2) ann(q1) is represented in such notation as a^{p1 p2 ... pN}_{q1 q2 ... qN},
/// hence it is natural to specify annihilators in the order of their particle index (i.e. as q1 q2, etc.) which is
/// reverse of the order of their appearance in Operator.
///
/// @note The tensor notation becomes less intuitive for number non-conserving operators, e.g. cre(p1) cre(p2) ann(q2) could
/// be represented as a^{p1 p2}_{q2 ⎵} or a^{p1 p2}_{⎵ q2}. To make it explicit that ann(q2) acts on same particle as
/// cre(p2) the latter notation is used; similarly, cre(p1) ann(q1) ann(q2) is represented as a^{⎵ p1}_{q1 q2}.
///
/// @tparam S specifies the particle statistics
template<Statistics S = Statistics::FermiDirac>
class NormalOperator : public Operator<S> {
 public:
  static constexpr Statistics statistics = S;

  /// constructs an identity operator
  NormalOperator(Vacuum v = get_default_context().vacuum()) {}

  /// @param creators sequence of creators
  /// @param annihilators sequence of annihilators (in order of particle indices, see the class documentation for more info).
  NormalOperator(std::initializer_list<Op<S>> creators,
                 std::initializer_list<Op<S>> annihilators,
                 Vacuum v = get_default_context().vacuum())
      : Operator<S>{}, vacuum_(v), ncreators_(size(creators)) {
    for(const auto& op: creators) {
      assert(op.action() == Action::create);
    }
    for(const auto& op: annihilators) {
      assert(op.action() == Action::annihilate);
    }
    this->reserve(size(creators) + size(annihilators));
    this->insert(this->end(), cbegin(creators), cend(creators) );
    this->insert(this->end(), crbegin(annihilators), crend(annihilators) );
  }

  /// @param creator_indices sequence of creator indices
  /// @param annihilator_indices sequence of annihilator indices (in order of particle indices, see the class documentation for more info).
  NormalOperator(std::initializer_list<Index> creator_indices,
                 std::initializer_list<Index> annihilator_indices,
                 Vacuum v = get_default_context().vacuum())
      : Operator<S>{}, vacuum_(v), ncreators_(size(creator_indices)) {
    this->reserve(size(creator_indices) + size(annihilator_indices));
    for(const auto& i: creator_indices) {
      this->emplace_back(i, Action::create);
    }
    for(const auto& i: annihilator_indices|ranges::view::reverse) {
      this->emplace_back(i, Action::annihilate);
    }
  }

  /// @param creator_index_labels sequence of creator index labels
  /// @param annihilator_index_labels sequence of annihilator index labels (in order of particle indices, see the class documentation for more info).
  NormalOperator(std::initializer_list<std::wstring_view> creator_index_labels,
                 std::initializer_list<std::wstring_view> annihilator_index_labels,
                 Vacuum v = get_default_context().vacuum())
      : Operator<S>{}, vacuum_(v), ncreators_(size(creator_index_labels)) {
    this->reserve(size(creator_index_labels) + size(annihilator_index_labels));
    for(const auto& l: creator_index_labels) {
      this->emplace_back(Index{l}, Action::create);
    }
    for(const auto& l: annihilator_index_labels|ranges::view::reverse) {
      this->emplace_back(Index{l}, Action::annihilate);
    }
  }

  /// @return the vacuum state with respect to which the operator is normal-ordered.
  Vacuum vacuum() const { return vacuum_; }
  /// @return the range of creators
  auto creators() const { return ranges::view::counted(this->cbegin(), ncreators()); }
  /// @return the range of annihilators
  auto annihilators() const { return ranges::view::counted(this->crbegin(), nannihilators()); }
  /// @return the number of creators
  auto ncreators() const { return ncreators_; }
  /// @return the number of annihilators
  auto nannihilators() const { return this->size() - ncreators(); }

  NormalOperator &adjoint() {
    static_cast<Operator<S>&>(*this).adjoint();
    ncreators_ = this->size() - ncreators_;
    return *this;
  }

 private:
  Vacuum vacuum_;
  std::size_t ncreators_ = 0;
};

template<Statistics S>
bool operator==(const NormalOperator<S> &op1, const NormalOperator<S> &op2) {
  return op1.vacuum() == op2.vacuum() && ranges::equal(op1,op2);
}

/// @brief NormalOperatorSequence is a sequence NormalOperator objects, all ordered with respect to same vacuum
///
/// @tparam S specifies the particle statistics
template<Statistics S = Statistics::FermiDirac>
class NormalOperatorSequence : public std::vector<NormalOperator<S>> {
 public:
  static constexpr Statistics statistics = S;

  NormalOperatorSequence(std::initializer_list<NormalOperator<S>> operators)
      : std::vector<NormalOperator<S>>(operators) {
    check_vacuum();
  }

  Vacuum vacuum() const { return vacuum_; }

  operator const std::vector<NormalOperator<S>> &() const & { return *this; }
  operator std::vector<NormalOperator<S>> &&() && { return *this; }

  /// applies (Hermitian) adjoint operation to this
  /// @return reference to @c *this , for daisy-chaining
  NormalOperatorSequence &adjoint() {
    std::reverse(this->begin(), this->end());
    std::for_each(this->begin(), this->end(), [](NormalOperator<S> &op) { op.adjoint(); });
    return *this;
  }

 private:
  Vacuum vacuum_ = Vacuum::Invalid;
  /// ensures that all operators use same vacuum, and sets vacuum_
  void check_vacuum() {
    vacuum_ =
        std::accumulate(this->cbegin(), this->cend(), Vacuum::Invalid, [](Vacuum v1, const NormalOperator<S> &v2) {
          if (v1 == Vacuum::Invalid) {
            return v2.vacuum();
          } else {
            if (v1 != v2.vacuum())
              throw std::invalid_argument(
                  "NormalOperatorSequence expects all constituent NormalOperator objects to use same vacuum");
            else
              return v1;
          }
        });
  }
};

using BOp = Op<Statistics::BoseEinstein>;
using BOperator = Operator<Statistics::BoseEinstein>;
using BNOperator = NormalOperator<Statistics::BoseEinstein>;
using BNOperatorSeq = NormalOperatorSequence<Statistics::BoseEinstein>;
using FOp = Op<Statistics::FermiDirac>;
using FOperator = Operator<Statistics::FermiDirac>;
using FNOperator = NormalOperator<Statistics::FermiDirac>;
using FNOperatorSeq = NormalOperatorSequence<Statistics::FermiDirac>;

inline BOp bcre(Index i) { return BOp(i, Action::create); }
inline BOp bcre(std::wstring_view i) { return BOp(Index{i}, Action::create); }
inline BOp bcre(std::wstring_view i,
                std::initializer_list<std::wstring_view> pi) {
  return BOp(Index(i, pi), Action::create);
}
inline BOp bann(Index i) { return BOp(i, Action::annihilate); }
inline BOp bann(std::wstring_view i) {
  return BOp(Index{i}, Action::annihilate);
}
inline BOp bann(std::wstring_view i,
                std::initializer_list<std::wstring_view> pi) {
  return BOp(Index(i, pi), Action::annihilate);
}
inline FOp fcre(Index i) { return FOp(i, Action::create); }
inline FOp fcre(std::wstring_view i) { return FOp(Index{i}, Action::create); }
inline FOp fcre(std::wstring_view i,
                std::initializer_list<std::wstring_view> pi) {
  return FOp(Index(i, pi), Action::create);
}
inline FOp fann(Index i) { return FOp(i, Action::annihilate); }
inline FOp fann(std::wstring_view i) {
  return FOp(Index{i}, Action::annihilate);
}
inline FOp fann(std::wstring_view i,
                std::initializer_list<std::wstring_view> pi) {
  return FOp(Index(i, pi), Action::annihilate);
}

template<Statistics S>
std::wstring to_latex(const Op<S> &op) {
  std::wstring result;
  result = L"{";
  result += (S == Statistics::FermiDirac ? L"a" : L"b");
  result += (op.action() == Action::create ? L"^{\\dagger}_" : L"_");
  result += to_latex(op.index());
  result += L"}";
  return result;
}

template<Statistics S>
std::wstring to_latex(const Operator<S> &op) {
  std::wstring result;
  result = L"{";
  for (const auto &o : op)
    result += to_latex(o);
  result += L"}";
  return result;
}

template<Statistics S>
std::wstring to_latex(const NormalOperator<S> &op) {
  std::wstring result;
  result = L"{";
  result += (S == Statistics::FermiDirac
             ? (op.vacuum() == Vacuum::Physical ? L"a" : L"\\tilde{a}")
             : (op.vacuum() == Vacuum::Physical ? L"b" : L"\\tilde{b}"));
  result += L"^{";
  const auto ncreators = op.ncreators();
  const auto nannihilators = op.nannihilators();
  if (ncreators < nannihilators) // pad on the left with square underbrackets, i.e. ⎵
    for(auto i = 0; i!=(nannihilators-ncreators); ++i)
      result += L"\\textvisiblespace\\,";
  for (const auto &o : op.creators())
    result += to_latex(o.index());
  result += L"}_{";
  if (ncreators > nannihilators) // pad on the left with square underbrackets, i.e. ⎵
    for(auto i = 0; i!=(ncreators-nannihilators); ++i)
      result += L"\\textvisiblespace\\,";
  for (const auto &o : op.annihilators())
    result += to_latex(o.index());
  result += L"}}";
  return result;
}

template<Statistics S>
std::wstring to_latex(const NormalOperatorSequence<S> &opseq) {
  std::wstring result;
  result = L"{";
  for (const auto &op: opseq)
    result += to_latex(op);
  result += L"}";
  return result;
}

} // namespace sequant2

#endif // SEQUANT2_OP_H
