//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT2_OP_H
#define SEQUANT2_OP_H

#include <numeric>

#include "index.hpp"
#include "sequant.hpp"
#include "vacuum.hpp"

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

/// @brief Operator is a list of Op objects
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

/// @brief NormalOperator is an Operator normal-ordered with respect to vacuum
///
/// @tparam S specifies the particle statistics
template<Statistics S = Statistics::FermiDirac>
class NormalOperator {
 public:
  static constexpr Statistics statistics = S;

  NormalOperator(Vacuum v = get_default_context().vacuum()) {}

  NormalOperator(std::initializer_list<Operator<S>> creators,
                 std::initializer_list<Operator<S>> annihilators,
                 Vacuum v = get_default_context().vacuum())
      : vacuum_(v), creators_(creators), annihilators_(annihilators) {}

  NormalOperator(std::initializer_list<Index> creator_indices,
                 std::initializer_list<Index> annihilator_indices,
                 Vacuum v = get_default_context().vacuum())
      : vacuum_(v), creators_(Operator<S>(Action::create, creator_indices)),
        annihilators_(Operator<S>(Action::annihilate, annihilator_indices)) {}

  NormalOperator(std::initializer_list<std::wstring_view> creator_index_labels,
                 std::initializer_list<std::wstring_view> annihilator_index_labels,
                 Vacuum v = get_default_context().vacuum())
      : vacuum_(v), creators_(Operator<S>(Action::create, creator_index_labels)),
        annihilators_(Operator<S>(Action::annihilate, annihilator_index_labels)) {}

  operator Operator<S>() const & {
    std::vector<Op<S>> grandlist(creators_.size() + annihilators_.size());
    std::copy(cbegin(creators_), cend(creators_), begin(grandlist));
    std::copy(crbegin(annihilators_), crend(annihilators_), begin(grandlist) + creators_.size());
    return Operator<S>(std::move(grandlist));
  }

  operator Operator<S>() && {
    std::vector<Op<S>> grandlist(std::move(creators_));
    const auto ncre = grandlist.size();
    grandlist.resize(ncre + annihilators_.size());
    std::copy(crbegin(annihilators_), crend(annihilators_), begin(grandlist) + ncre);
    return Operator<S>(std::move(grandlist));
  }

  Vacuum vacuum() const { return vacuum_; }
  const Operator<S> &creators() const { return creators_; }
  const Operator<S> &annihilators() const { return annihilators_; }

  /// applies (Hermitian) adjoint operation to this
  /// @return reference to @c *this , for daisy-chaining
  NormalOperator &adjoint() {
    std::swap(creators_, annihilators_);
    std::for_each(creators_.begin(), creators_.end(), [](Op<S> &op) { op.adjoint(); });
    std::for_each(annihilators_.begin(), annihilators_.end(), [](Op<S> &op) { op.adjoint(); });
    return *this;
  }

 private:
  Vacuum vacuum_;
  Operator<S> creators_;
  Operator<S> annihilators_;
};

template<Statistics S>
bool operator==(const NormalOperator<S> &op1, const NormalOperator<S> &op2) {
  return op1.vacuum() == op2.vacuum() && op1.creators() == op2.creators() && op1.annihilators() == op2.annihilators();
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
  for (const auto &o : op.creators())
    result += to_latex(o.index());
  result += L"}_{";
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
