//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT_CORE_OP_H
#define SEQUANT_CORE_OP_H

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/hugenholtz.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/ranges.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/strong.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

#include <range/v3/all.hpp>

namespace sequant {

// strong type wrapper for objects associated with creation operators
DEFINE_STRONG_TYPES_FOR_RANGE_AND_RANGESIZE(cre);
// strong type wrapper for objects associated with annihilation operators
DEFINE_STRONG_TYPES_FOR_RANGE_AND_RANGESIZE(ann);

/// @brief Op is a creator/annihilator operator
///
/// Op = Index + Action
/// @tparam S specifies the particle statistics
template <Statistics S = Statistics::FermiDirac>
class Op {
 public:
  static constexpr Statistics statistics = S;

  Op() = default;
  Op(Index index, Action action) noexcept
      : index_(std::move(index)), action_(action) {}

  const Index &index() const { return index_; }
  Index &index() { return index_; }
  const Action &action() const { return action_; }

  /// @brief changes this to its (Hermitian) adjoint
  void adjoint() { action_ = sequant::adjoint(action_); }

  static std::wstring core_label() {
    return get_default_context(S).spbasis() == SPBasis::spinorbital
               ? (S == Statistics::FermiDirac ? L"a" : L"b")
               : L"E";
  }

  /// @return the string representation of @c this in LaTeX format
  std::wstring to_latex() const {
    std::wstring result;
    result = L"{";
    result += core_label();
    result += (action() == Action::create ? L"^{\\dagger}_" : L"_");
    result += index().to_latex();
    result += L"}";
    return result;
  }

  /// compares 2 Op objects while ignoring Index label (see Index::TypeCompare)
  struct TypeCompare {
    bool operator()(const Op<S> &first, const Op<S> &second) {
      bool result;
      if (first.action() == second.action()) {
        result = Index::TypeCompare{}(first.index(), second.index());
      } else
        result = first.action() < second.action();
      return result;
    }
  };

  /// tests equality of 2 Op objects while ignoring Index label (see
  /// Index::TypeEquality)
  struct TypeEquality {
    bool operator()(const Op<S> &first, const Op<S> &second) {
      bool result = (first.action() == second.action()) &&
                    Index::TypeEquality{}(first.index(), second.index());
      return result;
    }
  };

 private:
  Index index_;
  Action action_ = Action::invalid;
};

/// @brief The ordering operator

/// @return true if @c op1 preceeds @c op2 in the canonical order; Op objects
/// are ordered lexicographically, first by statistics (fermions < bosons)
/// then by action (cre < ann), then by their indices
template <Statistics S1, Statistics S2>
inline bool operator<(const Op<S1> &op1, const Op<S2> &op2) {
  if constexpr (S1 == S2) {
    if (op1.action() == op2.action()) {
      if (op1.index() == op2.index()) {
        return false;
      } else {
        return op1.index() < op2.index();
      }
    } else {
      return op1.action() < op2.action();
    }
  } else
    return S1 < S2;
}

/// @brief hashing function

/// @tparam S a Statistics value specifying the operator statistics
/// @param[in] op a const reference to an `Op<S>` object
/// @return the hash value of the object referred to by @c op
template <Statistics S>
inline auto hash_value(const Op<S> &op) {
  auto val = hash_value(op.index());
  hash::combine(val, op.action());
  hash::combine(val, S);
  return val;
}

template <Statistics S>
bool operator==(const Op<S> &op1, const Op<S> &op2) {
  return op1.index() == op2.index() && op1.action() == op2.action();
}

template <Statistics S>
bool operator!=(const Op<S> &op1, const Op<S> &op2) {
  return !(op1 == op2);
}

using BOp = Op<Statistics::BoseEinstein>;
using FOp = Op<Statistics::FermiDirac>;

inline BOp bcre(Index i) { return BOp(i, Action::create); }
template <typename I>
inline BOp bcre(Index i, std::initializer_list<I> pi) {
  return BOp(Index(i, pi), Action::create);
}
inline BOp bann(Index i) { return BOp(i, Action::annihilate); }
template <typename I>
inline BOp bann(Index i, std::initializer_list<I> pi) {
  return BOp(Index(i, pi), Action::annihilate);
}
inline FOp fcre(Index i) { return FOp(i, Action::create); }
template <typename I>
inline FOp fcre(Index i, std::initializer_list<I> pi) {
  return FOp(Index(i, pi), Action::create);
}
inline FOp fann(Index i) { return FOp(i, Action::annihilate); }
template <typename I>
inline FOp fann(Index i, std::initializer_list<I> pi) {
  return FOp(Index(i, pi), Action::annihilate);
}

/// @return true if this is a particle creator, false otherwise
template <Statistics S>
bool is_creator(const Op<S> &op) {
  return op.action() == Action::create;
};

/// @return true if this is a particle annihilator, false otherwise
template <Statistics S>
bool is_annihilator(const Op<S> &op) {
  return op.action() == Action::annihilate;
};

/// @return true if this is a pure quasiparticle creator with respect to the
/// given vacuum, false otherwise
template <Statistics S>
bool is_pure_qpcreator(const Op<S> &op,
                       Vacuum vacuum = get_default_context(S).vacuum()) {
  const auto &isr = get_default_context(S).index_space_registry();
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::create;
    case Vacuum::SingleProduct: {
      return (isr->is_pure_occupied(op.index().space()) &&
              op.action() == Action::annihilate) ||
             (isr->is_pure_unoccupied(op.index().space()) &&
              op.action() == Action::create);
    }
    default:
      throw std::logic_error(
          "is_pure_qpcreator: cannot handle MultiProduct vacuum");
  }
};

/// @return true if this is a quasiparticle creator with respect to the given
/// vacuum, false otherwise
template <Statistics S>
bool is_qpcreator(const Op<S> &op,
                  Vacuum vacuum = get_default_context(S).vacuum()) {
  const auto &isr = get_default_context(S).index_space_registry();
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::create;
    case Vacuum::SingleProduct: {
      return (isr->contains_occupied(op.index().space()) &&
              op.action() == Action::annihilate) ||
             (isr->contains_unoccupied(op.index().space()) &&
              op.action() == Action::create);
    }
    default:
      throw std::logic_error("is_qpcreator: cannot handle MultiProduct vacuum");
  }
};

template <Statistics S>
IndexSpace qpcreator_space(const Op<S> &op,
                           Vacuum vacuum = get_default_context(S).vacuum()) {
  const auto &isr = get_default_context(S).index_space_registry();
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::create ? op.index().space()
                                           : IndexSpace::null;
    case Vacuum::SingleProduct:
      return op.action() == Action::annihilate
                 ? isr->intersection(
                       op.index().space(),
                       isr->vacuum_occupied_space(op.index().space().qns()))
                 : isr->intersection(
                       op.index().space(),
                       isr->vacuum_unoccupied_space(op.index().space().qns()));
    default:
      throw std::logic_error(
          "qpcreator_space: cannot handle MultiProduct vacuum");
  }
}

/// @return true if this is a pure quasiparticle annihilator with respect to
/// the given vacuum, false otherwise
template <Statistics S>
bool is_pure_qpannihilator(const Op<S> &op,
                           Vacuum vacuum = get_default_context(S).vacuum()) {
  const auto &isr = get_default_context(S).index_space_registry();
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::annihilate;
    case Vacuum::SingleProduct: {
      return (isr->is_pure_unoccupied(op.index().space()) &&
              op.action() == Action::annihilate) ||
             (isr->is_pure_occupied(op.index().space()) &&
              op.action() == Action::create);
    }
    default:
      throw std::logic_error(
          "is_pure_qpannihilator: cannot handle MultiProduct vacuum");
  }
};

/// @return true if this is a quasiparticle annihilator with respect to the
/// given vacuum, false otherwise
template <Statistics S>
bool is_qpannihilator(const Op<S> &op,
                      Vacuum vacuum = get_default_context(S).vacuum()) {
  const auto &isr = get_default_context(S).index_space_registry();
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::annihilate;
    case Vacuum::SingleProduct: {
      return (isr->contains_occupied(op.index().space()) &&
              op.action() == Action::create) ||
             (isr->contains_unoccupied(op.index().space()) &&
              op.action() == Action::annihilate);
    }
    default:
      throw std::logic_error(
          "is_qpannihilator: cannot handle MultiProduct vacuum");
  }
};

template <Statistics S>
IndexSpace qpannihilator_space(
    const Op<S> &op, Vacuum vacuum = get_default_context(S).vacuum()) {
  const auto &isr = get_default_context(S).index_space_registry();
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::annihilate ? op.index().space()
                                               : IndexSpace::null;
    case Vacuum::SingleProduct:
      return op.action() == Action::create
                 ? isr->intersection(
                       op.index().space(),
                       isr->vacuum_occupied_space(op.index().space().qns()))
                 : isr->intersection(
                       op.index().space(),
                       isr->vacuum_unoccupied_space(op.index().space().qns()));
    default:
      throw std::logic_error(
          "qpcreator_space: cannot handle MultiProduct vacuum");
  }
}

template <Statistics S = Statistics::FermiDirac>
class NormalOperator;

/// @brief Operator is a sequence of Op objects
///
/// @tparam S specifies the particle statistics
template <Statistics S = Statistics::FermiDirac>
class Operator : public container::svector<Op<S>>, public Expr {
 public:
  using base_type = container::svector<Op<S>>;
  static constexpr Statistics statistics = S;

  // iterate over this using the base_type
  using iterator = typename base_type::iterator;
  using const_iterator = typename base_type::const_iterator;

  using base_type::at;
  using base_type::begin;
  using base_type::cbegin;
  using base_type::cend;
  using base_type::empty;
  using base_type::end;
  using base_type::size;
  using base_type::operator[];

  Operator() = default;
  explicit Operator(std::initializer_list<Op<S>> ops) : base_type(ops) {}
  explicit Operator(base_type &&ops) : base_type(std::move(ops)) {}
  template <typename I>
  Operator(Action action, std::initializer_list<I> indices)
      : base_type(make_ops(action, indices)) {}

  operator base_type &() const & { return *this; }
  operator base_type &&() && { return *this; }

  /// @brief adjoint of an Operator is a reversed string of the adjoints of its
  /// ops
  virtual void adjoint() override {
    std::reverse(this->begin(), this->end());
    std::for_each(this->begin(), this->end(), [](Op<S> &op) { op.adjoint(); });
    this->reset_hash_value();
  }

  /// @return the string representation of @c this in LaTeX format
  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    for (const auto &o : *this) result += o.to_latex();
    result += L"}";
    return result;
  }

  type_id_type type_id() const override { return get_type_id<Operator>(); };

  ExprPtr clone() const override { return std::make_shared<Operator>(*this); }

 private:
  base_type make_ops(Action action, IndexList indices) {
    base_type result;
    result.reserve(indices.size());
    for (const auto &idx : indices) result.emplace_back(idx, action);
    return result;
  }

  base_type make_ops(Action action, WstrList index_labels) {
    base_type result;
    result.reserve(index_labels.size());
    for (const auto &idx_label : index_labels)
      result.emplace_back(Index{idx_label}, action);
    return result;
  }

  bool static_equal(const Expr &that) const override;

  bool is_cnumber() const override { return false; }

  bool commutes_with_atom(const Expr &that) const override {
    bool result = true;
    // does not commute with Operator<S>
    // TODO implement checks of commutativity with Operator<S>
    if (that.is<Operator<S>>()) {
      result = false;
    } else if (that.is<NormalOperator<S>>()) {
      result = that.as<NormalOperator<S>>().commutes_with_atom(*this);
    }
    return result;
  }

  hash_type memoizing_hash() const override {
    using std::begin;
    using std::end;
    const auto &ops = static_cast<const base_type &>(*this);
    return hash::range(begin(ops), end(ops));
  }
};

template <Statistics S>
inline bool operator==(const Operator<S> &one, const Operator<S> &another) {
  using base_type = container::svector<Op<S>>;
  if (one.size() == another.size()) {
    if (one.empty()) return true;
    // compare hash values first
    if (one.hash_value() ==
        another.hash_value())  // hash values agree -> do full comparison
      return static_cast<const base_type &>(one) ==
             static_cast<const base_type &>(another);
    else
      return false;
  } else
    return false;
}

template <Statistics S>
bool Operator<S>::static_equal(const Expr &that) const {
  const auto &that_cast = static_cast<const Operator &>(that);
  return *this == that_cast;
}

/// @brief NormalOperator is an Operator normal-ordered with respect to vacuum.

/// @note Normal ordering means all creators are to the left of all
/// annihilators. It is natural to express at least number-conserving normal
/// operators (i.e. those with equal number of creators and annihilators) as
/// tensors with creators as superscripts and annihilators as subscripts.
/// Operator cre(p1) cre(p2) ... cre(pN) ann(qN) ... ann(q2) ann(q1) is
/// represented in such notation as a^{p1 p2 ... pN}_{q1 q2 ... qN}, hence it is
/// natural to specify annihilators in the order of their particle index (i.e.
/// as q1 q2, etc.) which is reverse of the order of their appearance in
/// Operator.
///
/// @note The tensor notation becomes less intuitive for number non-conserving
/// operators, e.g. cre(p1) cre(p2) ann(q2) could be represented as a^{p1
/// p2}_{q2 ⎵} or a^{p1 p2}_{⎵ q2}. To make it explicit that ann(q2) acts on
/// same particle as cre(p2) the latter notation is used; similarly, cre(p1)
/// ann(q1) ann(q2) is represented as a^{⎵ p1}_{q1 q2}.
///
/// @tparam S specifies the particle statistics
template <Statistics S>
class NormalOperator : public Operator<S>,
                       public AbstractTensor,
                       public Labeled {
 public:
  static constexpr Statistics statistics = S;
  using base_type = Operator<S>;
  using vector_type = typename Operator<S>::base_type;

  // iterate over this using the base_type
  using iterator = typename Operator<S>::iterator;
  using const_iterator = typename Operator<S>::const_iterator;

  using base_type::at;
  using base_type::begin;
  using base_type::cbegin;
  using base_type::cend;
  using base_type::empty;
  using base_type::end;
  using base_type::size;
  using base_type::operator[];

  /// constructs an identity operator
  NormalOperator(Vacuum v = get_default_context(S).vacuum()) {}

  /// @tparam IndexOrOpSequence1 type representing a sequence of objects that
  /// can be statically cast into Index or Op<S>
  /// @tparam IndexOrOpSequence2 type representing a sequence of objects that
  /// can be statically cast into Index or Op<S>
  /// @param creators sequence of creator indices or operators (in order of
  /// particle indices)
  /// @param annihilators sequence of annihilator indices or operators (in order
  /// of particle indices).
  /// @param v vacuum state with respect to which the operator is normal-ordered
  template <
      typename IndexOrOpSequence1, typename IndexOrOpSequence2,
      typename = std::enable_if_t<
          (meta::is_statically_castable_v<
               meta::range_value_t<IndexOrOpSequence1>, Index> ||
           meta::is_statically_castable_v<
               meta::range_value_t<IndexOrOpSequence1>,
               Op<S>>)&&(meta::
                             is_statically_castable_v<
                                 meta::range_value_t<IndexOrOpSequence2>,
                                 Index> ||
                         meta::is_statically_castable_v<
                             meta::range_value_t<IndexOrOpSequence2>, Op<S>>)>>
  NormalOperator(const cre<IndexOrOpSequence1> &creators,
                 const ann<IndexOrOpSequence2> &annihilators,
                 Vacuum v = get_default_context(S).vacuum())
      : Operator<S>{}, vacuum_(v), ncreators_(ranges::size(creators)) {
    this->reserve(ranges::size(creators) + ranges::size(annihilators));
    for (const auto &c : creators) {
      assert((!std::is_same_v<meta::remove_cvref_t<IndexOrOpSequence1>,
                              std::array<meta::castable_to_any, 0>>));
      if constexpr (meta::is_statically_castable_v<
                        meta::range_value_t<IndexOrOpSequence1>, Index>)
        this->emplace_back(c, Action::create);
      else {
        assert(c.action() == Action::create);
        this->emplace_back(c);
      }
    }
    for (const auto &a : annihilators | ranges::views::reverse) {
      assert((!std::is_same_v<meta::remove_cvref_t<IndexOrOpSequence2>,
                              std::array<meta::castable_to_any, 0>>));
      if constexpr (meta::is_statically_castable_v<
                        meta::range_value_t<IndexOrOpSequence2>, Index>)
        this->emplace_back(a, Action::annihilate);
      else {
        assert(a.action() == Action::annihilate);
        this->emplace_back(a);
      }
    }
  }

  NormalOperator(const NormalOperator &other)
      : Operator<S>(other),
        vacuum_(other.vacuum_),
        ncreators_(other.ncreators_),
        hug_(other.hug_ ? std::make_unique<hug_type>(*other.hug_) : nullptr) {}
  NormalOperator(NormalOperator &&) = default;
  NormalOperator &operator=(NormalOperator &&) = default;
  NormalOperator &operator=(const NormalOperator &other) {
    static_cast<base_type &>(*this) = static_cast<const base_type &>(other);
    vacuum_ = other.vacuum_;
    ncreators_ = other.ncreators_;
    hug_ = other.hug_ ? std::make_unique<hug_type>(*other.hug_) : nullptr;
    return *this;
  }

  /// @return the vacuum state with respect to which the operator is
  /// normal-ordered.
  Vacuum vacuum() const { return vacuum_; }
  /// @return the range of creators, in the order of increasing particle index
  auto creators() const {
    return ranges::views::counted(this->cbegin(), ncreators());
  }
  /// @return the range of annihilators, in the order of increasing particle
  /// index
  auto annihilators() const {
    return ranges::views::counted(this->crbegin(), nannihilators());
  }
  /// @return the number of creators
  auto ncreators() const { return ncreators_; }
  /// @return the number of annihilators
  auto nannihilators() const { return this->size() - ncreators(); }
  /// @return view of creators and annihilators as a single range
  auto creann() const {
    return ranges::views::concat(creators(), annihilators());
  }

  /// @return number of creators/annihilators
  /// @throw std::logic_error if the operator is not particle number conserving
  /// (i.e. if ncreators() != nannihilators() )
  auto rank() const {
    if (ncreators() != nannihilators()) {
      throw std::logic_error(
          "NormalOperator::rank(): ncreators != nannihilators");
    }
    return ncreators();
  }

  /// @return the representation of @c *this as a Hugenholtz vertex
  /// @sa HugenholtzVertex
  const auto &hug() const {
    if (!hug_) {
      hug_ = std::make_unique<hug_type>(
          static_cast<const typename base_type::base_type &>(*this));
    }
    return hug_;
  }

  /// @return all possible values returned by label() for this operator type
  static const container::svector<std::wstring> &labels();

  std::wstring_view label() const override;

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    if (vacuum() == Vacuum::Physical) {
      result += Op<S>::core_label();
    } else {
      result += L"\\tilde{";
      result += Op<S>::core_label();
      result += L"}";
    }
    const auto ncreators = this->ncreators();
    const auto nannihilators = this->nannihilators();
    if (ncreators > 0) {
      result += L"^{";
      if (ncreators <
          nannihilators) {  // if have more annihilators than creators pad on
                            // the left with square underbrackets, i.e. ⎵
        const auto iend = nannihilators - ncreators;
        if (iend > 0) result += L"\\textvisiblespace";
        for (size_t i = 1; i != iend; ++i) {
          result += L"\\,\\textvisiblespace";
        }
        result += L"\\,";
      }
      for (const auto &o : creators()) result += o.index().to_latex();
      result += L"}";
    }
    if (nannihilators > 0) {
      result += L"_{";
      if (ncreators >
          nannihilators) {  // if have more creators than annihilators pad on
                            // the left with square underbrackets, i.e. ⎵
        const auto iend = ncreators - nannihilators;
        if (iend > 0) result += L"\\textvisiblespace";
        for (size_t i = 1; i != iend; ++i) {
          result += L"\\,\\textvisiblespace";
        }
        result += L"\\,";
      }
      for (const auto &o : annihilators()) result += o.index().to_latex();
      result += L"}";
    }
    result += L"}";
    return result;
  }

  /// overload base_type::erase
  iterator erase(const_iterator it) {
    if (it->action() == Action::create) --ncreators_;
    if (hug_) hug_->erase(it - begin(), *it);
    return Operator<S>::erase(it);
  }

  /// overload base_type::insert
  template <typename T>
  iterator insert(const_iterator it, T &&value) {
    if (value.action() == Action::create) ++ncreators_;
    auto result = Operator<S>::insert(it, std::forward<T>(value));
    if (hug_) hug_->insert(result - begin(), *result);
    return result;
  }

  Expr::type_id_type type_id() const override {
    return Expr::get_type_id<NormalOperator>();
  };

  ExprPtr clone() const override {
    return std::make_shared<NormalOperator>(*this);
  }

  virtual void adjoint() override {
    // same as base adjoint(), but updates extra state
    Operator<S>::adjoint();
    hug_.reset();
    ncreators_ = nannihilators();
  }

  /// Replaces indices using the index map
  /// @param index_map maps Index to Index
  /// @return true if one or more indices changed
  /// @pre indices are not tagged, or (if want to protect them from replacement)
  /// tagged with (int)0
  /// @post indices that were replaced will be tagged with (int)0
  template <template <typename, typename, typename... Args> class Map,
            typename... Args>
  bool transform_indices(const Map<Index, Index, Args...> &index_map) {
    bool mutated = false;
    ranges::for_each(*this, [&](auto &&op) {
      if (op.index().transform(index_map)) mutated = true;
    });
    if (mutated) this->reset_hash_value();
    return mutated;
  }

 private:
  Vacuum vacuum_;
  std::size_t ncreators_ = 0;
  using hug_type = HugenholtzVertex<Op<S>, typename Op<S>::TypeEquality>;
  mutable std::unique_ptr<hug_type> hug_;  // only created if needed

  bool static_equal(const Expr &that) const override;

  bool static_less_than(const Expr &that) const override {
    auto range_hash = [](const auto &sized_range) {
      using ranges::begin;
      using ranges::size;
      auto b = begin(sized_range);
      auto e = b + size(sized_range);
      auto val = hash::range(b, e);
      return val;
    };
    auto range_compare = [](const auto &sized_range1,
                            const auto &sized_range2) {
      using ranges::begin;
      using ranges::size;
      auto b1 = begin(sized_range1);
      auto e1 = b1 + size(sized_range1);
      auto b2 = begin(sized_range2);
      auto e2 = b2 + size(sized_range2);
      auto val = std::lexicographical_compare(b1, e1, b2, e2);
      return val;
    };

    const auto &that_cast = static_cast<const NormalOperator<S> &>(that);
    if (this == &that) return false;
    if (this->ncreators() == that_cast.ncreators()) {
      if (this->nannihilators() == that_cast.nannihilators()) {
        // unlike Tensor comparison, we don't memoize hashes of creators and
        // annihilators separately
        auto cre_hash = range_hash(this->creators());
        auto that_cre_hash = range_hash(that_cast.creators());
        if (cre_hash == that_cre_hash) {
          auto ann_hash = range_hash(this->annihilators());
          auto that_ann_hash = range_hash(that_cast.annihilators());
          if (ann_hash == that_ann_hash)
            return false;
          else {
            return range_compare(this->annihilators(),
                                 that_cast.annihilators());
          }
        } else {
          return range_compare(this->creators(), that_cast.creators());
        }
      } else {
        return this->nannihilators() < that_cast.nannihilators();
      }
    } else {
      return this->ncreators() < that_cast.ncreators();
    }
  }

  template <Statistics>
  friend class Operator;
  bool commutes_with_atom(const Expr &that) const override {
    const auto &isr = get_default_context(S).index_space_registry();
    // same as WickTheorem::can_contract
    auto can_contract = [this, &isr](const Op<S> &left, const Op<S> &right) {
      if (is_qpannihilator<S>(left, vacuum_) &&
          is_qpcreator<S>(right, vacuum_)) {
        const auto qpspace_left = qpannihilator_space<S>(left, vacuum_);
        const auto qpspace_right = qpcreator_space<S>(right, vacuum_);
        const auto qpspace_common =
            isr->intersection(qpspace_left, qpspace_right);
        if (qpspace_common) return true;
      }
      return false;
    };

    bool result = true;
    /// does not commute with Operator<S>
    /// TODO implement checks of commutativity with Operator<S>
    if (that.is<Operator<S>>()) {
      result = false;
    } else if (that.is<NormalOperator<S>>()) {
      const auto &op_that = that.as<NormalOperator<S>>();
      if (vacuum() ==
          op_that.vacuum()) {  // can only check commutativity for same vacua
        for (auto &&op_l : *this) {
          for (auto &&op_r : op_that) {
            if (can_contract(op_l, op_r)) {
              return false;
            }
          }
        }
      } else {  // NormalOperators for different vacua do not commute
        result = false;
      }
    }
    return result;
  }

  // these implement the AbstractTensor interface
  AbstractTensor::const_any_view_randsz _bra() const override final {
    return annihilators() |
           ranges::views::transform(
               [](auto &&op) -> const Index & { return op.index(); });
  }
  AbstractTensor::const_any_view_randsz _ket() const override final {
    return creators() |
           ranges::views::transform(
               [](auto &&op) -> const Index & { return op.index(); });
  }
  AbstractTensor::const_any_view_rand _braket() const override final {
    return ranges::views::concat(annihilators(), creators()) |
           ranges::views::transform(
               [](auto &&op) -> const Index & { return op.index(); });
  }
  std::size_t _bra_rank() const override final { return nannihilators(); }
  std::size_t _ket_rank() const override final { return ncreators(); }
  Symmetry _symmetry() const override final {
    return (S == Statistics::FermiDirac
                ? (get_default_context(S).spbasis() == SPBasis::spinorbital
                       ? Symmetry::antisymm
                       : Symmetry::nonsymm)
                : (Symmetry::symm));
  }
  BraKetSymmetry _braket_symmetry() const override final {
    return BraKetSymmetry::nonsymm;
  }
  ParticleSymmetry _particle_symmetry() const override final {
    return ParticleSymmetry::symm;
  }
  std::size_t _color() const override final {
    return S == Statistics::FermiDirac ? 1 : 2;
  }
  bool _is_cnumber() const override final { return false; }
  std::wstring_view _label() const override final { return label(); }
  std::wstring _to_latex() const override final { return to_latex(); }
  bool _transform_indices(
      const container::map<Index, Index> &index_map) override final {
    return transform_indices(index_map);
  }
  void _reset_tags() override final {
    ranges::for_each(*this, [](const auto &op) { op.index().reset_tag(); });
  }
  bool operator<(const AbstractTensor &other) const override final {
    auto *other_nop = dynamic_cast<const NormalOperator<S> *>(&other);
    if (other_nop) {
      const Expr *other_expr = static_cast<const Expr *>(other_nop);
      return this->static_less_than(*other_expr);
    } else
      return false;  // TODO do we compare typeid? labels? probably the latter
  }

  AbstractTensor::any_view_randsz _bra_mutable() override final {
    this->reset_hash_value();
    return ranges::views::counted(this->rbegin(), nannihilators()) |
           ranges::views::transform(
               [](auto &&op) -> Index & { return op.index(); });
  }
  AbstractTensor::any_view_randsz _ket_mutable() override final {
    this->reset_hash_value();
    return ranges::views::counted(this->begin(), ncreators()) |
           ranges::views::transform(
               [](auto &&op) -> Index & { return op.index(); });
  }
};

template <Statistics S>
bool operator==(const NormalOperator<S> &op1, const NormalOperator<S> &op2) {
  using base_type = Operator<S>;
  if (op1.vacuum() == op2.vacuum() && op1.ncreators() == op2.ncreators()) {
    return static_cast<const base_type &>(op1) ==
           static_cast<const base_type &>(op2);
  } else
    return false;
}

template <Statistics S>
bool NormalOperator<S>::static_equal(const Expr &that) const {
  const auto &that_cast = static_cast<const NormalOperator &>(that);
  return *this == that_cast;
}

/// @brief NormalOperatorSequence is a sequence NormalOperator objects, all
/// ordered with respect to same vacuum
///
/// @tparam S specifies the particle statistics
template <Statistics S = Statistics::FermiDirac>
class NormalOperatorSequence : public container::svector<NormalOperator<S>>,
                               public Expr {
 public:
  using base_type = container::svector<NormalOperator<S>>;
  static constexpr Statistics statistics = S;

  using base_type::at;
  using base_type::begin;
  using base_type::cbegin;
  using base_type::cend;
  using base_type::empty;
  using base_type::end;
  using base_type::size;
  using base_type::operator[];

  /// constructs an empty sequence
  NormalOperatorSequence() : vacuum_(get_default_context(S).vacuum()) {}

  /// constructs from a parameter pack
  template <typename... NOps,
            typename = std::enable_if_t<
                (std::is_convertible_v<NOps, NormalOperator<S>> && ...)>>
  NormalOperatorSequence(NOps &&...operators)
      : base_type({std::forward<NOps>(operators)...}) {
    check_vacuum();
  }

  /// constructs from an initializer list
  NormalOperatorSequence(std::initializer_list<NormalOperator<S>> operators)
      : base_type(operators) {
    check_vacuum();
  }

  Vacuum vacuum() const { return vacuum_; }

  operator const base_type &() const & { return *this; }
  operator base_type &&() && { return *this; }

  /// @return the total number of `Op<S>` objects in this
  /// @warning not to be confused with `NormalOperatorSequence<S>::size()` that
  /// returns the number of `NormalOperator<S>` objects
  auto opsize() const {
    size_t opsz = 0;
    for (auto &&nop : *this) {
      opsz += nop.size();
    }
    return opsz;
  }

  /// @brief adjoint of a NormalOperatorSequence is a reversed sequence of
  /// adjoints
  virtual void adjoint() override {
    std::reverse(this->begin(), this->end());
    std::for_each(this->begin(), this->end(),
                  [](NormalOperator<S> &op) { op.adjoint(); });
    reset_hash_value();
  }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    for (const auto &op : *this) result += op.to_latex();
    result += L"}";
    return result;
  }

  Expr::type_id_type type_id() const override {
    return Expr::get_type_id<NormalOperatorSequence>();
  };

 private:
  Vacuum vacuum_ = Vacuum::Invalid;
  /// ensures that all operators use same vacuum, and sets vacuum_
  void check_vacuum() {
    vacuum_ = std::accumulate(
        this->cbegin(), this->cend(), Vacuum::Invalid,
        [](Vacuum v1, const NormalOperator<S> &v2) {
          if (v1 == Vacuum::Invalid) {
            return v2.vacuum();
          } else {
            if (v1 != v2.vacuum())
              throw std::invalid_argument(
                  "NormalOperatorSequence expects all constituent "
                  "NormalOperator objects to use same vacuum");
            else
              return v1;
          }
        });
  }

  bool static_equal(const Expr &that) const override {
    const auto &that_cast = static_cast<const NormalOperatorSequence &>(that);
    if (this->vacuum() == that_cast.vacuum()) {
      if (this->empty()) return true;
      if (this->hash_value() == that.hash_value())
        return static_cast<const base_type &>(*this) ==
               static_cast<const base_type &>(*this);
      else
        return false;
    } else
      return false;
  }
};

using BOperator = Operator<Statistics::BoseEinstein>;
using BNOperator = NormalOperator<Statistics::BoseEinstein>;
using BNOperatorSeq = NormalOperatorSequence<Statistics::BoseEinstein>;
using FOperator = Operator<Statistics::FermiDirac>;
using FNOperator = NormalOperator<Statistics::FermiDirac>;
using FNOperatorSeq = NormalOperatorSequence<Statistics::FermiDirac>;

template <typename... Attr>
inline ExprPtr bcrex(Index i, Attr &&...attr) {
  return ex<BNOperator>(cre({bcre(i, std::forward<Attr>(attr)...)}), ann());
}
template <typename... Attr>
inline ExprPtr bannx(Index i, Attr &&...attr) {
  return ex<BNOperator>(cre({}), ann({bann(i, std::forward<Attr>(attr)...)}));
}
template <typename... Attr>
inline ExprPtr fcrex(Index i, Attr &&...attr) {
  return ex<FNOperator>(cre({fcre(i, std::forward<Attr>(attr)...)}), ann());
}
template <typename... Attr>
inline ExprPtr fannx(Index i, Attr &&...attr) {
  return ex<FNOperator>(cre(), ann({fann(i, std::forward<Attr>(attr)...)}));
}

template <Statistics S>
std::wstring to_latex(const NormalOperator<S> &op) {
  return op.to_latex();
}

template <Statistics S>
std::wstring to_latex(const NormalOperatorSequence<S> &opseq) {
  return opseq.to_latex();
}

namespace detail {
struct OpIdRegistrar {
  OpIdRegistrar();
};
}  // namespace detail

/// converts NormalOperatorSequence to NormalOperator
/// @tparam S Statistics
/// @param[in] opseq a `NormalOperatorSequence<S>` object
/// @param[in] target_partner_indices ptr to sequence of Index pairs whose
/// `Op<S>` will act on same particle, if possible; if null, will not be used
/// @return @c {phase,normal_operator} , where @c phase is +1 or -1, and @c
/// normal_operator is a `NormalOperator<S>` object
/// @note will try to ensure that `Op<S>` objects for each there is a pairs of
/// Indices in @p target_index_columns will act on the same particle in the
/// result
template <Statistics S = Statistics::FermiDirac>
std::tuple<int, std::shared_ptr<NormalOperator<S>>> normalize(
    const NormalOperatorSequence<S> &opseq,
    const container::svector<std::pair<Index, Index>> &target_partner_indices =
        {}) {
  int phase = 1;
  container::svector<Op<S>> creators, annihilators;

  using opseq_view_type = flattened_rangenest<const NormalOperatorSequence<S>>;
  auto opseq_view = opseq_view_type(&opseq);
  using std::begin;
  using std::end;

  auto opseq_view_iter = begin(opseq_view);
  auto opseq_view_end = end(opseq_view);
  std::size_t pos = 0;
  const auto nops = opseq.opsize();  // # of Op<S>
  const auto nnops = opseq.size();   // # of NormalOperator<S>
  assert(nnops > 0);
  auto vacuum = opseq.vacuum();
  container::svector<container::svector<Op<S>>> annihilator_groups(nnops);
  while (opseq_view_iter != opseq_view_end) {
    // creators: since they go left, concat them to preserve the order of
    // NormalOperators and creators within Operators ...
    if (is_creator(*opseq_view_iter)) {
      creators.push_back(*opseq_view_iter);
    }  // annihilators: rev-concat the groups corresponding to NormalOperators
       // to obtain the desired order of normalized NormalOperators, but within
       // each group simply concat to preserve the original order
    else {
      assert(is_annihilator(*opseq_view_iter));

      // distance from the given Op to the end of NOpSeq (to rev-concat NOps we
      // are pushing to the first available group ... since we are using vector
      // we are filling Op groups from the back, but this does not create any
      // extra phase)
      const auto distance_to_end = nops - pos - 1;
      const auto op_group =
          nnops - ranges::get_cursor(opseq_view_iter).range_ordinal() -
          1;  // group idx = reverse of the NOp index within NopSeq
      assert(op_group < nnops);
      const auto total_distance =
          distance_to_end + annihilator_groups[op_group]
                                .size();  // we are fwd-concating Ops within
                                          // groups, hence extra phase
      annihilator_groups[op_group].push_back(*opseq_view_iter);
      if (S == Statistics::FermiDirac && total_distance % 2 == 1) {
        phase *= -1;
      }
    }
    ++pos;
    ++opseq_view_iter;
  }

  // convert annihilator_groups to annihilators
  // N.B. the NormalOperator ctor expects annihilators in reverse order!
  ranges::for_each(annihilator_groups | ranges::views::reverse,
                   [&annihilators](const auto &agrp) {
                     ranges::for_each(agrp | ranges::views::reverse,
                                      [&annihilators](const auto &op) {
                                        annihilators.push_back(op);
                                      });
                   });

  // if particle_ops given, reorder to preserve the "original" pairs of Op<S>,
  // if possible. Max number of such pairs (annihilator/creator columns,
  // in tensor notation) = min(#cre,#ann)
  if (!target_partner_indices.empty()) {
    const auto ncre = ranges::size(creators);
    const auto nann = ranges::size(annihilators);
    const auto rank = std::min(ncre, nann);

    // for every creator/annihilator, in reverse/forward order, if it is in
    // target_index_columns, locate matching annihilator/creator, if any, and
    // place in the
    if (ncre == rank) {
      const auto nann_extra = nann - rank;
      for (std::int64_t p = rank - 1; p >= 0;
           --p) {  // "reverse" particle index
        const auto &cre_op = *(creators.begin() + p);
        auto cre_it_in_index_cols = ranges::find_if(
            target_partner_indices, [&cre_op](const auto &idx_pair) {
              return idx_pair.first == cre_op.index();
            });
        if (cre_it_in_index_cols != target_partner_indices.end()) {
          const auto &ann_idx = cre_it_in_index_cols->second;
          auto source_ann_it = ranges::find_if(
              annihilators,
              [&ann_idx](const auto &v) { return v.index() == ann_idx; });
          if (source_ann_it != annihilators.end()) {
            assert(p + nann_extra < annihilators.size());
            const auto target_ann_it = annihilators.begin() + p + nann_extra;
            if (target_ann_it != source_ann_it) {
              if (S == Statistics::FermiDirac) phase *= -1.;  // 1 swap
              std::swap(*source_ann_it, *target_ann_it);
            }
          }
        }
      }
    }       // ncre == rank
    else {  // nann == rank
      assert(nann == rank);
      const auto ncre_extra = ncre - rank;
      for (std::int64_t p = rank - 1; p >= 0;
           --p) {  // "reverse" particle index
        const auto &ann_op = *(annihilators.begin() + p);
        auto ann_it_in_index_cols = ranges::find_if(
            target_partner_indices, [&ann_op](const auto &idx_pair) {
              return idx_pair.second == ann_op.index();
            });
        if (ann_it_in_index_cols != target_partner_indices.end()) {
          const auto &cre_idx = ann_it_in_index_cols->first;
          auto source_cre_it = ranges::find_if(
              creators,
              [&cre_idx](const auto &v) { return v.index() == cre_idx; });
          if (source_cre_it != creators.end()) {
            assert(p + ncre_extra < creators.size());
            const auto target_cre_it = creators.begin() + p + ncre_extra;
            if (source_cre_it != target_cre_it) {
              if (S == Statistics::FermiDirac) phase *= -1.;  // 1 swap
              std::swap(*source_cre_it, *target_cre_it);
            }
          }
        }
      }
    }  // nann == rank
  }

  return std::make_tuple(phase, std::make_shared<NormalOperator<S>>(
                                    cre(std::move(creators)),
                                    ann(std::move(annihilators)), vacuum));
}

template <typename T>
std::decay_t<T> adjoint(
    T &&t, std::void_t<decltype(std::declval<T &>().adjoint())> * = nullptr) {
  if constexpr (std::is_reference_v<T>) {
    std::decay_t<T> t_copy(t);
    t_copy.adjoint();
    return t_copy;
  } else {
    t.adjoint();
    return std::move(t);
  }
}

}  // namespace sequant

#endif  // SEQUANT_CORE_OP_H
