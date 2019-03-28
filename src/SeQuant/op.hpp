//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT_OP_H
#define SEQUANT_OP_H

#include <numeric>

#include <range/v3/all.hpp>

#include <boost/functional/hash_fwd.hpp>

#include "abstract_tensor.hpp"
#include "attr.hpp"
#include "container.hpp"
#include "expr.hpp"
#include "hugenholtz.hpp"
#include "index.hpp"
#include "ranges.hpp"
#include "sequant.hpp"

namespace sequant {

/// @brief Op is a creator/annihilator operator
///
/// Op = Index + Action
/// @tparam S specifies the particle statistics
template<Statistics S = Statistics::FermiDirac>
class Op {
 public:
  static constexpr Statistics statistics = S;

  Op() = default;
  Op(Index index, Action action) noexcept : index_(std::move(index)), action_(action) {}

  const Index &index() const { return index_; }
  Index &index() { return index_; }
  const Action &action() const { return action_; }
  /// applies (Hermitian) adjoint to this operator
  Op &adjoint() {
    action_ = sequant::adjoint(action_);
    return *this;
  }

  /// @return the string representation of @c this in LaTeX format
  std::wstring to_latex() const {
    std::wstring result;
    result = L"{";
    result += (S == Statistics::FermiDirac ? L"a" : L"b");
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
      }
      else {
        return op1.index() < op2.index();
      }
    } else {
      return op1.action() < op2.action();
    }
  }
  else
    return S1 < S2;
}

/// @brief hashing function

/// @tparam S a Statistics value specifying the operator statistics
/// @paramp[in] op a const reference to an Op<S> object
/// @return the hash value of the object referred to by @c op
/// @note uses boost::hash_value
template <Statistics S>
inline auto hash_value(const Op<S> &op) {
  auto val = hash_value(op.index());
  boost::hash_combine(val, op.action());
  return val;
}

template<Statistics S>
bool operator==(const Op<S> &op1, const Op<S> &op2) {
  return op1.index() == op2.index() && op1.action() == op2.action();
}

template<Statistics S>
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

/// @return true if this is a pure quasdiparticle creator with respect to the
/// given vacuum, false otherwise
template <Statistics S>
bool is_pure_qpcreator(const Op<S> &op,
                       Vacuum vacuum = get_default_context().vacuum()) {
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::create;
    case Vacuum::SingleProduct: {
      const auto occ_class = occupancy_class(op.index().space());
      return (occ_class < 0 && op.action() == Action::annihilate) ||
             (occ_class > 0 && op.action() == Action::create);
    }
    default:
      throw std::logic_error(
          "is_pure_qpcreator: cannot handle MultiProduct vacuum");
  }
};

/// @return true if this is a quasdiparticle creator with respect to the given
/// vacuum, false otherwise
template <Statistics S>
bool is_qpcreator(const Op<S> &op,
                  Vacuum vacuum = get_default_context().vacuum()) {
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::create;
    case Vacuum::SingleProduct: {
      const auto occ_class = occupancy_class(op.index().space());
      return (occ_class <= 0 && op.action() == Action::annihilate) ||
             (occ_class >= 0 && op.action() == Action::create);
    }
    default:
      throw std::logic_error("is_qpcreator: cannot handle MultiProduct vacuum");
  }
};

template <Statistics S>
IndexSpace qpcreator_space(const Op<S> &op,
                           Vacuum vacuum = get_default_context().vacuum()) {
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::create ? op.index().space()
                                           : IndexSpace::null_instance();
    case Vacuum::SingleProduct:
      return op.action() == Action::annihilate
                 ? intersection(op.index().space(),
                                IndexSpace::instance(IndexSpace::occupied))
                 : intersection(
                       op.index().space(),
                       IndexSpace::instance(IndexSpace::complete_unoccupied));
    default:
      throw std::logic_error(
          "qpcreator_space: cannot handle MultiProduct vacuum");
  }
}

/// @return true if this is a pure quasdiparticle annihilator with respect to
/// the given vacuum, false otherwise
template <Statistics S>
bool is_pure_qpannihilator(const Op<S> &op,
                           Vacuum vacuum = get_default_context().vacuum()) {
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::annihilate;
    case Vacuum::SingleProduct: {
      const auto occ_class = occupancy_class(op.index().space());
      return (occ_class > 0 && op.action() == Action::annihilate) ||
             (occ_class < 0 && op.action() == Action::create);
    }
    default:
      throw std::logic_error(
          "is_pure_qpannihilator: cannot handle MultiProduct vacuum");
  }
};

/// @return true if this is a quasdiparticle annihilator with respect to the
/// given vacuum, false otherwise
template <Statistics S>
bool is_qpannihilator(const Op<S> &op,
                      Vacuum vacuum = get_default_context().vacuum()) {
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::annihilate;
    case Vacuum::SingleProduct: {
      const auto occ_class = occupancy_class(op.index().space());
      return (occ_class >= 0 && op.action() == Action::annihilate) ||
             (occ_class <= 0 && op.action() == Action::create);
    }
    default:
      throw std::logic_error(
          "is_qpannihilator: cannot handle MultiProduct vacuum");
  }
};

template <Statistics S>
IndexSpace qpannihilator_space(const Op<S> &op,
                               Vacuum vacuum = get_default_context().vacuum()) {
  switch (vacuum) {
    case Vacuum::Physical:
      return op.action() == Action::annihilate ? op.index().space()
                                               : IndexSpace::null_instance();
    case Vacuum::SingleProduct:
      return op.action() == Action::create
                 ? intersection(op.index().space(),
                                IndexSpace::instance(IndexSpace::occupied))
                 : intersection(
                       op.index().space(),
                       IndexSpace::instance(IndexSpace::complete_unoccupied));
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
template<Statistics S = Statistics::FermiDirac>
class Operator : public container::svector<Op<S>>, public Expr {
 public:
  using base_type = container::svector<Op<S>>;
  static constexpr Statistics statistics = S;

  // iterate over this using the base_type
  using iterator = typename base_type::iterator;
  using const_iterator = typename base_type::const_iterator;

  using base_type::begin;
  using base_type::end;
  using base_type::cbegin;
  using base_type::cend;
  using base_type::empty;
  using base_type::size;
  using base_type::at;
  using base_type::operator[];

  Operator() = default;
  explicit Operator(std::initializer_list<Op<S>> ops)
      : base_type(ops) {}
  explicit Operator(base_type &&ops)
      : base_type(std::move(ops)) {}
  template <typename I>
  Operator(Action action, std::initializer_list<I> indices)
      : base_type(make_ops(action, indices)) {}

  operator base_type &() const & { return *this; }
  operator base_type &&() && { return *this; }

  /// applies (Hermitian) adjoint operation to this
  /// @return reference to @c *this , for daisy-chaining
  Operator &adjoint() {
    std::reverse(this->begin(), this->end());
    std::for_each(this->begin(), this->end(), [](Op<S> &op) { op.adjoint(); });
    return *this;
  }

  /// @return the string representation of @c this in LaTeX format
  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    for (const auto &o : *this)
      result += o.to_latex();
    result += L"}";
    return result;
  }

  type_id_type type_id() const override {
    return get_type_id<Operator>();
  };

  std::shared_ptr<Expr> clone() const override {
    return std::make_shared<Operator>(*this);
  }

 private:
  base_type make_ops(Action action,
                     IndexList indices) {
    base_type result;
    result.reserve(indices.size());
    for (const auto &idx : indices)
      result.emplace_back(idx, action);
    return result;
  }

  base_type
  make_ops(Action action,
           WstrList index_labels) {
    base_type result;
    result.reserve(index_labels.size());
    for (const auto &idx_label : index_labels)
      result.emplace_back(Index{idx_label}, action);
    return result;
  }

  bool static_equal(const Expr &that) const override {
    const auto &that_cast = static_cast<const Operator &>(that);
    if (this->size() == that_cast.size()) {
      if (this->empty()) return true;
      // compare hash values first
      if (this->hash_value() ==
          that.hash_value())  // hash values agree -> do full comparison
        return static_cast<const base_type &>(*this) ==
               static_cast<const base_type &>(that_cast);
      else
        return false;
    } else
      return false;
  }

  bool is_cnumber() const override {
    return false;
  }

  bool commutes_with_atom(const Expr& that) const override {
    bool result = true;
    /// does not commute with Operator<S>
    /// TODO implement checks of commutativity with Operator<S>
    if (that.is<Operator<S>>()) {
      result = false;
    }
    else if (that.is<NormalOperator<S>>()) {
      result = that.as<NormalOperator<S>>().commutes_with_atom(*this);
    }
    return result;
  }
};

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
template<Statistics S>
class NormalOperator : public Operator<S>, public AbstractTensor {
 public:
  static constexpr Statistics statistics = S;
  using base_type = Operator<S>;

  // iterate over this using the base_type
  using iterator = typename Operator<S>::iterator;
  using const_iterator = typename Operator<S>::const_iterator;

  using base_type::begin;
  using base_type::end;
  using base_type::cbegin;
  using base_type::cend;
  using base_type::empty;
  using base_type::size;
  using base_type::at;
  using base_type::operator[];

  /// constructs an identity operator
  NormalOperator(Vacuum v = get_default_context().vacuum()) {}

  /// @param creators sequence of creators
  /// @param annihilators sequence of annihilators (in order of particle indices, see the class documentation for more info).
  template <
      typename IndexContainer,
      typename = std::enable_if_t<
          !meta::is_initializer_list_v<std::decay_t<IndexContainer>> &&
          std::is_same_v<typename std::decay_t<IndexContainer>::value_type, Index>>>
  NormalOperator(IndexContainer &&creator_indices,
                 IndexContainer &&annihilator_indices,
                 Vacuum v = get_default_context().vacuum())
      : Operator<S>{}, vacuum_(v), ncreators_(creator_indices.size()) {
    this->reserve(creator_indices.size() + annihilator_indices.size());
    for (const auto &i : creator_indices) {
      this->emplace_back(i, Action::create);
    }
    for (const auto &i : annihilator_indices | ranges::view::reverse) {
      this->emplace_back(i, Action::annihilate);
    }
  }

  /// @param creators sequence of creators
  /// @param annihilators sequence of annihilators (in order of particle indices, see the class documentation for more info).
  NormalOperator(std::initializer_list<Op<S>> creators,
                 std::initializer_list<Op<S>> annihilators,
                 Vacuum v = get_default_context().vacuum())
      : Operator<S>{}, vacuum_(v), ncreators_(size(creators)) {
    for (const auto &op: creators) {
      assert(op.action() == Action::create);
    }
    for (const auto &op: annihilators) {
      assert(op.action() == Action::annihilate);
    }
    this->reserve(size(creators) + size(annihilators));
    this->insert(this->end(), cbegin(creators), cend(creators));
    this->insert(this->end(), crbegin(annihilators), crend(annihilators));
  }

//  /// @param creator_indices sequence of creator indices
//  /// @param annihilator_indices sequence of annihilator indices (in order of particle indices, see the class documentation for more info).
  template <typename I, typename = std::enable_if_t<!std::is_same_v<std::decay_t<I>,Op<S>>>>
  NormalOperator(std::initializer_list<I> creator_indices,
                 std::initializer_list<I> annihilator_indices,
                 Vacuum v = get_default_context().vacuum())
      : Operator<S>{}, vacuum_(v), ncreators_(creator_indices.size()) {
    this->reserve(creator_indices.size() + annihilator_indices.size());
    for (const auto &i: creator_indices) {
      this->emplace_back(i, Action::create);
    }
    for (const auto &i : annihilator_indices | ranges::view::reverse) {
      this->emplace_back(i, Action::annihilate);
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

  /// @return the vacuum state with respect to which the operator is normal-ordered.
  Vacuum vacuum() const { return vacuum_; }
  /// @return the range of creators, in the order of increasing particle index
  auto creators() const { return ranges::view::counted(this->cbegin(), ncreators()); }
  /// @return the range of annihilators, in the order of increasing particle index
  auto annihilators() const { return ranges::view::counted(this->crbegin(), nannihilators()); }
  /// @return the number of creators
  auto ncreators() const { return ncreators_; }
  /// @return the number of annihilators
  auto nannihilators() const { return this->size() - ncreators(); }
  /// @return view of creators and annihilators as a single range
  auto creann() const {
    return ranges::view::concat(creators(), annihilators());
  }

  /// @return number of creators/annihilators
  /// @throw std::logic_error if the operator is not particle number conserving (i.e. if ncreators() != nannihilators() )
  auto rank() const {
    if (ncreators() != nannihilators()) {
      throw std::logic_error("NormalOperator::rank(): ncreators != nannihilators");
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

  NormalOperator &adjoint() {
    static_cast<Operator<S> &>(*this).adjoint();
    ncreators_ = this->size() - ncreators_;
    // TODO rebuild the Hug
    hug_.reset();
    return *this;
  }

  /// @return all possible values returned by label() for this operator type
  static const auto& labels() {
    using namespace std::literals;
    static container::vector<std::wstring> labels_(
        S == Statistics::FermiDirac
            ? std::initializer_list<std::wstring>{L"a"s, L"ã"s}
            : std::initializer_list<std::wstring>{L"b"s, L"ᵬ"s});
    return labels_;
  }

  std::wstring label() const {
    return (S == Statistics::FermiDirac
            ? (vacuum() == Vacuum::Physical ? L"a" : L"ã")
            : (vacuum() == Vacuum::Physical ? L"b" : L"ᵬ"));
  }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    result += (S == Statistics::FermiDirac
               ? (vacuum() == Vacuum::Physical ? L"a" : L"\\tilde{a}")
               : (vacuum() == Vacuum::Physical ? L"b" : L"\\tilde{b}"));
    result += L"^{";
    const auto ncreators = this->ncreators();
    const auto nannihilators = this->nannihilators();
    if (ncreators <
        nannihilators) {  // pad on the left with square underbrackets, i.e. ⎵
      const auto iend = nannihilators - ncreators;
      if (iend > 0) result += L"\\textvisiblespace";
      for (size_t i = 1; i != iend; ++i) {
        result += L"\\,\\textvisiblespace";
      }
      if (ncreators > 0) {
        result += L"\\,";
      }
    }
    for (const auto &o : creators())
      result += o.index().to_latex();
    result += L"}_{";
    if (ncreators >
        nannihilators) {  // pad on the left with square underbrackets, i.e. ⎵
      const auto iend = ncreators - nannihilators;
      if (iend > 0) result += L"\\textvisiblespace";
      for (size_t i = 1; i != iend; ++i) {
        result += L"\\,\\textvisiblespace";
      }
      if (nannihilators > 0) {
        result += L"\\,";
      }
    }
    for (const auto &o : annihilators())
      result += o.index().to_latex();
    result += L"}}";
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

  std::shared_ptr<Expr> clone() const override {
    return std::make_shared<NormalOperator>(*this);
  }

  /// Replaces indices using the index map
  /// @param index_map maps Index to Index
  /// @param tag_transformed_indices if true, will tag transformed indices with
  /// integer 0 and skip any tagged indices that it encounters
  /// @return true if one or more indices changed
  template <template <typename, typename, typename... Args> class Map,
      typename... Args>
  bool transform_indices(const Map<Index, Index, Args...> &index_map,
                         bool tag_transformed_indices = false) {
    bool mutated = false;
    ranges::for_each(*this, [&](auto &&op) {
      if (op.index().transform(index_map, tag_transformed_indices)) mutated = true;
    });
    if (mutated)
      this->reset_hash_value();
    return mutated;
  }

 private:
  Vacuum vacuum_;
  std::size_t ncreators_ = 0;
  using hug_type = HugenholtzVertex<Op<S>, typename Op<S>::TypeEquality>;
  mutable std::unique_ptr<hug_type> hug_;  // only created if needed

  bool static_equal(const Expr &that) const override {
    const auto &that_cast = static_cast<const NormalOperator &>(that);
    if (this->vacuum() == that_cast.vacuum() &&
        this->ncreators() == that_cast.ncreators()) {
      return static_cast<const base_type &>(*this) ==
             static_cast<const base_type &>(*this);
    } else
      return false;
  }

  bool static_less_than(const Expr &that) const override {

    auto range_hash = [](const auto& sized_range) {
      using ranges::begin;
      using ranges::size;
      auto b = begin(sized_range);
      auto e = b + size(sized_range);
      auto val = boost::hash_range(b, e);
      return val;
    };
    auto range_compare = [](const auto& sized_range1, const auto& sized_range2) {
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
        // unlike Tensor comparison, we don't memoize hashes of creators and annihilators separately
        auto cre_hash = range_hash(this->creators());
        auto that_cre_hash = range_hash(that_cast.creators());
        if (cre_hash == that_cre_hash) {
          auto ann_hash = range_hash(this->annihilators());
          auto that_ann_hash = range_hash(that_cast.annihilators());
          if (ann_hash == that_ann_hash)
            return false;
          else {
            return range_compare(this->annihilators(), that_cast.annihilators());
          }
        }
        else {
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
    // same as WickTheorem::can_contract
    auto can_contract = [this](const Op<S> &left, const Op<S> &right) {
      if (is_qpannihilator<S>(left, vacuum_) &&
          is_qpcreator<S>(right, vacuum_)) {
        const auto qpspace_left = qpannihilator_space<S>(left, vacuum_);
        const auto qpspace_right = qpcreator_space<S>(right, vacuum_);
        const auto qpspace_common = intersection(qpspace_left, qpspace_right);
        if (qpspace_common != IndexSpace::null_instance()) return true;
      }
      return false;
    };

    bool result = true;
    /// does not commute with Operator<S>
    /// TODO implement checks of commutativity with Operator<S>
    if (that.is<Operator<S>>()) {
      result = false;
    } else if (that.is<NormalOperator<S>>()) {
      const auto& op_that = that.as<NormalOperator<S>>();
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
           ranges::view::transform(
               [](auto &&op) -> const Index & { return op.index(); });
  }
  AbstractTensor::const_any_view_randsz _ket() const override final {
    return creators() | ranges::view::transform([](auto &&op) -> const Index & {
             return op.index();
           });
  }
  AbstractTensor::const_any_view_rand _braket() const override final {
    return ranges::view::concat(annihilators(), creators()) |
           ranges::view::transform(
               [](auto &&op) -> const Index & { return op.index(); });
  }
  std::size_t _bra_rank() const override final {
    return nannihilators();
  }
  std::size_t _ket_rank() const override final {
    return ncreators();
  }
  Symmetry _symmetry() const override final {
    return S == Statistics::FermiDirac ? Symmetry::antisymm : Symmetry::symm;
  }
  BraKetSymmetry _braket_symmetry() const override final {
    return BraKetSymmetry::nonsymm;
  }
  std::size_t _color() const override final {
    return S == Statistics::FermiDirac ? 1 : 2;
  }
  bool _is_cnumber() const override final {
    return false;
  }
  std::wstring _label() const override final {
    return label();
  }
  std::wstring _to_latex() const override final {
    return to_latex();
  }
  bool _transform_indices(const container::map<Index, Index>& index_map,
                          bool tag_tranformed_indices) override final {
    return transform_indices(index_map, tag_tranformed_indices);
  }
  void _reset_tags() override final {
    ranges::for_each(*this, [](const auto &op) { op.index().reset_tag(); });
  }
  bool operator<(const AbstractTensor& other) const override final {
    auto* other_nop = dynamic_cast<const NormalOperator<S>*>(&other);
    if (other_nop) {
      const Expr* other_expr = static_cast<const Expr*>(other_nop);
      return this->static_less_than(*other_expr);
    }
    else
      return false; // TODO do we compare typeid? labels? probably the latter
  }

  AbstractTensor::any_view_randsz _bra_mutable() override final {
    this->reset_hash_value();
    return ranges::view::counted(this->rbegin(), nannihilators()) |
           ranges::view::transform(
               [](auto &&op) -> Index & { return op.index(); });
  }
  AbstractTensor::any_view_randsz _ket_mutable() override final {
    this->reset_hash_value();
    return ranges::view::counted(this->begin(), ncreators()) |
           ranges::view::transform(
               [](auto &&op) -> Index & { return op.index(); });
  }

};

template<Statistics S>
bool operator==(const NormalOperator<S> &op1, const NormalOperator<S> &op2) {
  return op1.vacuum() == op2.vacuum() && ranges::equal(op1, op2);
}

/// @brief NormalOperatorSequence is a sequence NormalOperator objects, all
/// ordered with respect to same vacuum
///
/// @tparam S specifies the particle statistics
template<Statistics S = Statistics::FermiDirac>
class NormalOperatorSequence : public container::svector<NormalOperator<S>>, public Expr {
 public:
  using base_type = container::svector<NormalOperator<S>>;
  static constexpr Statistics statistics = S;

  using base_type::begin;
  using base_type::end;
  using base_type::cbegin;
  using base_type::cend;
  using base_type::empty;
  using base_type::size;
  using base_type::at;
  using base_type::operator[];

  /// constructs an empty sequence
  NormalOperatorSequence() : vacuum_(get_default_context().vacuum()) {}

  NormalOperatorSequence(std::initializer_list<NormalOperator<S>> operators)
      : base_type(operators) {
    check_vacuum();
  }

  Vacuum vacuum() const { return vacuum_; }

  operator const base_type &() const & { return *this; }
  operator base_type &&() && { return *this; }

  /// @return the total number of Op<S> objects in this
  /// @warning not to be confused with NormalOperatorSequence::size() that
  /// returns the number of NormalOperator<S> objects
  auto opsize() const {
    size_t opsz = 0;
    for (auto &&nop : *this) {
      opsz += nop.size();
    }
    return opsz;
  }

  /// applies (Hermitian) adjoint operation to this
  /// @return reference to @c *this , for daisy-chaining
  NormalOperatorSequence &adjoint() {
    std::reverse(this->begin(), this->end());
    std::for_each(this->begin(), this->end(), [](NormalOperator<S> &op) { op.adjoint(); });
    return *this;
  }

  std::wstring to_latex() const override {
    std::wstring result;
    result = L"{";
    for (const auto &op: *this)
      result += op.to_latex();
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

  bool static_equal(const Expr &that) const override {
    const auto& that_cast = static_cast<const NormalOperatorSequence&>(that);
    if (this->vacuum() == that_cast.vacuum()) {
      if (this->empty()) return true;
      if (this->hash_value() == that.hash_value())
        return static_cast<const base_type&>(*this) == static_cast<const base_type&>(*this);
       else return false;
    } else return false;
  }

};

using BOperator = Operator<Statistics::BoseEinstein>;
using BNOperator = NormalOperator<Statistics::BoseEinstein>;
using BNOperatorSeq = NormalOperatorSequence<Statistics::BoseEinstein>;
using FOperator = Operator<Statistics::FermiDirac>;
using FNOperator = NormalOperator<Statistics::FermiDirac>;
using FNOperatorSeq = NormalOperatorSequence<Statistics::FermiDirac>;

template<Statistics S>
std::wstring to_latex(const NormalOperator<S> &op) {
  return op.to_latex();
}

template<Statistics S>
std::wstring to_latex(const NormalOperatorSequence<S> &opseq) {
  return opseq.to_latex();
}

namespace detail {
  struct OpIdRegistrar {
    OpIdRegistrar();
  };
}  // namespace detail
}  // namespace sequant

#endif // SEQUANT_OP_H
