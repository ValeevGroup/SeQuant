//
// Created by Eduard Valeyev on 2019-03-26.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_HPP
#define SEQUANT_DOMAIN_MBPT_OP_HPP

#include <string>
#include <vector>

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/interval.hpp>

namespace sequant {
namespace mbpt {

/// enumerates the known Operator types
enum class OpType {
  h,       //!< 1-body Hamiltonian
  f,       //!< Fock operator
  g,       //!< 2-body Coulomb
  t,       //!< cluster amplitudes
  lambda,  //!< deexcitation cluster amplitudes
  A,       //!< antisymmetrizer
  L,       //!< left-hand eigenstate
  R,       //!< right-hand eigenstate
  R12,     //!< geminal kernel
  GR,      //!< GR kernel from f12 theory
  C        //!< cabs singles op
};
inline const std::map<OpType, std::wstring> optype2label{
    {OpType::h, L"h"}, {OpType::f, L"f"},      {OpType::g, L"g"},
    {OpType::t, L"t"}, {OpType::lambda, L"λ"}, {OpType::A, L"A"},
    {OpType::L, L"L"}, {OpType::R, L"R"},      {OpType::R12, L"F"}};
inline const std::map<std::wstring, OpType> label2optype =
    ranges::views::zip(ranges::views::values(optype2label),
                       ranges::views::keys(optype2label)) |
    ranges::to<std::map<std::wstring, OpType>>();

/// Operator character relative to Fermi vacuum
enum class OpClass { ex, deex, gen };

/// @return the tensor labels in the cardinal order
std::vector<std::wstring> cardinal_tensor_labels();

std::wstring to_wstring(OpType op);
OpClass to_class(OpType op);

//////////////////////////////

/// mbpt::Operator = Tensor * NormalOperator

/// It is an abstract (index-free) representation of a many-body operator
/// that tracks one or more quantum numbers (e.g., excitation level for
/// single-reference MBPT, particle numbers in each subspace for multi-reference
/// MBPT, etc.) and can be used for operator-level optimizations.
template <typename QuantumNumbers, Statistics S = Statistics::FermiDirac>
class Operator;

template <typename QuantumNumbers>
using FOperator = Operator<QuantumNumbers, Statistics::FermiDirac>;
template <typename QuantumNumbers>
using BOperator = Operator<QuantumNumbers, Statistics::BoseEinstein>;

/// Operator<void> does not track any quantum numbers, this is a base for
/// Operator
template <Statistics S>
class Operator<void, S> : public Expr, public Labeled {
 protected:
  Operator() = default;
  Operator(std::function<std::wstring_view()> label_generator,
           std::function<ExprPtr()> tensor_form_generator)
      : label_generator_(std::move(label_generator)),
        tensor_form_generator_(tensor_form_generator) {}

 public:
  virtual ~Operator() = default;

  std::wstring_view label() const override {
    assert(label_generator_);
    return label_generator_();
  }

  virtual ExprPtr tensor_form() const {
    assert(tensor_form_generator_);
    return tensor_form_generator_();
  }

  /// general n-body operator is not a c-number, so treat all of them as
  /// q-numbers
  bool is_cnumber() const override { return false; }

 private:
  std::function<std::wstring_view()> label_generator_;
  std::function<ExprPtr()> tensor_form_generator_;
};

using FOperatorBase = FOperator<void>;
using BOperatorBase = BOperator<void>;

/// tracks changes in particle number across \c N subspaces

/// \tparam N the number of subspaces; e.g., `N=2` for single-reference
/// theories, `N=3` for multi-reference theories \note In practice, we deal with
/// operators that are change the total # of particles by fixed `k` (`k=0` for
/// particle conserving operators, `k=1` for ionizing, `k=-1` for
/// electron-attaching operators, etc. This instead of tracking the number of
/// particles in each subspace we track the total number of particles and the
/// number of particles in `N-1` subspaces
template <std::size_t N>
class ParticleNumberChange
    : public std::array<boost::numeric::interval<int64_t>, N> {
 public:
  using interval_t = boost::numeric::interval<int64_t>;
  using base_type = std::array<boost::numeric::interval<int64_t>, N>;
  using this_type = ParticleNumberChange<N>;

  ParticleNumberChange() {
    std::fill(this->begin(), this->end(), interval_t{});
  }

  template <typename I, typename = std::enable_if_t<std::is_integral_v<I>>>
  explicit ParticleNumberChange(std::initializer_list<I> i) {
    if (i.size() == N) {
      std::copy(i.begin(), i.end(), this->begin());
    } else {
      if (i.size() == 2)
        this->front() =
            boost::numeric::interval<int64_t>(*(i.begin()), *(i.begin() + 1));
      else
        throw std::runtime_error(
            "ParticleNumberChange<N>(initializer_list i): i.size() must be N; "
            "if N==1 then i.size() can be also 2");
    }
  }

  template <typename I, typename = std::enable_if_t<std::is_integral_v<I>>>
  explicit ParticleNumberChange(
      std::initializer_list<std::initializer_list<I>> i) {
    assert(i.size() == N);
#ifndef NDEBUG
    if (std::find_if(i.begin(), i.end(),
                     [](const auto& ii) { return ii.size() != 2; }) != i.end())
      throw std::runtime_error(
          "ParticleNumberChange<N>(initializer_list<initializer_list> i): each "
          "element of i must contain 2 elements");
#endif
    for (std::size_t c = 0; c != N; ++c) {
      this->operator[](c) = interval_t{*((i.begin() + c)->begin()),
                                       *((i.begin() + c)->begin() + 1)};
    }
  }

  ParticleNumberChange& operator+=(const ParticleNumberChange& other) {
    for (std::size_t c = 0; c != N; ++c) this->operator[](c) += other[c];
    return *this;
  }

  bool operator==(const ParticleNumberChange<N>& b) const {
    return std::equal(
        this->begin(), this->end(), b.begin(),
        [](const auto& ia, const auto& ib) { return equal(ia, ib); });
  }
  bool operator!=(const ParticleNumberChange<N>& b) const {
    return !this->operator==(b);
  }
  template <typename I, std::size_t N_ = N,
            typename = std::enable_if_t<N_ == 1 && std::is_integral_v<I>>>
  bool operator==(I i) const {
    return boost::numeric::interval_lib::compare::possible::operator==(
        this->front(), static_cast<int64_t>(i));
  }
  template <typename I, std::size_t N_ = N,
            typename = std::enable_if_t<N_ == 1 && std::is_integral_v<I>>>
  bool operator!=(I i) const {
    return !this->operator==(i);
  }

  /// @param i an integer
  /// @return true if \p i is in `*this[0]`
  template <typename I, std::size_t N_ = N,
            typename = std::enable_if_t<N_ == 1 && std::is_integral_v<I>>>
  bool in(I i) {
    return boost::numeric::in(static_cast<int64_t>(i), this->front());
  }

  /// @param i an array of N integers
  /// @return true if `i[k]` is in `*this[k]` for all `k`
  template <typename I, typename = std::enable_if_t<std::is_integral_v<I>>>
  bool in(std::array<I, N> i) {
    for (std::size_t c = 0; c != N; ++c) {
      if (!boost::numeric::in(static_cast<int64_t>(i[c]), this->operator[](c)))
        return false;
    }
    return true;
  }

  /// @param i an array of N intervals
  /// @return true if `i[k]` overlaps with `*this[k]` for all `k`
  bool overlaps(std::array<interval_t, N> i) {
    for (std::size_t c = 0; c != N; ++c) {
      if (!boost::numeric::overlap(i[c], this->operator[](c))) return false;
    }
    return true;
  }

  auto hash_value() const {
    static_assert(N > 0);
    auto val = sequant::hash::value(this->operator[](0));
    for (std::size_t c = 1; c != N; ++c) {
      sequant::hash::combine(val, sequant::hash::value(this->operator[](c)));
    }
    return val;
  }
};

template <std::size_t N>
inline bool operator==(const ParticleNumberChange<N>& a,
                       const ParticleNumberChange<N>& b) {
  return a.operator==(b);
}

template <std::size_t N, typename I,
          typename = std::enable_if_t<N == 1 && std::is_integral_v<I>>>
inline bool operator==(const ParticleNumberChange<N>& a, I b) {
  return a.operator==(b);
}

template <std::size_t N>
inline bool equal(const ParticleNumberChange<N>& a,
                  const ParticleNumberChange<N>& b) {
  return operator==(a, b);
}

template <std::size_t N, typename I,
          typename = std::enable_if_t<N == 1 && std::is_integral_v<I>>>
inline bool operator!=(const ParticleNumberChange<N>& a, I b) {
  return a.operator!=(b);
}

template <typename QuantumNumbers, Statistics S>
class Operator : public Operator<void> {
  using this_type = Operator;
  using base_type = Operator<void>;

 protected:
  Operator() = default;

 public:
  Operator(std::function<std::wstring_view()> label_generator,
           std::function<ExprPtr()> tensor_form_generator,
           std::function<void(QuantumNumbers&)> qn_action)
      : base_type(std::move(label_generator), std::move(tensor_form_generator)),
        qn_action_(std::move(qn_action)) {}
  virtual ~Operator() = default;

  QuantumNumbers operator()(const QuantumNumbers& qns) const {
    QuantumNumbers result(qns);
    this->apply_to(result);
    return result;
  }

  virtual QuantumNumbers& apply_to(QuantumNumbers& qns) const {
    assert(qn_action_);
    qn_action_(qns);
    return qns;
  }

  bool static_less_than(const Expr& that) const override {
    assert(that.is<this_type>());
    auto& that_op = that.as<this_type>();

    // compare cardinal tensor labels first, then QN ranks
    auto& cardinal_tensor_labels =
        TensorCanonicalizer::cardinal_tensor_labels();
    const auto this_label = this->label();
    const auto that_label = that_op.label();
    if (this_label == that_label) return this->less_than_rank_of(that_op);
    const auto this_is_cardinal_it = ranges::find_if(
        cardinal_tensor_labels,
        [&this_label](const std::wstring& l) { return l == this_label; });
    const auto this_is_cardinal =
        this_is_cardinal_it != ranges::end(cardinal_tensor_labels);
    const auto that_is_cardinal_it = ranges::find_if(
        cardinal_tensor_labels,
        [&that_label](const std::wstring& l) { return l == that_label; });
    const auto that_is_cardinal =
        that_is_cardinal_it != ranges::end(cardinal_tensor_labels);
    if (this_is_cardinal && that_is_cardinal) {
      if (this_is_cardinal_it != that_is_cardinal_it)
        return this_is_cardinal_it < that_is_cardinal_it;
      else
        return this->less_than_rank_of(that_op);
    } else if (this_is_cardinal && !that_is_cardinal)
      return true;
    else if (!this_is_cardinal && that_is_cardinal)
      return false;
    else {  // !this_is_cardinal && !that_is_cardinal
      if (this_label == that_label)
        return this->less_than_rank_of(that_op);
      else
        return this_label < that_label;
    }
  }

  bool commutes_with_atom(const Expr& that) const override {
    assert(that.is_cnumber() || that.is<this_type>());
    if (that.is_cnumber())
      return true;
    else {
      auto& that_op = that.as<this_type>();
      // if this has annihilators/creators in same space as that has
      // creator/annihilators return false
      const auto delta_this = (*this)(QuantumNumbers{});
      const auto delta_that = (that_op)(QuantumNumbers{});
      auto gtzero = [](const auto& interval) {
        using interval_type = std::decay_t<decltype(interval)>;
        using base_type = typename interval_type::base_type;
        return boost::numeric::overlap(
            interval, interval_type{1, std::numeric_limits<base_type>::max()});
      };
      auto ltzero = [](const auto& interval) {
        using interval_type = std::decay_t<decltype(interval)>;
        using base_type = typename interval_type::base_type;
        return boost::numeric::overlap(
            interval, interval_type{std::numeric_limits<base_type>::min(), -1});
      };
      bool this_has_cre_0 = gtzero(delta_this[1]);
      bool this_has_ann_0 = ltzero(delta_this[1]);
      bool that_has_cre_0 = gtzero(delta_that[1]);
      bool that_has_ann_0 = ltzero(delta_that[1]);
      bool this_has_cre_1 = gtzero(delta_this[0] - delta_this[1]);
      bool this_has_ann_1 = ltzero(delta_this[0] - delta_this[1]);
      bool that_has_cre_1 = gtzero(delta_that[0] - delta_that[1]);
      bool that_has_ann_1 = ltzero(delta_that[0] - delta_that[1]);
      auto result = !((this_has_cre_0 && that_has_ann_0) ||
                      (that_has_cre_0 && this_has_ann_0) ||
                      (this_has_cre_1 && that_has_ann_1) ||
                      (that_has_cre_1 && this_has_ann_1));
      if (result == true)
        std::wcout << this->to_latex() << " commutes with " << that.to_latex()
                   << std::endl;
      else
        std::wcout << this->to_latex() << " does NOT commute with "
                   << that.to_latex() << std::endl;
      return result;
    }
  }

  void adjoint() override {
    if constexpr (std::is_same_v<QuantumNumbers,
                                 mbpt::ParticleNumberChange<2>>) {
      const auto dN = (*this)(QuantumNumbers{});
      using qns_t = std::decay_t<decltype(dN)>;

      const auto lbl = this->label();
      const auto tnsr = this->tensor_form();
      *this = Operator{
          [=]() -> std::wstring_view { return lbl; },  // return same label
          [=]() -> ExprPtr {
            return sequant::adjoint(tnsr);  // return adjoint of tensor form
          },
          [=](qns_t& qn) {
            qn += qns_t{{dN[0].upper() * -1, dN[0].lower() * -1},
                        {dN[1].upper() * -1, dN[1].lower() * -1}};
            return qn;  // return modified qns
          }};
    } else
      throw std::runtime_error(
          "mbpt::Operator::adjoint: only implemented for the single-reference "
          "case");
  }

 private:
  std::function<void(QuantumNumbers&)> qn_action_;

  bool less_than_rank_of(const this_type& that) const {
    return (*this)(QuantumNumbers{}) < that(QuantumNumbers{});
  }

  Expr::type_id_type type_id() const override {
    return get_type_id<this_type>();
  };

  ExprPtr clone() const override { return ex<this_type>(*this); }

  std::wstring to_latex() const override {
    const auto dN = (*this)(QuantumNumbers{});
    if constexpr (std::is_same_v<QuantumNumbers,
                                 mbpt::ParticleNumberChange<2>>) {
      assert(dN[0].lower() == dN[0].upper());
      const auto dN_total = dN[0].lower();

      auto lbl = std::wstring(this->label());
      std::wstring result = L"{\\hat{" + lbl + L"}";
      auto it = label2optype.find(lbl);
      if (it != label2optype.end()) {  // handle special cases
        const auto optype = it->second;
        if (optype == OpType::lambda) {  // λ -> \lambda
          result = L"{\\hat{\\lambda}";
        }
        if (to_class(optype) == OpClass::gen) {
          result += L"}";
          return result;
        }
      }
      // generic operator
      const auto dN_second = dN[1];
      const auto dN_definite = dN_second.lower() == dN_second.upper();
      decltype(dN_second) dN_first(dN_total - dN_second.upper(),
                                   dN_total - dN_second.lower());
      if (dN_definite) {
        if (dN_total == 0) {  // N-conserving
          using std::abs;
          result += L"_{" + std::to_wstring(abs(dN_first.lower())) + L"}";
        } else {  // N-nonconserving
          result += L"_{" + std::to_wstring(dN_first.lower()) + L"}^{" +
                    std::to_wstring(dN_second.lower()) + L"}";
        }
      } else {
        result += L"_{[" + std::to_wstring(dN_first.lower()) + L"," +
                  std::to_wstring(dN_first.upper()) + L"]}^{[" +
                  std::to_wstring(dN_second.lower()) + L"," +
                  std::to_wstring(dN_second.upper()) + L"]}";
      }
      result += L"}";
      return result;
    } else
      throw std::runtime_error(
          "mbpt::Operator::to_latex: only implemented for the single-reference "
          "case");
  }

  hash_type memoizing_hash() const override {
    using std::begin;
    using std::end;
    auto qns = (*this)(QuantumNumbers{});
    auto val = sequant::hash::value(qns);
    sequant::hash::combine(val, std::wstring(this->label()));
    hash_value_ = val;
    return *hash_value_;
  }

};  // Operator

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_OP_HPP
