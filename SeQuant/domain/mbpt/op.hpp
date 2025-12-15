//
// Created by Eduard Valeyev on 2019-03-26.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_HPP
#define SEQUANT_DOMAIN_MBPT_OP_HPP

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/interval.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/strong.hpp>

#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/range/primitives.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <range/v3/view/zip.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace sequant {
namespace mbpt {

DEFINE_STRONG_TYPE_FOR_INTEGER(nₚ, std::int64_t);  // define nₚ
DEFINE_STRONG_TYPE_FOR_INTEGER(nₕ, std::int64_t);  // define nₕ

#ifndef DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT
#define DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(OP)                      \
  inline ExprPtr OP(std::int64_t Rank) { return OP(nₚ(Rank), nₕ(Rank)); } \
  inline ExprPtr OP(nₚ Rank) { return OP(Rank, nₕ(Rank)); }
#endif  // DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT

template <typename QuantumNumbers>
bool is_vacuum(QuantumNumbers qns);

/// converts an IndexSpace::Type to IndexSpace using default quantum number set
inline IndexSpace make_space(const IndexSpace::Type& type) {
  return get_default_context().index_space_registry()->retrieve(type,
                                                                Spin::any);
}

/// enumerates the known Operator types
enum class OpType {
  h,    //!< 1-body Hamiltonian
  f,    //!< Fock operator
  f̃,    //!< closed Fock operator (i.e. Fock operator due to fully-occupied
        //!< orbitals)
  s,    //!< 1-body overlap
  θ,    //!< general fock space operator
  g,    //!< 2-body Coulomb
  t,    //!< cluster amplitudes
  λ,    //!< deexcitation cluster amplitudes
  A,    //!< antisymmetrizer
  S,    //!< particle symmetrizer
  L,    //!< left-hand eigenstate
  R,    //!< right-hand eigenstate
  R12,  //!< geminal kernel
  GR,   //!< GR kernel from f12 theory
  C,    //!< cabs singles op
  RDM,  //!< RDM
  RDMCumulant,  //!< RDM cumulant
  δ,            //!< Kronecker delta (=identity) operator
  h_1,          //!< Hamiltonian perturbation
  t_1,          //!< first order perturbed excitation cluster operator
  λ_1,          //!< first order perturbed deexcitation cluster operator
};

/// maps operator types to their labels
inline const container::map<OpType, std::wstring> optype2label{
    {OpType::h, L"h"},
    {OpType::f, L"f"},
    {OpType::f̃, L"f̃"},
    {OpType::s, overlap_label()},
    {OpType::δ, kronecker_label()},
    {OpType::g, L"g"},
    {OpType::θ, L"θ"},
    {OpType::t, L"t"},
    {OpType::λ, L"λ"},
    {OpType::A, L"A"},
    {OpType::S, L"S"},
    {OpType::L, L"L"},
    {OpType::R, L"R"},
    {OpType::R12, L"F"},
    {OpType::GR, L"GR"},
    {OpType::C, L"C"},
    {OpType::RDM, L"γ"},
    // see https://en.wikipedia.org/wiki/Cumulant
    {OpType::RDMCumulant, L"κ"},
    {OpType::h_1, L"h¹"},
    {OpType::t_1, L"t¹"},
    {OpType::λ_1, L"λ¹"}};

/// maps operator labels to their types
inline const container::map<std::wstring, OpType> label2optype =
    ranges::views::zip(ranges::views::values(optype2label),
                       ranges::views::keys(optype2label)) |
    ranges::to<container::map<std::wstring, OpType>>();

/// Operator character relative to Fermi vacuum
enum class OpClass { ex, deex, gen };

/// @brief Struct which holds parameters for operator construction
struct OpParams {
  std::size_t order = 1;  ///< perturbation order (for _pt operators)
  std::optional<size_t> nbatch = std::nullopt;  ///< number of batching indices
  container::svector<std::size_t> batch_ordinals{};
  ///< custom batching index ordinals (empty = no batching)
  bool skip1 = false;  ///< skip single excitations (for sum operators)

  /// @brief Validates the parameters for consistency and correctness
  void validate() const {
    SEQUANT_ASSERT(!(nbatch && !batch_ordinals.empty()) &&
                   "OpParams: Cannot specify both nbatch and batch_ordinals");
    // ensure batch ordinals are unique
    if (!batch_ordinals.empty()) {
#ifdef SEQUANT_ASSERT_ENABLED
      SEQUANT_ASSERT(ranges::is_sorted(batch_ordinals) &&
                     "OpParams: batch_ordinals must be sorted");
      auto duplicate = ranges::adjacent_find(batch_ordinals);
      SEQUANT_ASSERT(duplicate == batch_ordinals.end() &&
                     "OpParams: batch_ordinals must contain unique values");
#endif
    }
  }
};

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

  /// @brief Constructs an Operator with a label generator and a tensor form
  /// generator
  /// @param label_generator A function that generates a label for the operator
  /// @param tensor_form_generator A function that generates the tensor form of
  /// the operator
  Operator(std::function<std::wstring_view()> label_generator,
           std::function<ExprPtr()> tensor_form_generator)
      : label_generator_(std::move(label_generator)),
        tensor_form_generator_(tensor_form_generator) {}

 public:
  virtual ~Operator() = default;

  /// @return label of the operator
  std::wstring_view label() const override {
    SEQUANT_ASSERT(label_generator_);
    return label_generator_();
  }

  /// @return tensor form of the operator
  virtual ExprPtr tensor_form() const {
    SEQUANT_ASSERT(tensor_form_generator_);
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

struct default_qns_tag {};

// clang-format off
/// Tracks changes in \c N quantum numbers

/// Represents changes in a set of quantum numbers; this is useful for
/// tracking the change of quantum numbers induced by a many-body operator, such as the number of particles,
/// the number of quasiparticles, the number of ops (creators/annihilators) in each subspace, etc.
/// For example, to operator products expressed in normal order with respect to physical vacuum it is sufficient to track
/// the number of creators and annihilators; For the fermi vacuum case, the number of creators and annihilators in each
/// subspace becomes important, hence the number of ops is tracked for each base space (determined by the IndexSpaceRegistry object in Context).
/// The interval representation is necessary to dictate how many creators or annihilators could be in each subspace.
/// This is pertinent when user defined hole_space or particle_space are NOT base spaces.
/// Since the choice of space partitioning is up to the user, the base class must be a dynamic container.
/// \tparam Tag a tag type to distinguish different instances of QuantumNumberChange<N>
/// \tparam QNV the quantum number value type, defaults to \c std::int64_t
// clang-format on
template <typename QNV = std::int64_t, typename Tag = default_qns_tag>
class QuantumNumberChange
    : public container::svector<
          boost::numeric::interval<std::make_signed_t<QNV>>, 8> {
 public:
  using QNC = std::make_signed_t<QNV>;  // change in quantum numbers
  using interval_t = boost::numeric::interval<QNC>;
  using base_type =
      container::svector<boost::numeric::interval<std::make_signed_t<QNV>>, 8>;
  using this_type = QuantumNumberChange<QNV, Tag>;

  std::size_t size() const {
    if (get_default_context().vacuum() == Vacuum::Physical) {
      return 2;
    } else if (get_default_context().vacuum() == Vacuum::SingleProduct) {
      auto isr = get_default_context().index_space_registry();
      const auto& isr_base_spaces = isr->base_spaces();
      SEQUANT_ASSERT(isr_base_spaces.size() > 0);
      return isr_base_spaces.size() * 2;
    } else {
      throw std::logic_error("unknown Vacuum type");
    }
  }

  /// initializes all values with zeroes
  QuantumNumberChange() {
    this->resize(this->size());
    SEQUANT_ASSERT(this->base().size() != 0);
    std::fill(this->begin(), this->end(), interval_t{});
  }

  /// constructs QuantumNumberChange from a sequence of elements convertible to
  /// QNV
  /// \tparam I the type of the initializer_list elements
  /// \param i the sequence of QNV-convertible elements
  template <typename I, typename Range,
            typename = std::enable_if_t<
                meta::is_range_v<std::remove_reference_t<Range>> &&
                std::is_convertible_v<I, interval_t>>>
  explicit QuantumNumberChange(Range&& i) : QuantumNumberChange() {
    SEQUANT_ASSERT(i.size() == size());
    std::copy(i.begin(), i.end(), this->begin());
  }

  /// constructs QuantumNumberChange from a sequence of elements convertible to
  /// QNV
  /// \tparam I the type of the initializer_list elements
  /// \param i the sequence of QNV-convertible elements
  template <typename I, typename = std::enable_if_t<std::is_convertible_v<
                            std::initializer_list<I>, interval_t>>>
  explicit QuantumNumberChange(
      std::initializer_list<std::initializer_list<I>> i)
      : QuantumNumberChange() {
    SEQUANT_ASSERT(i.size() == size());
#ifdef SEQUANT_ASSERT_ENABLED
    SEQUANT_ASSERT(
        std::find_if(i.begin(), i.end(),
                     [](const auto& ii) { return ii.size() != 2; }) ==
            i.end() &&
        "QuantumNumberChange<N>(initializer_list<initializer_list> i): each "
        "element of i must contain 2 elements");
#endif
    for (std::size_t c = 0; c != size(); ++c) {
      this->operator[](c) = interval_t{*((i.begin() + c)->begin()),
                                       *((i.begin() + c)->begin() + 1)};
    }
  }

  QuantumNumberChange& operator+=(const QuantumNumberChange& other) {
    for (std::size_t c = 0; c != size(); ++c) this->operator[](c) += other[c];
    return *this;
  }

  bool operator==(const this_type& b) const {
    return std::equal(
        this->begin(), this->end(), b.begin(),
        [](const auto& ia, const auto& ib) { return equal(ia, ib); });
  }
  bool operator!=(const this_type& b) const { return !this->operator==(b); }

  bool operator<(const this_type& that) const {
    return ranges::lexicographical_compare(
        *this, that, [](const interval_t& a, const interval_t& b) {
          if (a.lower() != b.lower()) return a.lower() < b.lower();
          return a.upper() < b.upper();
        });
  }

  // determines the number of physical vacuum creators and annihilators for the
  // active particle and hole space from the Context. for general operators this
  // is not defined. for example: O_{e_1}^{i_1 m_1} a_{i_1 m_1}^{e_1} asking for
  // the active particle annihilators in this example is nonsense and will
  // return -1.

  /// @brief determines the number of creators in the particle space defined in
  /// the current context
  interval_t ncre_particles() {
    const auto& qnvec = this->base();
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    interval_t result = 0;
    for (unsigned int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      const auto intersect_type =
          base_space.attr()
              .intersection(isr->particle_space(base_space.qns()).attr())
              .type();
      if (IndexSpace::Type{} != intersect_type) {
        result += qnvec[2 * i];
      }
    }
    return result;
  }

  /// @brief determines the number of annihilators in the particle space defined
  /// in the current context
  interval_t nann_particles() {
    const auto& qnvec = this->base();
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    interval_t result = 0;
    for (unsigned int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      const auto intersect_type =
          base_space.attr()
              .intersection(isr->particle_space(base_space.qns()).attr())
              .type();
      if (IndexSpace::Type{} != intersect_type) {
        result += qnvec[2 * i + 1];
      }
    }
    return result;
  }

  /// @brief determines the number of creators in the hole space defined in the
  /// current context
  interval_t ncre_holes() {
    const auto& qnvec = this->base();
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    interval_t result = 0;
    for (unsigned int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      const auto intersect_type =
          base_space.attr()
              .intersection(isr->hole_space(base_space.qns()).attr())
              .type();
      if (IndexSpace::Type{} != intersect_type) {
        result += qnvec[2 * i];
      }
    }
    return result;
  }

  /// @brief determines the number of annihilators in the hole space defined in
  /// the current context
  interval_t nann_holes() {
    const auto& qnvec = this->base();
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    interval_t result = 0;
    for (unsigned int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      const auto intersect_type =
          base_space.attr()
              .intersection(isr->hole_space(base_space.qns()).attr())
              .type();
      if (IndexSpace::Type{} != intersect_type) {
        result += qnvec[2 * i + 1];
      }
    }
    return result;
  }

  /// tests whether particular changes in quantum number change
  /// @param i an integer
  /// @return true if \p i is in `*this[0]`
  template <typename I>
  bool in(I i) {
    return boost::numeric::in(static_cast<int64_t>(i), this->front());
  }

  ///
  /// @param i an integer
  /// @return true if \p i is in `*this[0]`
  template <typename I, typename = std::enable_if_t<std::is_integral_v<I>>>
  bool in(std::initializer_list<I> i) {
    SEQUANT_ASSERT(i.size() == size());
    std::array<I, 4> i_arr;
    std::copy(i.begin(), i.end(), i_arr.begin());
    return this->in(i_arr);
  }

  /// @param i an array of N intervals
  /// @return true if `i[k]` overlaps with `*this[k]` for all `k`
  bool overlaps_with(base_type i) {
    for (std::size_t c = 0; c != this->size(); ++c) {
      if (!boost::numeric::overlap(i[c], this->operator[](c))) {
        return false;
      }
    }
    return true;
  }

  auto hash_value() const {
    SEQUANT_ASSERT(size() > 0);
    auto val = sequant::hash::value(this->operator[](0));
    for (std::size_t c = 1; c != size(); ++c) {
      sequant::hash::combine(val, sequant::hash::value(this->operator[](c)));
    }
    return val;
  }

 private:
  auto& base() { return static_cast<base_type&>(*this); }
};

template <std::size_t N, typename Tag, typename QNV>
inline bool operator==(const QuantumNumberChange<Tag, QNV>& a,
                       const QuantumNumberChange<Tag, QNV>& b) {
  return a.operator==(b);
}

template <std::size_t N, typename Tag, typename QNV>
inline bool operator!=(const QuantumNumberChange<Tag, QNV>& a,
                       const QuantumNumberChange<Tag, QNV>& b) {
  return !(a == b);
}

template <std::size_t N, typename Tag, typename QNV, typename I,
          typename = std::enable_if_t<N == 1 && std::is_integral_v<I>>>
inline bool operator==(const QuantumNumberChange<Tag, QNV>& a, I b) {
  return a.operator==(b);
}

template <std::size_t N, typename Tag, typename QNV>
inline bool equal(const QuantumNumberChange<Tag, QNV>& a,
                  const QuantumNumberChange<Tag, QNV>& b) {
  return operator==(a, b);
}

template <std::size_t N, typename Tag, typename QNV, typename I,
          typename = std::enable_if_t<N == 1 && std::is_integral_v<I>>>
inline bool operator!=(const QuantumNumberChange<Tag, QNV>& a, I b) {
  return a.operator!=(b);
}

// clang-format off
/// algebra of operators normal order with respect to physical vacuum
/// can be screened by tracking the number of creators and annihilators.
/// the order of of elements is {# of creators, # of annihilators}
/// \note use signed integer, although could use unsigned in this case,
/// so that can represent quantum numbers and their changes by the same type
// clang-format on
using qns_t = mbpt::QuantumNumberChange<>;
using qninterval_t = typename qns_t::interval_t;
/// changes in quantum number represented by quantum numbers themselves
using qnc_t = qns_t;
using op_t = mbpt::Operator<qnc_t>;

/// combines 2 sets of quantum numbers using Wick's theorem
qns_t combine(qns_t, qns_t);

/// @brief Constructs quantum numbers for an excitation operator based on the
/// defined context
/// @param k the rank of the operator
/// @param SQN the spin quantum number
qns_t excitation_type_qns(std::size_t k,
                          IndexSpace::QuantumNumbers SQN = Spin::any);

/// @brief Constructs quantum numbers for an excitation operator based on the
/// defined context. Sometimes we want to guarantee that a qns has an interval
/// from 0 to \p k regardless of base spaces
/// @param k the rank of the operator, QN has the interval from 0 to \p k
/// @param SQN the spin quantum number
qns_t interval_excitation_type_qns(std::size_t k,
                                   IndexSpace::QuantumNumbers SQN = Spin::any);

/// @brief Constructs quantum numbers for an deexcitation operator based on the
/// defined context
/// @param k the rank of the operator
/// @param SQN the spin quantum number
qns_t deexcitation_type_qns(std::size_t k,
                            IndexSpace::QuantumNumbers SQN = Spin::any);

/// @brief Constructs quantum numbers for an deexcitation operator based on the
/// defined context. Sometimes we want to guarantee that a qns has an interval
/// from 0 to \p k regardless of base spaces
/// @param k the rank of the operator, QN has the interval from 0 to \p k
/// @param SQN the spin quantum number
qns_t interval_deexcitation_type_qns(
    std::size_t k, IndexSpace::QuantumNumbers SQN = Spin::any);

/// @brief Constructs quantum numbers for a general operator based on the
/// defined context
/// @param k the rank of the operator
qns_t general_type_qns(std::size_t k);

/// @brief Constructs quantum numbers for a generic excitation operator
/// @param particle_rank number of operators in the particle space
/// @param hole_rank number operators in the hole space
/// @param particle_space particle space within the defined context
/// @param hole_space hole space within the defined context
qns_t generic_excitation_qns(std::size_t particle_rank, std::size_t hole_rank,
                             IndexSpace particle_space, IndexSpace hole_space,
                             IndexSpace::QuantumNumbers SQN = Spin::any);

/// @brief Constructs quantum numbers for a generic deexcitation operator
/// @param particle_rank number of operators in the particle space
/// @param hole_rank number operators in the hole space
/// @param particle_space particle space within the defined context
/// @param hole_space hole space within the defined context
qns_t generic_deexcitation_qns(std::size_t particle_rank, std::size_t hole_rank,
                               IndexSpace particle_space, IndexSpace hole_space,
                               IndexSpace::QuantumNumbers SQN = Spin::any);

inline namespace op {
namespace tensor {
namespace detail {
ExprPtr expectation_value_impl(ExprPtr expr,
                               std::vector<std::pair<int, int>> nop_connections,
                               bool use_top, bool full_contractions);
}  // namespace detail

/// @brief computes the reference expectation value of a tensor-level expression
/// @param expr input expression
/// @param nop_connections connectivity information
/// @param use_top if true, WickTheorem uses topological equivalence of terms
ExprPtr ref_av(ExprPtr expr,
               std::vector<std::pair<int, int>> nop_connections = {},
               bool use_top = true);

/// @brief computes the vacuum expectation value of a tensor-level expression,
/// forces full contractions in WickTheorem
/// @param expr input expression
/// @param nop_connections connectivity information
/// @param use_top if true, WickTheorem uses topological equivalence of terms
ExprPtr vac_av(ExprPtr expr,
               std::vector<std::pair<int, int>> nop_connections = {},
               bool use_top = true);
}  // namespace tensor
}  // namespace op

}  // namespace mbpt

/// @param qns the quantum numbers to adjoint
/// @return the adjoint of \p qns
mbpt::qns_t adjoint(mbpt::qns_t qns);

namespace mbpt {

// clang-format off
/// @brief makes a tensor-level many-body operator

/// A many-body operator has the following generic form:
/// \f$ \frac{1}{P} T_{b_1 b_2 \dots b_B}^{k_1 k_2 \dots k_K} A^{b_1 b_2 \dots b_B}_{k_1 k_2 \dots k_K} \f$
/// where \f$ \{B,K\} \f$ are number of bra/ket indices of \f$ T \f$ or, equivalently, the number of creators/annihilators
/// of normal-ordered (w.r.t. the default vacuum) operator \f$ A \f$.
/// Hence \f$ \{ b_i \} \f$ / \f$ \{ k_i \} \f$ are (quasi)particle creation/annihilation indices.
/// For example, for fermionic operators relative to Fermi vacuum these are:
/// - (pure) _excitation_: (active) unoccupied/occupied, respectively;
/// - (pure) _deexcitation_: occupied/unoccupied, respectively;
/// For _generic_ operators (neither excitation nor deexcitation) complete basis indices are assumed by default.
///
/// \f$ P \f$ is the "normalization" factor and depends on the vacuum used to define \f$ A \f$,
/// and indices \f$ \{ b_i \} \f$ / \f$ \{ k_i \} \f$.
/// @note The choice of computational basis can be controlled by the default Context.
///       See `SeQuant/core/context.hpp` and `SeQuant/mbpt/context.hpp`
/// @warning Tensor \f$ T \f$ will be antisymmetrized if `get_default_context().spbasis() == SPBasis::Spinor`, else it will be particle-symmetric; the latter is only valid if # of bra and ket indices coincide.
/// @internal bless the maker and his water
// clang-format on
template <Statistics S>
class OpMaker {
  using IndexContainer = container::svector<Index>;
  using IndexSpaceContainer = container::svector<IndexSpace>;

 public:
  /// @param[in] op the operator type
  /// @param[in] cre_list list of creator indices
  /// @param[in] ann_list list of annihilator indices
  template <typename IndexSpaceTypeRange1, typename IndexSpaceTypeRange2>
  OpMaker(OpType op, const cre<IndexSpaceTypeRange1>& cre_list,
          const ann<IndexSpaceTypeRange2>& ann_list)
      : op_(op),
        cre_spaces_(cre_list.begin(), cre_list.end()),
        ann_spaces_(ann_list.begin(), ann_list.end()) {
    SEQUANT_ASSERT(ncreators() > 0 || nannihilators() > 0);
  }

  /// @param[in] op the operator type
  /// @param[in] nc number of bra indices/creators
  /// @param[in] na number of ket indices/annihilators
  OpMaker(OpType op, ncre nc, nann na);

  /// @brief creates a particle-conserving replacement operator
  /// @param[in] op the operator type
  /// @param[in] rank particle rank of the operator (# of creators = # of
  /// annihilators = @p rank )
  OpMaker(OpType op, std::size_t rank);

  /// @param[in] op the operator type
  /// @param[in] nc number of bra indices/creators
  /// @param[in] na number of ket indices/annihilators
  /// @param[in] cre_space IndexSpace referred to be the creator
  /// @param[in] ann_space IndexSpace referred to be the annihilators
  OpMaker(OpType op, ncre nc, nann na, const cre<IndexSpace>& cre_space,
          const ann<IndexSpace>& ann_space);

  /// @brief Creates operator from OpParams
  /// @param[in] op the operator type
  /// @param[in] nc number of bra indices/creators
  /// @param[in] na number of ket indices/annihilators
  /// @param[in] params named parameters for operator construction
  OpMaker(OpType op, ncre nc, nann na, const OpParams& params);

  enum class UseDepIdx {
    /// bra/cre indices depend on ket
    Bra,
    /// ket/ann indices depend on bra
    Ket,
    /// use plain indices
    None
  };

  /// struct to hold the information about the operator
  struct OpInfo {
    container::svector<Index> creidxs;  //!< creator indices
    container::svector<Index> annidxs;  //!< annihilator indices
    sequant::intmax_t mult;             //!< normalization factor
    Symmetry opsymm;                    //!< symmetry of the operator
    UseDepIdx dep;                      //!< dependency of the bra/ket indices
  };

  // clang-format off
  /// @param[in] dep_opt if given, controls whether bra (`*dep_opt == UseDepIdx::Bra`)
  /// / ket (`*dep_opt == UseDepIdx::Ket`) indices
  /// are dependent on the, respectively, ket/bra indices
  /// (i.e., use them as protoindices);
  /// if (`*dep_opt == UseDepIdx::None`) then plain indices are used; if
  /// \p dep_opt is not given then the default is determined by the MBPT context.
  /// @param[in] opsymm_opt if given, controls whether (anti)symmetric
  /// tensor is returned; if \p opsymm_opt is not given then the default is
  /// determined by the MBPT context.
  // clang-format on
  ExprPtr operator()(std::optional<UseDepIdx> dep_opt = {},
                     std::optional<Symmetry> opsymm_opt = {}) const;

  /// @brief Creates an OpInfo struct containing creator and annihilator
  /// indices, normalization factor, symmetry, and dependency information.
  /// @param cre_spaces A container of IndexSpace objects representing the
  /// creator indices
  /// @param ann_spaces A container of IndexSpace objects representing the
  /// annihilator indices
  /// @param dep An optional parameter specifying the dependency of indices.
  /// @return An OpInfo struct containing the created indices, normalization
  /// factor, symmetry, and dependency information.
  static OpInfo build_op_info(const IndexSpaceContainer& cre_spaces,
                              const IndexSpaceContainer& ann_spaces,
                              UseDepIdx dep = UseDepIdx::None) {
    const bool symm = get_default_context().spbasis() ==
                      SPBasis::Spinor;  // antisymmetrize if spin-orbital basis
    const auto dep_bra = dep == UseDepIdx::Bra;
    const auto dep_ket = dep == UseDepIdx::Ket;

    // not sure what it means to use nonsymmetric operator if nbra != nket
    if (!symm)
      SEQUANT_ASSERT(ranges::size(cre_spaces) == ranges::size(ann_spaces));

    auto make_idx_vector = [](const auto& spaces) {
      return spaces | ranges::views::transform([](const IndexSpace& space) {
               return Index::make_tmp_index(space);
             }) |
             ranges::to<container::svector<Index>>();
    };

    auto make_depidx_vector = [](const auto& spaces, auto&& protoidxs) {
      return spaces |
             ranges::views::transform([&protoidxs](const IndexSpace& space) {
               return Index::make_tmp_index(space, protoidxs, true);
             }) |
             ranges::to<container::svector<Index>>();
    };

    container::svector<Index> creidxs, annidxs;
    if (dep_bra) {
      annidxs = make_idx_vector(ann_spaces);
      creidxs = make_depidx_vector(cre_spaces, annidxs);
    } else if (dep_ket) {
      creidxs = make_idx_vector(cre_spaces);
      annidxs = make_depidx_vector(ann_spaces, creidxs);
    } else {
      creidxs = make_idx_vector(cre_spaces);
      annidxs = make_idx_vector(ann_spaces);
    }

    using ranges::size;
    const auto mult =
        symm ? factorial(size(cre_spaces)) * factorial(size(ann_spaces))
             : factorial(size(cre_spaces));
    const auto opsymm = symm ? (S == Statistics::FermiDirac ? Symmetry::Antisymm
                                                            : Symmetry::Symm)
                             : Symmetry::Nonsymm;

    return OpInfo{creidxs, annidxs, mult, opsymm, dep};
  }

  /// @tparam TensorGenerator callable with signature
  /// `TensorGenerator(range<Index>, range<Index>, Symmetry)` that returns a
  /// Tensor with the respective bra/cre and ket/ann indices and of the given
  /// symmetry
  /// @param[in] cre_spaces creator IndexSpaces
  /// @param[in] ann_spaces annihilator IndexSpaces
  /// @param[in] tensor_generator the callable that generates the tensor
  /// @param[in] dep whether to use dependent indices
  template <typename TensorGenerator>
  static ExprPtr make(const IndexSpaceContainer& cre_spaces,
                      const IndexSpaceContainer& ann_spaces,
                      TensorGenerator&& tensor_generator,
                      UseDepIdx dep = UseDepIdx::None) {
    const auto op_info = build_op_info(cre_spaces, ann_spaces, dep);

    const auto t =
        tensor_generator(op_info.creidxs, op_info.annidxs, op_info.opsymm);
    return ex<Constant>(rational{1, op_info.mult}) * t *
           ex<NormalOperator<S>>(cre(op_info.creidxs), ann(op_info.annidxs),
                                 get_default_context().vacuum());
  }

  /// @tparam TensorGenerator callable with signature
  /// `TensorGenerator(range<Index>, range<Index>, Symmetry)` that returns a
  /// Tensor with the respective bra/cre and ket/ann indices and of the given
  /// symmetry
  /// @param[in] cre_spaces creator IndexSpaces as an initializer list
  /// @param[in] ann_spaces annihilator IndexSpaces as an initializer list
  /// @param[in] tensor_generator the callable that generates the tensor
  /// @param[in] csv whether to use dependent indices
  template <typename TensorGenerator>
  static ExprPtr make(std::initializer_list<IndexSpace::Type> cre_spaces,
                      std::initializer_list<IndexSpace::Type> ann_spaces,
                      TensorGenerator&& tensor_generator,
                      UseDepIdx csv = UseDepIdx::None) {
    IndexSpaceContainer cre_vec(cre_spaces.begin(), cre_spaces.end());
    IndexSpaceContainer ann_vec(ann_spaces.begin(), ann_spaces.end());
    return OpMaker::make(cre_vec, ann_vec,
                         std::forward<TensorGenerator>(tensor_generator), csv);
  }

  /// @tparam TensorGenerator callable with signature
  /// `TensorGenerator(range<Index>, range<Index>, range<Index>, Symmetry)` that
  /// returns a Tensor with the respective bra/cre, ket/ann, and batch indices
  /// and of the given symmetry
  /// @param[in] cre_spaces creator IndexSpaces
  /// @param[in] ann_spaces annihilator IndexSpaces
  /// @param[in] batch_indices batch indices
  /// @param[in] tensor_generator the callable that generates the tensor
  /// @param[in] dep whether to use dependent indices
  template <typename TensorGenerator>
  static ExprPtr make(const IndexSpaceContainer& cre_spaces,
                      const IndexSpaceContainer& ann_spaces,
                      const IndexContainer& batch_indices,
                      TensorGenerator&& tensor_generator,
                      UseDepIdx dep = UseDepIdx::None) {
    mbpt::check_for_batching_space();
    SEQUANT_ASSERT(!batch_indices.empty());
    [[maybe_unused]] auto batch_space =
        get_default_context().index_space_registry()->retrieve(L"z");
    // assumes that there are no more than one type of batch space
    for ([[maybe_unused]] const auto& idx : batch_indices) {
      SEQUANT_ASSERT(idx.space() == batch_space);
    }

    const auto op_info = build_op_info(cre_spaces, ann_spaces, dep);
    const auto t = tensor_generator(op_info.creidxs, op_info.annidxs,
                                    batch_indices, op_info.opsymm);

    return ex<Constant>(rational{1, op_info.mult}) * t *
           ex<NormalOperator<S>>(cre(op_info.creidxs), ann(op_info.annidxs),
                                 get_default_context().vacuum());
  }

  /// @tparam TensorGenerator callable with signature
  /// `TensorGenerator(range<Index>, range<Index>, range<Index>, Symmetry)` that
  /// returns a Tensor with the respective bra/cre, ket/ann, and batch indices
  /// and of the given symmetry
  /// @param[in] creators creator IndexSpaces as an initializer list
  /// @param[in] annihilators annihilator IndexSpaces as an initializer list
  /// @param[in] batch_indices batch indices as an initializer list
  /// @param[in] tensor_generator the callable that generates the tensor
  /// @param[in] csv whether to use dependent indices
  template <typename TensorGenerator>
  static ExprPtr make(std::initializer_list<IndexSpace::Type> creators,
                      std::initializer_list<IndexSpace::Type> annihilators,
                      std::initializer_list<Index> batch_indices,
                      TensorGenerator&& tensor_generator,
                      UseDepIdx csv = UseDepIdx::None) {
    IndexSpaceContainer cre_vec(creators.begin(), creators.end());
    IndexSpaceContainer ann_vec(annihilators.begin(), annihilators.end());
    IndexContainer batchidx_vec(batch_indices.begin(), batch_indices.end());
    return OpMaker::make(cre_vec, ann_vec, batchidx_vec,
                         std::forward<TensorGenerator>(tensor_generator), csv);
  }

 protected:
  OpType op_;
  IndexSpaceContainer cre_spaces_;
  IndexSpaceContainer ann_spaces_;
  std::optional<IndexContainer> batch_indices_ = std::nullopt;

  OpMaker(OpType op);

  [[nodiscard]] auto ncreators() const { return cre_spaces_.size(); };
  [[nodiscard]] auto nannihilators() const { return ann_spaces_.size(); };
};

extern template class OpMaker<Statistics::FermiDirac>;
extern template class OpMaker<Statistics::BoseEinstein>;

/// \tparam QuantumNumbers a sequence of quantum numbers, must be
/// default-initializable
template <typename QuantumNumbers, Statistics S>
class Operator : public Operator<void, S> {
  using this_type = Operator<QuantumNumbers, S>;
  using base_type = Operator<void, S>;

 protected:
  Operator();

 public:
  /// @brief Constructs an operator with the given label and tensor form and
  /// quantum number action
  /// @param label_generator a function that generates the label for the
  /// operator
  /// @param tensor_form_generator a function that generates the tensor form of
  /// the operator
  /// @param qn_action a function that modifies the quantum numbers
  Operator(std::function<std::wstring_view()> label_generator,
           std::function<ExprPtr()> tensor_form_generator,
           std::function<void(QuantumNumbers&)> qn_action);

  /// @brief Constructs an operator with the given label and tensor form and
  /// quantum number action
  /// @param label_generator a function that generates the label for the
  /// operator
  /// @param tensor_form_generator a function that generates the tensor form of
  /// the operator
  /// @param qn_action a function that modifies the quantum numbers
  /// @param batch_idx_rank the rank of the batch index, must be non-zero
  Operator(std::function<std::wstring_view()> label_generator,
           std::function<ExprPtr()> tensor_form_generator,
           std::function<void(QuantumNumbers&)> qn_action,
           size_t batch_idx_rank);

  /// @brief Constructs an operator with the given label and tensor form and
  /// quantum number action
  /// @param label_generator a function that generates the label for the
  /// operator
  /// @param tensor_form_generator a function that generates the tensor form of
  /// the operator
  /// @param qn_action a function that modifies the quantum numbers
  /// @param batch_ordinals the unique, sorted ordinals of the batch indices
  Operator(std::function<std::wstring_view()> label_generator,
           std::function<ExprPtr()> tensor_form_generator,
           std::function<void(QuantumNumbers&)> qn_action,
           const container::svector<std::size_t>& batch_ordinals);

  virtual ~Operator();

  /// evaluates the result of applying this operator to \p qns
  /// \param qns the quantum numbers of the state to which this operator is
  /// applied; if not provided, the default-constructed \c QuantumNumbers are
  /// used \return the quantum numbers after applying this operator to \p qns
  QuantumNumbers operator()(const QuantumNumbers& qns = {}) const;

  /// evaluates the result of applying this operator to initializer-list-encoded
  /// \p qns
  /// \param qns the quantum numbers of the state to which this operator
  /// is applied; if not provided, the default-constructed \c QuantumNumbers are
  /// used
  /// \return the quantum numbers after applying this operator to \p qns
  template <typename I, typename = std::enable_if_t<std::is_integral_v<I>>>
  QuantumNumbers operator()(std::initializer_list<I> qns) const {
    QuantumNumbers result(qns);
    this->apply_to(result);
    return result;
  }

  /// evaluates the result of applying this operator to \p qns
  /// \param[in,out] qns the quantum numbers of the state to which this operator
  /// is applied; on return contains the quantum numbers after applying this
  /// operator
  /// \return a reference to `*this`
  virtual QuantumNumbers& apply_to(QuantumNumbers& qns) const;

  bool static_less_than(const Expr& that) const override;

  bool commutes_with_atom(const Expr& that) const override;

  void adjoint() override;

  /// @brief returns the batch ordinals if any
  std::optional<container::svector<std::size_t>> batch_ordinals() const {
    return batch_ordinals_;
  }

 private:
  std::function<void(QuantumNumbers&)> qn_action_;

  bool is_adjoint_ = false;

  std::optional<container::svector<std::size_t>> batch_ordinals_ = std::nullopt;

  bool less_than_rank_of(const this_type& that) const;

  Expr::type_id_type type_id() const override;

  ExprPtr clone() const override;

  std::wstring to_latex() const override;

  Expr::hash_type memoizing_hash() const override;

  bool static_equal(const Expr& other) const override;

};  // class Operator

extern template class Operator<qns_t, Statistics::FermiDirac>;
extern template class Operator<qns_t, Statistics::BoseEinstein>;

inline namespace op {
namespace tensor {

using mbpt::nₕ;
using mbpt::nₚ;

// clang-format off
/// @brief `k`-body contribution to the "generic" Hamiltonian (in normal order relative to the default vacuum)
/// @param[in] k the rank of the particle interactions; only `k<=2` is
/// supported
// clang-format on
ExprPtr H_(std::size_t k);

/// @brief Total Hamiltonian including up to `k`-body interactions
/// @param[in] k the maximum rank of the particle interactions; only `k<=2` is
/// supported
ExprPtr H(std::size_t k = 2);

/// @brief Fock operator implied one-body operator, optional explicit
/// construction requires user to specify the IndexSpace corresponding to all
/// orbitals which may have non-zero density.
ExprPtr F(bool use_tensor = true, IndexSpace reference_occupied = {L"", 0});

/// A general operator of rank \p K
ExprPtr θ(std::size_t K);

/// Makes particle-conserving excitation operator of rank \p K based on the
/// defined context
ExprPtr T_(std::size_t K);

/// Makes sum of particle-conserving excitation operators of all ranks up to \p
/// K based on the defined context
ExprPtr T(std::size_t K, bool skip1 = false);

/// Makes particle-conserving deexcitation operator of rank \p K based on the
/// defined context
ExprPtr Λ_(std::size_t K);

/// Makes sum of particle-conserving deexcitation operators of all ranks up to
/// \p K based on the defined context
ExprPtr Λ(std::size_t K, bool skip1 = false);

/// @brief Makes generic right-hand replacement operator
/// @param na number of annihilators
/// @param nc number of creators
/// @param cre_space IndexSpace on which creators act
/// @param ann_space IndexSpace on which annihilators act
ExprPtr R_(
    nann na, ncre nc,
    const cre<IndexSpace>& cre_space = cre(get_particle_space(Spin::any)),
    const ann<IndexSpace>& ann_space = ann(get_hole_space(Spin::any)));

/// @brief Makes generic excitation operator
/// @param np number of particle creators
/// @param nh number of hole creators
ExprPtr R_(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(R_);

/// @brief Makes generic left-hand replacement operator
/// @param na number of annihilators
/// @param nc number of creators
/// @param cre_space IndexSpace on which creators act
/// @param ann_space IndexSpace on which annihilators act
ExprPtr L_(
    nann na, ncre nc,
    const cre<IndexSpace>& cre_space = cre(get_hole_space(Spin::any)),
    const ann<IndexSpace>& ann_space = ann(get_particle_space(Spin::any)));

/// @brief Makes generic deexcitation operator
/// @param np number of particle annihilators
/// @param nh number of hole annihilators
ExprPtr L_(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(L_);

// clang-format off
/// makes projector onto excited bra (if \p np > 0 && \p nh > 0) or ket (if \p np < 0 && \p nh <0) manifold
/// @param np number of particle creators (if > 0) or annihilators (< 0)
/// @param nh number of hole creators (if > 0) or annihilators (< 0); if omitted, will use \p np
/// @note if using spin-free basis, only supports particle-symmetric operators `K = Kh = Kp`, returns `S(-K)`
/// else supports particle non-conserving operators and returns `A(-np, -nh)`
// clang-format on
ExprPtr P(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(P);

// clang-format off
/// @brief makes generic bra/ket-antisymmetric excitation (if \p nh > 0 && \p np > 0) or deexcitation (if \p nh < 0 && \p np < 0) operator
/// @param np number of particle creators (if > 0) or annihilators (< 0)
/// @param nh number of hole creators (if > 0) or annihilators (< 0); if omitted, will use \p np
/// (default is to set \p np to \p nh)
/// @note supports particle non-conserving operators
// clang-format on
ExprPtr A(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(A);

/// @brief makes generic particle-symmetric excitation (if \p K > 0) or
/// deexcitation (if \p K < 0) operator of rank `|K|`
ExprPtr S(std::int64_t K);

/// @brief Makes perturbation operator
/// @param R rank of the perturbation operator
/// @param params OpParams for operator construction. Default: order=1
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Hʼ(std::size_t R, const OpParams& params = {.order = 1});

/// @brief Makes perturbed particle-conserving excitation operator
/// @param K rank of the excitation operator
/// @param params OpParams for operator construction. Default: order=1
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Tʼ_(std::size_t K, const OpParams& params = {.order = 1});

/// @brief Makes sum of perturbed particle-conserving excitation operators
/// @param K rank up to which the sum is to be formed
/// @param params OpParams for operator construction. Default: order=1,
/// skip1=false
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Tʼ(std::size_t K,
           const OpParams& params = {.order = 1, .skip1 = false});

/// @brief Makes perturbed particle-conserving deexcitation operator
/// @param K rank of the deexcitation operator
/// @param params OpParams for operator construction. Default: order=1
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Λʼ_(std::size_t K, const OpParams& params = {.order = 1});

/// @brief Makes sum of perturbed particle-conserving deexcitation operators
/// @param K rank up to which the sum is to be formed
/// @param params OpParams for operator construction. Default: order=1,
/// skip1=false
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Λʼ(std::size_t K,
           const OpParams& params = {.order = 1, .skip1 = false});
}  // namespace tensor
}  // namespace op

inline namespace op {
// clang-format off
/// @brief `k`-body contribution to the "generic" Hamiltonian (in normal order relative to the default vacuum)
/// @param[in] k the rank of the particle interactions; only `k<=2` is
/// supported
// clang-format on
ExprPtr H_(std::size_t k);

/// @brief Total Hamiltonian including up to `k`-body interactions
/// @param[in] k the maximum rank of the particle interactions; only `k<=2` is
/// supported
ExprPtr H(std::size_t k = 2);

/// @brief Fock operator implied one-body operator, optional explicit
/// construction requires user to specify the IndexSpace corresponding to all
/// orbitals which may have non-zero density.
ExprPtr F(bool use_tensor = true, IndexSpace reference_occupied = {L"", 0});

/// A general operator of rank \p K
ExprPtr θ(std::size_t K);

/// Makes particle-conserving excitation operator of rank \p K
ExprPtr T_(std::size_t K);

/// Makes sum of particle-conserving excitation operators of all ranks up to \p
/// K
ExprPtr T(std::size_t K, bool skip1 = false);

/// Makes particle-conserving deexcitation operator of rank \p K
ExprPtr Λ_(std::size_t K);

/// Makes sum of particle-conserving deexcitation operators of all ranks up to
/// \p K
ExprPtr Λ(std::size_t K, bool skip1 = false);

/// @brief Makes generic excitation operator
/// @param na number of annihilators
/// @param nc number of creators
/// @param cre_space IndexSpace on which creators act
/// @param ann_space IndexSpace on which annihilators act
ExprPtr R_(
    nann na, ncre nc,
    const cre<IndexSpace>& cre_space = cre(get_particle_space(Spin::any)),
    const ann<IndexSpace>& ann_space = ann(get_hole_space(Spin::any)));

/// @brief Makes generic excitation operator
/// @param np number of particle creators
/// @param nh number of hole creators
ExprPtr R_(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(R_);

/// @brief Makes generic deexcitation operator
/// @param na number of annihilators
/// @param nc number of creators
/// @param cre_space IndexSpace on which creators act
/// @param ann_space IndexSpace on which annihilators act
ExprPtr L_(
    nann na, ncre nc,
    const cre<IndexSpace>& cre_space = cre(get_hole_space(Spin::any)),
    const ann<IndexSpace>& ann_space = ann(get_particle_space(Spin::any)));

/// @brief Makes generic deexcitation operator
/// @param np number of particle annihilators
/// @param nh number of hole annihilators
ExprPtr L_(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(L_);

/// @brief Makes sum of generic right-hand replacement operators up to max rank
/// @param na number of annihilators
/// @param nc number of creators
/// @param cre_space IndexSpace on which creators act
/// @param ann_space IndexSpace on which annihilators act
/// @return `R_(na,nc) + R_(na-1,nc-1) + ...`
ExprPtr R(nann na, ncre nc,
          const cre<IndexSpace>& cre_space = cre(get_particle_space(Spin::any)),
          const ann<IndexSpace>& ann_space = ann(get_hole_space(Spin::any)));

/// @brief Makes sum of generic excitation operators up to max rank
/// @param np max number of particle creators
/// @param nh max number of hole creators
/// @return `R_(np,nh) + R_(np-1,nh-1) + ...`
ExprPtr R(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(R);

/// @brief Makes sum of generic "left-hand" replacement operators up to max rank
/// @param na number of annihilators
/// @param nc number of creators
/// @param cre_space IndexSpace on which creators act
/// @param ann_space IndexSpace on which annihilators act
/// @return `L_(na,nc) + L_(na-1,nc-1) + ...`
ExprPtr L(
    nann na, ncre nc,
    const cre<IndexSpace>& cre_space = cre(get_hole_space(Spin::any)),
    const ann<IndexSpace>& ann_space = ann(get_particle_space(Spin::any)));

/// @brief Makes sum of deexcitation operators up to max rank
/// @param np max number of particle annihilators
/// @param nh max number of hole annihilators
/// @return `L_(np,nh) + L_(np-1,nh-1) + ...`
ExprPtr L(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(L);

// clang-format off
/// makes projector onto excited bra (if \p np > 0 && \p nh > 0) or ket (if \p np < 0 && \p nh <0) manifold
/// @param np number of particle creators (if > 0) or annihilators (< 0)
/// @param nh number of hole creators (if > 0) or annihilators (< 0); if omitted, will use \p np
/// @note if using spin-free basis, only supports particle-symmetric operators `K = Kh = Kp`, returns `S(-K)`
/// else supports particle non-conserving operators and returns `A(-np, -nh)`
// clang-format on
ExprPtr P(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(P);

// clang-format off
/// @brief makes generic bra/ket-antisymmetric excitation (if \p nh > 0 && \p np > 0) or deexcitation (if \p nh < 0 && \p np < 0) operator
/// @param np number of particle creators (if > 0) or annihilators (< 0)
/// @param nh number of hole creators (if > 0) or annihilators (< 0); if omitted, will use \p np
/// @note supports particle non-conserving operators
// clang-format on
ExprPtr A(nₚ np, nₕ nh);
DEFINE_SINGLE_SIGNED_ARGUMENT_OP_VARIANT(A);

/// @brief makes generic particle-symmetric excitation (if \p K > 0) or
/// deexcitation (if \p K < 0) operator of rank `|K|`
ExprPtr S(std::int64_t K);

/// @brief Makes perturbation operator
/// @param R rank of the perturbation operator
/// @param params OpParams for operator construction. Default: order=1
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Hʼ(std::size_t R, const OpParams& params = {.order = 1});

/// @brief Makes perturbed particle-conserving excitation operator from OpParams
/// @param K rank of the excitation operator
/// @param params OpParams for operator construction. Default: order=1
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Tʼ_(std::size_t K, const OpParams& params = {.order = 1});

/// @brief Makes sum of perturbed particle-conserving excitation operators
/// @param K rank up to which the sum is to be formed
/// @param params OpParams for operator construction. Default: order=1,
/// skip1=false
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Tʼ(std::size_t K,
           const OpParams& params = {.order = 1, .skip1 = false});

/// @brief Makes perturbed particle-conserving deexcitation operator
/// @param K rank of the deexcitation operator
/// @param params OpParams for operator construction. Default: order=1
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Λʼ_(std::size_t K, const OpParams& params = {.order = 1});

/// @brief Makes sum of perturbed particle-conserving deexcitation operators
/// @param K rank up to which the sum is to be formed
/// @param params OpParams for operator construction. Default: order=1,
/// skip1=false
/// @pre `params.order==1`, only first order perturbation is supported now
/// @pre If batching is used, ISR must contain batching space
ExprPtr Λʼ(std::size_t K,
           const OpParams& params = {.order = 1, .skip1 = false});

/// @brief computes the quantum number change effected by a given Operator or
/// Operator Product when applied to the vacuum state
/// @param expr the operator or operator product whose quantum number change is
/// to be computed
/// @return the quantum number change effected by \p expr
/// @pre \p expr must be an Operator or an Operator Product
qns_t apply_to_vac(const ExprPtr& expr);

/// @brief Checks if a given Operator or Operator Product can change the quantum
/// numbers from \p source_qns to \p target_qns
/// @param op_or_op_product the operator or operator product to check
/// @param target_qns the target quantum numbers
/// @param source_qns the source quantum numbers
bool can_change_qns(const ExprPtr& op_or_op_product, const qns_t& target_qns,
                    const qns_t& source_qns = {});

/// @brief Checks if a given Operator or Operator Product can raise the vacuum
/// quantum numbers up to rank \p k
bool raises_vacuum_up_to_rank(const ExprPtr& op_or_op_product,
                              const unsigned long k);

/// @brief Checks if a given Operator or Operator Product can lower quantum
/// numbers of rank \p down up to vacuum
bool lowers_rank_or_lower_to_vacuum(const ExprPtr& op_or_op_product,
                                    const unsigned long k);

/// @brief Checks if a given Operator or Operator Product can raise the vacuum
/// quantum numbers to rank \p k
bool raises_vacuum_to_rank(const ExprPtr& op_or_op_product,
                           const unsigned long k);

/// @brief Checks if a given Operator or Operator Product can lower quantum
/// numbers of rank \p k down to vacuum
bool lowers_rank_to_vacuum(const ExprPtr& op_or_op_product,
                           const unsigned long k);

#include <SeQuant/domain/mbpt/vac_av.hpp>
}  // namespace op
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_OP_HPP
