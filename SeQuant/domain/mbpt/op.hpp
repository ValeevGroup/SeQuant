//
// Created by Eduard Valeyev on 2019-03-26.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_HPP
#define SEQUANT_DOMAIN_MBPT_OP_HPP

#include <SeQuant/domain/mbpt/fwd.hpp>

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

#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/range/primitives.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>
#include <range/v3/view/zip.hpp>

#include <algorithm>
#include <array>
#include <cassert>
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
  invalid       //!< invalid operator
};

/// maps operator types to their labels
inline const std::map<OpType, std::wstring> optype2label{
    {OpType::h, L"h"},
    {OpType::f, L"f"},
    {OpType::f̃, L"f̃"},
    {OpType::g, L"g"},
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
    {OpType::δ, L"δ"},
    {OpType::h_1, L"h¹"},
    {OpType::t_1, L"t¹"},
    {OpType::λ_1, L"λ¹"}};

/// maps operator labels to their types
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

struct default_qns_tag {};

// clang-format off
/// Tracks changes in \c N quantum numbers

/// Implements the concept of a quantum number change; this is useful for
/// tracking the quantum numbers of a many-body operator, such as the number of particles,
/// the number of quasiparticles, the number of ops (creators/annihilators) in each subspace, etc.
/// For example, to operator products expressed in normal order with respect to physical vacuum it is sufficient to track
/// the number of creators and annihilators; For the fermi vacuum case, the number of creators and annihilators in each
/// subspace becomes important. the number of ops is tracked for each base space (determined by the IndexSpaceRegistry object in Context).
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

  const std::size_t size() const {
    if (get_default_context().vacuum() == Vacuum::Physical) {
      return 2;
    } else if (get_default_context().vacuum() == Vacuum::SingleProduct) {
      auto isr = get_default_context().index_space_registry();
      const auto& isr_base_spaces = isr->base_spaces();
      assert(isr_base_spaces.size() > 0);
      return isr_base_spaces.size() * 2;
    } else {
      throw std::logic_error("unknown Vacuum type");
    }
  }

  /// initializes all values with zeroes
  QuantumNumberChange() {
    this->resize(this->size());
    assert(this->base().size() != 0);
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
    assert(i.size() == size());
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
    assert(i.size() == size());
#ifndef NDEBUG
    if (std::find_if(i.begin(), i.end(),
                     [](const auto& ii) { return ii.size() != 2; }) != i.end())
      throw std::invalid_argument(
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
    assert(i.size() == size());
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
    assert(size() > 0);
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
/// @note The choice of computational basis can be controlled by the default Formalism:
/// - if `get_default_formalism().sum_over_uocc() == SumOverUocc::Complete` IndexSpace::complete_unoccupied will be used instead of IndexSpace::active_unoccupied
/// - if `get_default_formalism().csv() == CSVFormalism::CSV` will use cluster-specific (e.g., PNO) unoccupied indices
/// @warning Tensor \f$ T \f$ will be antisymmetrized if `get_default_context().spbasis() == SPBasis::spinorbital`, else it will be particle-symmetric; the latter is only valid if # of bra and ket indices coincide.
/// @internal bless the maker and his water
// clang-format on
template <Statistics S>
class OpMaker {
 public:
  /// @param[in] op the operator type
  /// @param[in] bras the bra indices/creators
  /// @param[in] kets the ket indices/annihilators
  OpMaker(OpType op, std::initializer_list<IndexSpace> bras,
          std::initializer_list<IndexSpace> kets);

  /// @param[in] op the operator type
  /// @param[in] bras the bra indices/creators
  /// @param[in] kets the ket indices/annihilators
  template <typename IndexSpaceTypeRange1, typename IndexSpaceTypeRange2>
  OpMaker(OpType op, IndexSpaceTypeRange1&& bras, IndexSpaceTypeRange2&& kets)
      : op_(op),
        bra_spaces_(bras.begin(), bras.end()),
        ket_spaces_(kets.begin(), kets.end()) {
    assert(nbra() > 0 || nket() > 0);
  }

  OpMaker(OpType op, std::size_t Nbra, std::size_t Nket,
          IndexSpace particle_space = get_default_context()
                                          .index_space_registry()
                                          ->particle_space(Spin::any),
          IndexSpace hole_space = get_default_context()
                                      .index_space_registry()
                                      ->hole_space(Spin::any));

  // construct a particle conserving operator when constructing this way.
  OpMaker(OpType op, std::size_t nparticle);

  enum class UseDepIdx {
    /// bra/cre indices depend on ket
    Bra,
    /// ket/ann indices depend on bra
    Ket,
    /// use plain indices
    None
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

  /// @tparam TensorGenerator callable with signature
  /// `TensorGenerator(range<Index>, range<Index>, Symmetry)` that returns a
  /// Tensor with the respective bra and ket indices and of the given symmetry
  /// @param[in] bras the bra indices/creators
  /// @param[in] kets the ket indices/annihilators
  /// @param[in] tensor_generator the callable that generates the tensor
  /// @param[in] dep whether to use dependent indices
  template <typename TensorGenerator>
  static ExprPtr make(const container::svector<IndexSpace>& bras,
                      const container::svector<IndexSpace>& kets,
                      TensorGenerator&& tensor_generator,
                      UseDepIdx dep = UseDepIdx::None) {
    const bool symm =
        get_default_context().spbasis() ==
        SPBasis::spinorbital;  // antisymmetrize if spin-orbital basis
    const auto dep_bra = dep == UseDepIdx::Bra;
    const auto dep_ket = dep == UseDepIdx::Ket;

    // not sure what it means to use nonsymmetric operator if nbra != nket
    if (!symm) assert(ranges::size(bras) == ranges::size(kets));

    auto make_idx_vector = [](const auto& spacetypes) {
      container::svector<Index> result;
      const auto n = spacetypes.size();
      result.reserve(n);
      for (size_t i = 0; i != n; ++i) {
        auto space = spacetypes[i];
        result.push_back(Index::make_tmp_index(space));
      }
      return result;
    };

    auto make_depidx_vector = [](const auto& spacetypes, auto&& protoidxs) {
      const auto n = spacetypes.size();
      container::svector<Index> result;
      result.reserve(n);
      for (size_t i = 0; i != n; ++i) {
        auto space = spacetypes[i];
        result.push_back(Index::make_tmp_index(space, protoidxs, true));
      }
      return result;
    };

    container::svector<Index> braidxs, ketidxs;
    if (dep_bra) {
      ketidxs = make_idx_vector(kets);
      braidxs = make_depidx_vector(bras, ketidxs);
    } else if (dep_ket) {
      braidxs = make_idx_vector(bras);
      ketidxs = make_depidx_vector(kets, braidxs);
    } else {
      braidxs = make_idx_vector(bras);
      ketidxs = make_idx_vector(kets);
    }

    using ranges::size;
    const auto mult = symm ? factorial(size(bras)) * factorial(size(kets))
                           : factorial(size(bras));
    const auto opsymm = symm ? (S == Statistics::FermiDirac ? Symmetry::antisymm
                                                            : Symmetry::symm)
                             : Symmetry::nonsymm;
    return ex<Constant>(rational{1, mult}) *
           tensor_generator(braidxs, ketidxs, opsymm) *
           ex<NormalOperator<S>>(/* creators */ braidxs,
                                 /* annihilators */ ketidxs,
                                 get_default_context().vacuum());
  }

  template <typename TensorGenerator>
  static ExprPtr make(std::initializer_list<IndexSpace::Type> bras,
                      std::initializer_list<IndexSpace::Type> kets,
                      TensorGenerator&& tensor_generator,
                      UseDepIdx csv = UseDepIdx::None) {
    container::svector<IndexSpace> bra_vec(bras.begin(), bras.end());
    container::svector<IndexSpace> ket_vec(kets.begin(), kets.end());
    return OpMaker::make(bra_vec, ket_vec,
                         std::forward<TensorGenerator>(tensor_generator), csv);
  }

 protected:
  OpType op_;
  container::svector<IndexSpace> bra_spaces_;
  container::svector<IndexSpace> ket_spaces_;

  OpMaker(OpType op);

  const auto nbra() const { return bra_spaces_.size(); };
  const auto nket() const { return ket_spaces_.size(); };
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
  Operator(std::function<std::wstring_view()> label_generator,
           std::function<ExprPtr()> tensor_form_generator,
           std::function<void(QuantumNumbers&)> qn_action);

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

  bool is_adjoint_ = false;

 private:
  std::function<void(QuantumNumbers&)> qn_action_;

  bool less_than_rank_of(const this_type& that) const;

  Expr::type_id_type type_id() const override;

  ExprPtr clone() const override;

  std::wstring to_latex() const override;

  Expr::hash_type memoizing_hash() const override;

};  // class Operator

extern template class Operator<qns_t, Statistics::FermiDirac>;
extern template class Operator<qns_t, Statistics::BoseEinstein>;

inline namespace op {
namespace tensor {
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
ExprPtr Λ(std::size_t K);

/// @brief Makes generic excitation operator
/// @param nbra number of operators in bra space
/// @param nket number of operators in ket space
/// @param particle_space IndexSpace corresponding to the particle space
/// @param hole_space IndexSpace corresponding to the hole space
ExprPtr R_(
    std::size_t nbra, std::size_t nket,
    IndexSpace particle_space =
        get_default_context().index_space_registry()->particle_space(Spin::any),
    IndexSpace hole_space =
        get_default_context().index_space_registry()->hole_space(Spin::any));

/// @brief Makes generic deexcitation operator
/// @param nbra number of operators in bra space
/// @param nket number of operators in ket space
/// @param particle_space IndexSpace corresponding to the particle space
/// @param hole_space IndexSpace corresponding to the hole space
ExprPtr L_(
    std::size_t nbra, std::size_t nket,
    IndexSpace particle_space =
        get_default_context().index_space_registry()->particle_space(Spin::any),
    IndexSpace hole_space =
        get_default_context().index_space_registry()->hole_space(Spin::any));

ExprPtr P(std::int64_t K);

ExprPtr A(std::int64_t K);

ExprPtr S(std::int64_t K);

/// @brief Makes perturbation operator of rank \p R
/// \param order order of perturbation
/// \pre `order==1`, only first order perturbation is supported now
ExprPtr H_pt(std::size_t order, std::size_t R);

/// @brief Makes perturbed particle-conserving excitation operator of rank \p K
/// @param order order of perturbation
/// @pre `order==1`, only first order perturbation is supported now
ExprPtr T_pt_(std::size_t order, std::size_t K);

/// @brief Makes sum of perturbed particle-conserving excitation operators up to
/// rank \p K
/// @param order order of perturbation
/// @param skip1 if true, skips single excitations
/// @pre `order==1`, only first order perturbation is supported now
ExprPtr T_pt(std::size_t order, std::size_t K, bool skip1 = false);

/// @brief Makes perturbed particle-conserving deexcitation operator of rank \p
/// K
/// @param order order of perturbation
/// @pre `order==1`, only first order perturbation is supported now
ExprPtr Λ_pt_(std::size_t order, std::size_t K);

/// @brief Makes sum of perturbed particle-conserving deexcitation operators up
/// to rank \p K
/// @param order order of perturbation
/// @param skip1 if true, skips single deexcitations
/// @pre `order==1`, only first order perturbation is supported now
ExprPtr Λ_pt(std::size_t order, std::size_t K, bool skip1 = false);
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

/// Makes particle-conserving excitation operator of rank \p K
ExprPtr T_(std::size_t K);

/// Makes sum of particle-conserving excitation operators of all ranks up to \p
/// K
ExprPtr T(std::size_t K, bool skip1 = false);

/// Makes particle-conserving deexcitation operator of rank \p K
ExprPtr Λ_(std::size_t K);

/// Makes sum of particle-conserving deexcitation operators of all ranks up to
/// \p K
ExprPtr Λ(std::size_t K);

/// @brief Makes generic excitation operator
/// @param nbra number of operators in bra space
/// @param nket number of operators in ket space
/// @param particle_space IndexSpace corresponding to the particle space
/// @param hole_space IndexSpace corresponding to the hole space
ExprPtr R_(
    std::size_t nbra, std::size_t nket,
    IndexSpace particle_space =
        get_default_context().index_space_registry()->particle_space(Spin::any),
    IndexSpace hole_space =
        get_default_context().index_space_registry()->hole_space(Spin::any));

/// @brief Makes generic deexcitation operator
/// @param nbra number of operators in bra space
/// @param nket number of operators in ket space
/// @param particle_space IndexSpace corresponding to the particle space
/// @param hole_space IndexSpace corresponding to the hole space

ExprPtr L_(
    std::size_t nbra, std::size_t nket,
    IndexSpace particle_space =
        get_default_context().index_space_registry()->particle_space(Spin::any),
    IndexSpace hole_space =
        get_default_context().index_space_registry()->hole_space(Spin::any));

/// @brief Makes sum of generic excitation operators
/// @param nbra number of operators in bra space
/// @param nket number of operators in ket space
/// @param particle_space IndexSpace corresponding to the particle space
/// @param hole_space IndexSpace corresponding to the hole space
ExprPtr R(
    std::size_t nbra, std::size_t nket,
    IndexSpace particle_space =
        get_default_context().index_space_registry()->particle_space(Spin::any),
    IndexSpace hole_space =
        get_default_context().index_space_registry()->hole_space(Spin::any));

/// @brief Makes sum of generic deexcitation operators
/// @param nbra number of operators in bra space
/// @param nket number of operators in ket space
/// @param particle_space IndexSpace corresponding to the particle space
/// @param hole_space IndexSpace corresponding to the hole space
ExprPtr L(
    std::size_t nbra, std::size_t nket,
    IndexSpace particle_space =
        get_default_context().index_space_registry()->particle_space(Spin::any),
    IndexSpace hole_space =
        get_default_context().index_space_registry()->hole_space(Spin::any));

ExprPtr P(std::int64_t K);

ExprPtr A(std::int64_t K);

ExprPtr S(std::int64_t K);

ExprPtr H_pt(std::size_t order, std::size_t R);

ExprPtr T_pt_(std::size_t order, std::size_t K);

ExprPtr T_pt(std::size_t order, std::size_t K, bool skip1 = false);

ExprPtr Λ_pt_(std::size_t order, std::size_t K);

ExprPtr Λ_pt(std::size_t order, std::size_t K, bool skip1 = false);

bool raises_vacuum_up_to_rank(const ExprPtr& op_or_op_product,
                              const unsigned long k);

bool lowers_rank_or_lower_to_vacuum(const ExprPtr& op_or_op_product,
                                    const unsigned long k);

bool raises_vacuum_to_rank(const ExprPtr& op_or_op_product,
                           const unsigned long k);

bool lowers_rank_to_vacuum(const ExprPtr& op_or_op_product,
                           const unsigned long k);
#include <SeQuant/domain/mbpt/vac_av.hpp>
}  // namespace op
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_OP_HPP
