#ifndef SEQUANT_EVAL_RESULT_HPP
#define SEQUANT_EVAL_RESULT_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>

#include <TiledArray/einsum/tiledarray.h>
#include <btas/btas.h>
#include <tiledarray.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include <any>
#include <memory>
#include <utility>

namespace sequant {

namespace {

[[maybe_unused]] std::logic_error invalid_operand(
    std::string_view msg = "Invalid operand for binary op") noexcept {
  return std::logic_error{msg.data()};
}

[[maybe_unused]] std::logic_error unimplemented_method(
    std::string_view msg) noexcept {
  using namespace std::string_literals;
  return std::logic_error{"Not implemented in this derived class: "s +
                          msg.data()};
}

///
/// \brief This class represents a triplet of annotations used in a tensor
///        contraction or summation.
///
/// \tparam T Annotation type eg. TA::DistArray takes std::string.
///
template <typename T>
struct Annot {
  explicit Annot(std::array<std::any, 3> const& a)
      : lannot(std::any_cast<T>(a[0])),
        rannot(std::any_cast<T>(a[1])),
        this_annot(std::any_cast<T>(a[2])) {}

  /// Annotation of the left operand.
  T const lannot;

  /// Annotation of the right operand.
  T const rannot;

  /// Annotation of the result.
  T const this_annot;
};

// It is an iterator type
template <typename It>
struct IterPair {
  It first, second;
  IterPair(It beg, It end) noexcept : first{beg}, second{end} {};
};

template <typename It>
void swap(IterPair<It>& l, IterPair<It>& r) {
  using std::iter_swap;
  std::iter_swap(l.first, r.first);
  std::iter_swap(l.second, r.second);
}

template <typename It>
bool operator<(IterPair<It> const& l, IterPair<It> const& r) noexcept {
  return *l.first < *r.first;
}

using perm_t = container::svector<size_t>;

struct SymmetricParticleRange {
  perm_t::iterator bra_beg;
  perm_t::iterator ket_beg;
  size_t nparticles;
};

struct AntisymmetricParticleRange {
  perm_t::iterator bra_beg;
  size_t bra_size;
  perm_t::iterator ket_beg;
  size_t ket_size;
};

inline bool valid_particle_range(SymmetricParticleRange const& rng) {
  using std::distance;
  auto bra_end = rng.bra_beg + rng.nparticles;
  auto ket_end = rng.ket_beg + rng.nparticles;
  return std::is_sorted(rng.bra_beg, bra_end) &&
         std::is_sorted(rng.ket_beg, ket_end) &&
         distance(rng.bra_beg, bra_end) == distance(rng.ket_beg, ket_end);
}

inline bool valid_particle_range(AntisymmetricParticleRange const& rng) {
  auto bra_end = rng.bra_beg + rng.bra_size;
  auto ket_end = rng.ket_beg + rng.ket_size;
  return std::is_sorted(rng.bra_beg, bra_end) &&
         std::is_sorted(rng.ket_beg, ket_end);
}

inline auto iter_pairs(SymmetricParticleRange const& rng) {
  using ranges::views::iota;
  using ranges::views::transform;

  return iota(size_t{0}, rng.nparticles)  //
         | transform([b = rng.bra_beg, k = rng.ket_beg](auto i) {
             return IterPair{b + i, k + i};
           });
}

///
/// \brief This function permutes the given range of antisymmetric particles
///        in-place and calls the given callback function with the parity of the
///        permutation as the argument.
///
/// \param rng The antisymmetric particle range.
///
/// \param call_back A function object which is called with the parity of the
///                  particle antisymsmetric permutation.
///

template <typename F, typename = std::enable_if_t<std::is_invocable_v<F, int>>>
void antisymmetric_permutation(AntisymmetricParticleRange const& rng,
                               F call_back) {
  struct {
    int bra;
    int ket;
  } parity{0, 0};

  auto bra_end = rng.bra_beg + rng.bra_size;
  auto ket_end = rng.ket_beg + rng.ket_size;

  for (auto bra_yn = true;                            //
       bra_yn;                                        //
       bra_yn = next_permutation_parity(parity.bra,   //
                                        rng.bra_beg,  //
                                        bra_end)) {
    for (auto ket_yn = true;                            //
         ket_yn;                                        //
         ket_yn = next_permutation_parity(parity.ket,   //
                                          rng.ket_beg,  //
                                          ket_end)) {
      call_back((parity.bra + parity.ket) % 2);
    }
  }
}

///
/// \brief This function permutes the given range of symmetric particles
///        in-place and calls the given callback function.
///
/// \param rng The symmetric particle range.
///
/// \param call_back A function object which is called after permutation.
///

template <typename F, typename = std::enable_if_t<std::is_invocable_v<F>>>
void symmetric_permutation(SymmetricParticleRange const& rng, F call_back) {
  auto ips = iter_pairs(rng) | ranges::to_vector;
  do {
    call_back();
  } while (std::next_permutation(ips.begin(), ips.end()));
}

template <typename RngOfOrdinals>
std::string ords_to_annot(RngOfOrdinals const& ords) {
  using ranges::views::intersperse;
  using ranges::views::join;
  using ranges::views::transform;
  auto to_str = [](auto x) { return std::to_string(x); };
  return ords | transform(to_str) | intersperse(std::string{","}) | join |
         ranges::to<std::string>;
}

///
/// \param bk iterable of Index objects.
/// \return vector of long-type hash values
///         of the labels of indices in \c bk
///
template <typename Iterable>
auto index_hash(Iterable const& bk) {
  return ranges::views::transform(bk, [](auto const& idx) {
    //
    // WARNING!
    // The BTAS expects index types to be long by default.
    // There is no straight-forward way to turn the default.
    // Hence, here we explicitly cast the size_t values to long
    // Which is a potentially narrowing conversion leading to
    // integral overflow. Hence, the values in the returned
    // container are mixed negative and positive integers (long type)
    //
    return static_cast<long>(sequant::hash::value(Index{idx}.label()));
  });
}

///
/// \brief This function implements the symmetrization of TA::DistArray.
///
/// \param arr The array to be symmetrized
///
/// \pre The rank of the array must be even
///
/// \return The symmetrized TA::DistArray.
///
template <typename... Args>
auto particle_symmetrize_ta(TA::DistArray<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.trange().rank();
  if (rank % 2 != 0)
    throw std::domain_error("This function only supports even-ranked tensors");

  TA::DistArray<Args...> result;

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;

  auto const lannot = ords_to_annot(perm);

  auto call_back = [&result, &lannot, &arr, &perm = std::as_const(perm)]() {
    auto const rannot = ords_to_annot(perm);
    if (result.is_initialized()) {
      result(lannot) += arr(rannot);
    } else {
      result(lannot) = arr(rannot);
    }
  };

  auto const nparticles = rank / 2;
  symmetric_permutation(SymmetricParticleRange{perm.begin(),               //
                                               perm.begin() + nparticles,  //
                                               nparticles},
                        call_back);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

///
/// \brief This function implements the antisymmetrization of TA::DistArray.
///
/// \param arr The array to be antisymmetrized.
///
/// \param bra_rank The rank of the bra indices
///
/// \return The antisymmetrized TA::DistArray.
///
template <typename... Args>
auto particle_antisymmetrize_ta(TA::DistArray<Args...> const& arr,
                                size_t bra_rank) {
  using ranges::views::iota;

  size_t const rank = arr.trange().rank();
  assert(bra_rank <= rank);

  TA::DistArray<Args...> result;

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;

  auto const lannot = ords_to_annot(perm);

  auto call_back = [&lannot, &arr, &result,
                    &perm = std::as_const(perm)](int parity) {
    typename decltype(result)::numeric_type p_ = parity == 0 ? 1 : -1;
    if (result.is_initialized())
      result(lannot) += p_ * arr(ords_to_annot(perm));
    else
      result(lannot) = p_ * arr(ords_to_annot(perm));
  };

  auto const ket_rank = rank - bra_rank;

  antisymmetric_permutation(
      AntisymmetricParticleRange{perm.begin(),             //
                                 bra_rank,                 //
                                 perm.begin() + bra_rank,  //
                                 ket_rank},
      call_back);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

///
/// \brief This function implements the symmetrization of btas::Tensor.
///
/// \param arr The tensor to be symmetrized.
///
/// \pre The rank of the tensor must be even.
///
/// \return The symmetrized btas::Tensor.
///
template <typename... Args>
auto particle_symmetrize_btas(btas::Tensor<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.rank();

  if (rank % 2 != 0)
    throw std::domain_error("This function only supports even-ranked tensors");

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;

  auto const lannot = perm;

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  auto call_back = [&result, &lannot, &arr, &perm = std::as_const(perm)]() {
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, perm);
    result += temp;
  };

  auto const nparticles = rank / 2;
  symmetric_permutation(SymmetricParticleRange{perm.begin(),               //
                                               perm.begin() + nparticles,  //
                                               nparticles},
                        call_back);

  return result;
}

///
/// \brief This function implements the antisymmetrization of btas::Tensor.
///
/// \param arr The tensor to be antisymmetrized
///
/// \param bra_rank The rank of the bra indices
///
/// \return The antisymmetrized btas::Tensor.
///
template <typename... Args>
auto particle_antisymmetrize_btas(btas::Tensor<Args...> const& arr,
                                  size_t bra_rank) {
  using ranges::views::iota;

  size_t const rank = arr.rank();
  assert(bra_rank <= rank);

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;
  auto const lannot = perm;

  auto call_back = [&lannot, &arr, &result,
                    &perm = std::as_const(perm)](int parity) {
    typename decltype(result)::numeric_type p_ = parity == 0 ? 1 : -1;
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, perm);
    btas::scal(p_, temp);
    result += temp;
  };

  auto const ket_rank = rank - bra_rank;

  antisymmetric_permutation(
      AntisymmetricParticleRange{perm.begin(),             //
                                 bra_rank,                 //
                                 perm.begin() + bra_rank,  //
                                 ket_rank},
      call_back);

  return result;
}

template <typename... Args>
inline void log_result(Args const&... args) noexcept {
#ifdef SEQUANT_EVAL_TRACE
  auto l = Logger::instance();
  if (l.log_level_eval > 1) write_log(l, args...);
#endif
}

template <typename... Args>
inline void log_ta(Args const&... args) noexcept {
#ifdef SEQUANT_EVAL_TRACE
  log_result("[TA] ", args...);
#endif
}

template <typename... Args>
inline void log_constant(Args const&... args) noexcept {
#ifdef SEQUANT_EVAL_TRACE
  log_result("[CONST] ", args...);
#endif
}

}  // namespace

/// TA::Tensor memory use logger
/// If TiledArray was configured with TA_TENSOR_MEM_PROFILE set this
/// prints the current use of memory by TA::Tensor objects in host memory space
/// to \p os .
/// \param world the world object to use for logging
/// \param label string to prepend to the profile
void log_ta_tensor_host_memory_use(madness::World& world,
                                   std::string_view label = "");

class EvalResult;

using ERPtr = std::shared_ptr<EvalResult>;

///
/// \brief Factory function for EvalResult objects.
///
/// \tparam T A concrete type derived from EvalResult.
/// \tparam Args The argument types for the constructor of T.
/// \param args The arguments for the constructor of T.
/// \return ERPtr object.
///
template <typename T, typename... Args>
ERPtr eval_result(Args&&... args) noexcept {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

///
/// \brief A type-erased class for the result of an evaluation. An object of
///        this class can represent a tensor (eg. TA::TArrayD,
///        btas::Tesnor<double>, etc.) or a scalar (eg. double, complex<double>
///        etc.).
///
class EvalResult {
 public:
  using id_t = size_t;

  virtual ~EvalResult() noexcept = default;

  ///
  /// \return Returns true if the concrete type of the object is T.
  ///
  template <typename T>
  [[nodiscard]] bool is() const noexcept {
    return this->type_id() == id_for_type<std::decay_t<T>>();
  }

  ///
  /// \return T const& where T is the concrete type of the object.
  ///
  template <typename T>
  [[nodiscard]] T const& as() const {
    assert(this->is<std::decay_t<T>>());
    return static_cast<T const&>(*this);
  }

  ///
  /// \brief Sum other EvalResult object with this object.
  ///
  /// @note In std::array<std::any, 3> is expected to be [l,r,res] where
  ///       the elements are the annotations for left, right and result
  ///       respectively.
  ///
  [[nodiscard]] virtual ERPtr sum(EvalResult const&,
                                  std::array<std::any, 3> const&) const = 0;

  ///
  /// \brief Perform product binary operation with this object and other.
  ///
  /// @note In std::array<std::any, 3> is expected to be [l,r,res] where
  ///       the elements are the annotations for left, right and result
  ///       respectively.
  ///
  [[nodiscard]] virtual ERPtr prod(EvalResult const&,
                                   std::array<std::any, 3> const&,
                                   TA::DeNest DeNestFlag) const = 0;

  ///
  /// \brief Permute this object according to the annotations in the argument.
  ///
  /// @note In std::array<std::any, 2> is expected to be [pre,post] where
  ///       the elements are the annotations for the eval result before
  ///       permutation and after permutation respectively.
  ///
  [[nodiscard]] virtual ERPtr permute(std::array<std::any, 2> const&) const = 0;

  ///
  /// \brief Add other EvalResult object into this object.
  ///
  virtual void add_inplace(EvalResult const&) = 0;

  ///
  /// \brief Particle symmetrize this eval result according to the list of index
  ///        groups.
  ///
  /// @note vector<array<size_t,3>> represents list of particle
  ///       symmetry index groups. array<size_t, 3> is expected to be
  ///         [b1, b2, len]
  ///       where b1, and b2 are the zero-based positions of the tensor indices.
  ///       [b1, b1+len) and [b2, b2+len) are two ranges that will be permuted
  ///       simultaneously and at the equivalent positions.
  ///
  [[nodiscard]] virtual ERPtr symmetrize() const = 0;

  ///
  /// \brief Particle antisymmetrize this eval result according to the list of
  ///        index groups.
  ///
  /// @note vector<array<size_t,3>> represents list of antisymmetric index
  ///       groups. array<size_t,3> is expected to be
  ///         [b1,b2,len]
  ///       where b1, and b2 are the zero-based positions of the tensor indices.
  ///       [b1, b1+len) and [b2, b2+len) are two ranges that will be permuted
  ///       by keeping track of the parity -- (even/odd)-ness -- of the total
  ///       permutation.
  ///
  [[nodiscard]] virtual ERPtr antisymmetrize(size_t bra_rank) const = 0;

  [[nodiscard]] bool has_value() const noexcept;

  ///
  /// \return Cast the type-erased data to the type \tparam T, and return a ref.
  ///
  template <typename T>
  [[nodiscard]] T& get() {
    assert(has_value());
    return *std::any_cast<T>(&value_);
  }

  ///
  /// \return Cast the type-erased data to the type \tparam T, and return a
  ///         const ref.
  ///
  template <typename T>
  [[nodiscard]] T const& get() const {
    return const_cast<EvalResult&>(*this).get<T>();
  }

 protected:
  template <typename T,
            typename = std::enable_if_t<!std::is_convertible_v<T, EvalResult>>>
  explicit EvalResult(T&& arg) noexcept
      : value_{std::make_any<std::decay_t<T>>(std::forward<T>(arg))} {}

  [[nodiscard]] virtual id_t type_id() const noexcept = 0;

  template <typename T>
  [[nodiscard]] static id_t id_for_type() noexcept {
    static id_t id = next_id();
    return id;
  }

 private:
  std::any value_;

  [[nodiscard]] static id_t next_id() noexcept;
};

///
/// \brief EvalResult for a constant or a variable value.
///
/// \tparam T numeric type of the constant value (eg. double, complex<double>,
///         etc.)
template <typename T>
class EvalScalar final : public EvalResult {
 public:
  using EvalResult::id_t;

  explicit EvalScalar(T v) noexcept : EvalResult{std::move(v)} {}

  [[nodiscard]] T value() const noexcept { return get<T>(); }

  [[nodiscard]] ERPtr sum(EvalResult const& other,
                          std::array<std::any, 3> const&) const override {
    if (other.is<EvalScalar<T>>()) {
      auto const& o = other.as<EvalScalar<T>>();
      auto s = value() + o.value();

      log_constant(value(), " + ", o.value(), " = ", s, "\n");

      return eval_result<EvalScalar<T>>(s);
    } else {
      throw invalid_operand();
    }
  }

  [[nodiscard]] ERPtr prod(EvalResult const& other,
                           std::array<std::any, 3> const& maybe_empty,
                           TA::DeNest DeNestFlag) const override {
    if (other.is<EvalScalar<T>>()) {
      auto const& o = other.as<EvalScalar<T>>();
      auto p = value() * o.value();

      log_constant(value(), " * ", o.value(), " = ", p, "\n");

      return eval_result<EvalScalar<T>>(value() * o.value());
    } else {
      auto maybe_empty_ = maybe_empty;
      std::swap(maybe_empty_[0], maybe_empty_[1]);
      return other.prod(*this, maybe_empty_, DeNestFlag);
    }
  }

  [[nodiscard]] ERPtr permute(std::array<std::any, 2> const&) const override {
    throw unimplemented_method("permute");
  }

  void add_inplace(EvalResult const& other) override {
    assert(other.is<EvalScalar<T>>());
    log_constant(value(), " += ", other.get<T>(), "\n");
    auto& val = get<T>();
    val += other.get<T>();
  }

  [[nodiscard]] ERPtr symmetrize() const override {
    throw unimplemented_method("symmetrize");
  }

  [[nodiscard]] ERPtr antisymmetrize(size_t bra_rank) const override {
    throw unimplemented_method("antisymmetrize");
  }

 private:
  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<EvalScalar<T>>();
  }
};

///
/// \brief EvalResult for a tensor value of TA::DistArray type.
/// \tparam ArrayT TA::DistArray type. Tile type of ArrayT is regular tensor of
///                scalars (not a tensor of tensors)
///
template <typename ArrayT, typename = std::enable_if_t<TA::detail::is_tensor_v<
                               typename ArrayT::value_type>>>
class EvalTensorTA final : public EvalResult {
 public:
  using EvalResult::id_t;
  using numeric_type = typename ArrayT::numeric_type;

  explicit EvalTensorTA(ArrayT arr) : EvalResult{std::move(arr)} {}

 private:
  using this_type = EvalTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  [[nodiscard]] ERPtr sum(EvalResult const& other,
                          std::array<std::any, 3> const& annot) const override {
    assert(other.is<this_type>());
    auto const a = annot_wrap{annot};

    log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    result(a.this_annot) =
        get<ArrayT>()(a.lannot) + other.get<ArrayT>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ERPtr prod(EvalResult const& other,
                           std::array<std::any, 3> const& annot,
                           TA::DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<EvalScalar<numeric_type>>()) {
      auto result = get<ArrayT>();
      auto scalar = other.get<numeric_type>();

      log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

      result(a.this_annot) = scalar * result(a.lannot);

      decltype(result)::wait_for_lazy_cleanup(result.world());
      return eval_result<this_type>(std::move(result));
    }

    if (a.this_annot.empty()) {
      // DOT product
      assert(other.is<this_type>());
      numeric_type d =
          TA::dot(get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot));
      ArrayT::wait_for_lazy_cleanup(get<ArrayT>().world());
      ArrayT::wait_for_lazy_cleanup(other.get<ArrayT>().world());

      log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      return eval_result<EvalScalar<numeric_type>>(d);
    }

    if (!other.is<this_type>()) {
      // potential T * ToT
      auto annot_swap = annot;
      std::swap(annot_swap[0], annot_swap[1]);
      return other.prod(*this, annot_swap, DeNestFlag);
    }

    // confirmed: other.is<this_type>() is true

    log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;

    result = TA::einsum(get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot),
                        a.this_annot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ERPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    log_ta(pre_annot, " = ", post_annot, "\n");

    ArrayT result;
    result(post_annot) = get<ArrayT>()(pre_annot);
    ArrayT::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  void add_inplace(EvalResult const& other) override {
    assert(other.is<this_type>());

    auto& t = get<ArrayT>();
    auto const& o = other.get<ArrayT>();

    assert(t.trange() == o.trange());
    auto ann = TA::detail::dummy_annotation(t.trange().rank());

    log_ta(ann, " += ", ann, "\n");

    t(ann) += o(ann);
    ArrayT::wait_for_lazy_cleanup(t.world());
  }

  [[nodiscard]] ERPtr symmetrize() const override {
    return eval_result<this_type>(particle_symmetrize_ta(get<ArrayT>()));
  }

  [[nodiscard]] ERPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<this_type>(
        particle_antisymmetrize_ta(get<ArrayT>(), bra_rank));
  }
};

template <typename ArrayT,
          typename = std::enable_if_t<
              TA::detail::is_tensor_of_tensor_v<typename ArrayT::value_type>>>
class EvalTensorOfTensorTA final : public EvalResult {
 public:
  using EvalResult::id_t;
  using numeric_type = typename ArrayT::numeric_type;

  explicit EvalTensorOfTensorTA(ArrayT arr) : EvalResult{std::move(arr)} {}

 private:
  using this_type = EvalTensorOfTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  using _inner_tensor_type = typename ArrayT::value_type::value_type;

  using compatible_regular_distarray_type =
      TA::DistArray<_inner_tensor_type, typename ArrayT::policy_type>;

  // Only @c that_type type is allowed for ToT * T computation
  using that_type = EvalTensorTA<compatible_regular_distarray_type>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  [[nodiscard]] ERPtr sum(EvalResult const& other,
                          std::array<std::any, 3> const& annot) const override {
    assert(other.is<this_type>());
    auto const a = annot_wrap{annot};

    log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    result(a.this_annot) =
        get<ArrayT>()(a.lannot) + other.get<ArrayT>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ERPtr prod(EvalResult const& other,
                           std::array<std::any, 3> const& annot,
                           TA::DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<EvalScalar<numeric_type>>()) {
      auto result = get<ArrayT>();
      auto scalar = other.get<numeric_type>();

      log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

      result(a.this_annot) = scalar * result(a.lannot);

      decltype(result)::wait_for_lazy_cleanup(result.world());
      return eval_result<this_type>(std::move(result));
    } else if (a.this_annot.empty()) {
      // DOT product
      assert(other.is<this_type>());
      numeric_type d =
          TA::dot(get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot));
      ArrayT::wait_for_lazy_cleanup(get<ArrayT>().world());
      ArrayT::wait_for_lazy_cleanup(other.get<ArrayT>().world());

      log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      return eval_result<EvalScalar<numeric_type>>(d);
    }

    log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    if (other.is<that_type>()) {
      // ToT * T -> ToT
      auto result =
          TA::einsum(get<ArrayT>()(a.lannot),
                     other.get<compatible_regular_distarray_type>()(a.rannot),
                     a.this_annot);
      return eval_result<this_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == TA::DeNest::True) {
      // ToT * ToT -> T
      auto result = TA::einsum<TA::DeNest::True>(
          get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot), a.this_annot);
      return eval_result<that_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == TA::DeNest::False) {
      // ToT * ToT -> ToT
      auto result = TA::einsum(get<ArrayT>()(a.lannot),
                               other.get<ArrayT>()(a.rannot), a.this_annot);
      return eval_result<this_type>(std::move(result));
    } else {
      throw invalid_operand();
    }
  }

  [[nodiscard]] ERPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    log_ta(pre_annot, " = ", post_annot, "\n");

    ArrayT result;
    result(post_annot) = get<ArrayT>()(pre_annot);
    ArrayT::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  void add_inplace(EvalResult const& other) override {
    assert(other.is<this_type>());

    auto& t = get<ArrayT>();
    auto const& o = other.get<ArrayT>();

    assert(t.trange() == o.trange());
    auto ann = TA::detail::dummy_annotation(t.trange().rank());

    log_ta(ann, " += ", ann, "\n");

    t(ann) += o(ann);
    ArrayT::wait_for_lazy_cleanup(t.world());
  }

  [[nodiscard]] ERPtr symmetrize() const override { return nullptr; }

  [[nodiscard]] ERPtr antisymmetrize(size_t bra_rank) const override {
    return nullptr;
  }
};

///
/// \brief EvalResult for a tensor value of btas::Tensor type.
/// \tparam T btas::Tensor type. Must be a specialization of btas::Tensor.
///
template <typename T>
class EvalTensorBTAS final : public EvalResult {
 public:
  using EvalResult::id_t;
  using numeric_type = typename T::numeric_type;

  explicit EvalTensorBTAS(T arr) : EvalResult{std::move(arr)} {}

 private:
  // TODO make it same as that used by EvalExprBTAS class from eval.hpp file
  using annot_t = container::svector<long>;
  using annot_wrap = Annot<annot_t>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<EvalTensorBTAS<T>>();
  }

  [[nodiscard]] ERPtr sum(EvalResult const& other,
                          std::array<std::any, 3> const& annot) const override {
    assert(other.is<EvalTensorBTAS<T>>());
    auto const a = annot_wrap{annot};

    T lres, rres;
    btas::permute(get<T>(), a.lannot, lres, a.this_annot);
    btas::permute(other.get<T>(), a.rannot, rres, a.this_annot);
    return eval_result<EvalTensorBTAS<T>>(lres + rres);
  }

  [[nodiscard]] ERPtr prod(EvalResult const& other,
                           std::array<std::any, 3> const& annot,
                           TA::DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<EvalScalar<numeric_type>>()) {
      T result;
      btas::permute(get<T>(), a.lannot, result, a.this_annot);
      btas::scal(other.as<EvalScalar<numeric_type>>().value(), result);
      return eval_result<EvalTensorBTAS<T>>(std::move(result));
    }

    assert(other.is<EvalTensorBTAS<T>>());

    if (a.this_annot.empty()) {
      T rres;
      btas::permute(other.get<T>(), a.rannot, rres, a.lannot);
      return eval_result<EvalScalar<numeric_type>>(btas::dot(get<T>(), rres));
    }

    T result;
    btas::contract(numeric_type{1},           //
                   get<T>(), a.lannot,        //
                   other.get<T>(), a.rannot,  //
                   numeric_type{0},           //
                   result, a.this_annot);
    return eval_result<EvalTensorBTAS<T>>(std::move(result));
  }

  [[nodiscard]] ERPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<annot_t>(ann[0]);
    auto const post_annot = std::any_cast<annot_t>(ann[1]);
    T result;
    btas::permute(get<T>(), pre_annot, result, post_annot);
    return eval_result<EvalTensorBTAS<T>>(std::move(result));
  }

  void add_inplace(EvalResult const& other) override {
    auto& t = get<T>();
    auto const& o = other.get<T>();
    assert(t.range() == o.range());
    t += o;
  }

  [[nodiscard]] ERPtr symmetrize() const override {
    return eval_result<EvalTensorBTAS<T>>(particle_symmetrize_btas(get<T>()));
  }

  [[nodiscard]] ERPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<EvalTensorBTAS<T>>(
        particle_antisymmetrize_btas(get<T>(), bra_rank));
  }
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_RESULT_HPP
