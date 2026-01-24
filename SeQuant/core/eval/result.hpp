#ifndef SEQUANT_EVAL_RESULT_HPP
#define SEQUANT_EVAL_RESULT_HPP

#include <SeQuant/core/eval/fwd.hpp>

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/macros.hpp>

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

struct ParticleRange {
  perm_t::iterator beg;
  size_t size;
};

inline bool valid_particle_range(SymmetricParticleRange const& rng) {
  using std::distance;
  auto bra_end = rng.bra_beg + rng.nparticles;
  auto ket_end = rng.ket_beg + rng.nparticles;
  return std::is_sorted(rng.bra_beg, bra_end) &&
         std::is_sorted(rng.ket_beg, ket_end) &&
         distance(rng.bra_beg, bra_end) == distance(rng.ket_beg, ket_end);
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
///                  particle antisymmetric permutation.
///
template <typename F, typename = std::enable_if_t<std::is_invocable_v<F, int>>>
void antisymmetric_permutation(ParticleRange const& rng, F call_back) {
  // if the range has 1 or no elements, there is no permutation
  if (rng.size <= 1) {
    call_back(0);
    return;
  }
  int parity = 0;
  auto end = rng.beg + rng.size;
  for (auto yn = true; yn; yn = next_permutation_parity(parity, rng.beg, end)) {
    call_back(parity);
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

template <typename... Args>
inline void log_result(Args const&... args) noexcept {
  auto& l = Logger::instance();
  if (l.eval.level > 1) write_log(l, args...);
}

template <typename... Args>
inline void log_constant(Args const&... args) noexcept {
  log_result("[CONST] ", args...);
}

}  // namespace

/******************************************************************************/

///
/// \brief Factory function for Result objects.
///
/// \tparam T A concrete type derived from Result.
/// \tparam Args The argument types for the constructor of T.
/// \param args The arguments for the constructor of T.
/// \return ResultPtr object.
///
template <typename T, typename... Args>
ResultPtr eval_result(Args&&... args) noexcept {
  return std::make_shared<T>(std::forward<Args>(args)...);
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

///
/// \brief A type-erased class for the result of an evaluation. An object of
///        this class can represent a tensor (eg. TA::DistArray,
///        btas::Tensor, etc.) or a scalar (eg. double, std::complex
///        etc.).
///
class Result {
 public:
  using id_t = size_t;

  virtual ~Result() noexcept = default;

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
    SEQUANT_ASSERT(this->is<std::decay_t<T>>());
    return static_cast<T const&>(*this);
  }

  ///
  /// \brief Sum other Result object with this object.
  ///
  /// @note In std::array<std::any, 3> is expected to be [l,r,res] where
  ///       the elements are the annotations for left, right and result
  ///       respectively.
  ///
  [[nodiscard]] virtual ResultPtr sum(Result const&,
                                      std::array<std::any, 3> const&) const = 0;

  ///
  /// \brief Perform product binary operation with this object and other.
  ///
  /// @note In std::array<std::any, 3> is expected to be [l,r,res] where
  ///       the elements are the annotations for left, right and result
  ///       respectively.
  ///
  [[nodiscard]] virtual ResultPtr prod(Result const&,
                                       std::array<std::any, 3> const&,
                                       DeNest DeNestFlag) const = 0;

  ///
  /// \brief Permute this object according to the annotations in the argument.
  ///
  /// @note In std::array<std::any, 2> is expected to be [pre,post] where
  ///       the elements are the annotations for the eval result before
  ///       permutation and after permutation respectively.
  ///
  [[nodiscard]] virtual ResultPtr permute(
      std::array<std::any, 2> const&) const = 0;

  ///
  /// \brief Add other Result object into this object.
  ///
  virtual void add_inplace(Result const&) = 0;

  ///
  /// \brief Particle symmetrize the eval result
  ///
  [[nodiscard]] virtual ResultPtr symmetrize() const = 0;

  ///
  /// \brief Particle antisymmetrize the eval result
  ///
  [[nodiscard]] virtual ResultPtr antisymmetrize(size_t bra_rank) const = 0;

  [[nodiscard]] bool has_value() const noexcept;

  [[nodiscard]] virtual ResultPtr mult_by_phase(std::int8_t) const = 0;

  ///
  /// \return Cast the type-erased data to the type \tparam T, and return a ref.
  ///
  template <typename T>
  [[nodiscard]] T& get() {
    SEQUANT_ASSERT(has_value());
    return *std::any_cast<T>(&value_);
  }

  ///
  /// \return Cast the type-erased data to the type \tparam T, and return a
  ///         const ref.
  ///
  template <typename T>
  [[nodiscard]] T const& get() const {
    SEQUANT_ASSERT(has_value());
    return *std::any_cast<const T>(&value_);
  }

  /// @return the size of the object in bytes
  [[nodiscard]] virtual std::size_t size_in_bytes() const = 0;

 protected:
  template <typename T,
            typename = std::enable_if_t<!std::is_convertible_v<T, Result>>>
  explicit Result(T&& arg) noexcept
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
/// \brief Result for a constant or a variable value.
///
/// \tparam T numeric type of the constant value (eg. double, complex<double>,
///         etc.)
template <typename T>
class ResultScalar final : public Result {
 public:
  using Result::id_t;

  explicit ResultScalar(T v) noexcept : Result{std::move(v)} {}

  [[nodiscard]] T value() const noexcept { return get<T>(); }

  [[nodiscard]] ResultPtr sum(Result const& other,
                              std::array<std::any, 3> const&) const override {
    if (other.is<ResultScalar<T>>()) {
      auto const& o = other.as<ResultScalar<T>>();
      auto s = value() + o.value();

      log_constant(value(), " + ", o.value(), " = ", s, "\n");

      return eval_result<ResultScalar<T>>(s);
    } else {
      throw invalid_operand();
    }
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& maybe_empty,
                               DeNest DeNestFlag) const override {
    if (other.is<ResultScalar<T>>()) {
      auto const& o = other.as<ResultScalar<T>>();
      auto p = value() * o.value();

      log_constant(value(), " * ", o.value(), " = ", p, "\n");

      return eval_result<ResultScalar<T>>(value() * o.value());
    } else {
      auto maybe_empty_ = maybe_empty;
      std::swap(maybe_empty_[0], maybe_empty_[1]);
      return other.prod(*this, maybe_empty_, DeNestFlag);
    }
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const&) const override {
    throw unimplemented_method("permute");
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<ResultScalar<T>>());
    log_constant(value(), " += ", other.get<T>(), "\n");
    auto& val = get<T>();
    val += other.get<T>();
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    throw unimplemented_method("symmetrize");
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    throw unimplemented_method("antisymmetrize");
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    return eval_result<ResultScalar<T>>(value() * T(factor));
  }

 private:
  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<ResultScalar<T>>();
  }

  [[nodiscard]] std::size_t size_in_bytes() const final { return sizeof(T); }
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_RESULT_HPP
