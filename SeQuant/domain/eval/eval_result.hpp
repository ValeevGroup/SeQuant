#ifndef SEQUANT_EVAL_RESULT_HPP
#define SEQUANT_EVAL_RESULT_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/hash.hpp>

#include <btas/btas.h>
#include <tiledarray.h>
#include <range/v3/view.hpp>

#include <any>
#include <memory>
#include <utility>

namespace sequant::eval {

namespace {
std::logic_error invalid_operand(
    std::string_view msg = "Invalid operand for binary op") noexcept {
  return std::logic_error{msg.data()};
}

std::logic_error unimplemented_method(std::string_view msg) noexcept {
  using namespace std::string_literals;
  return std::logic_error{"Not implemented in this derived class: "s +
                          msg.data()};
}

template <typename T>
struct Annot {
  explicit Annot(std::array<std::any, 3> const& a)
      : lannot(std::any_cast<T>(a[0])),
        rannot(std::any_cast<T>(a[1])),
        this_annot(std::any_cast<T>(a[2])) {}
  T const lannot;
  T const rannot;
  T const this_annot;
};

}  // namespace

namespace {

using perm_type = container::svector<size_t>;

template <typename F,
          typename = std::enable_if_t<std::is_invocable_v<F, perm_type const&>>>
void symmetric_permutation(size_t half_rank, F&& symmetrizer) noexcept {
  using ranges::views::concat;
  // using ranges::views::iota;
  using ranges::views::transform;

  // this vector contains indices from 0...rank/2 where rank is that
  // of the tensor being symmetrized
  //
  // Caveat:
  // creating perm_vec this way is not allowed in gcc-11,
  // sth to do with container::svector (boost::svector)
  // clang-format off
  // auto perm_vec = iota(size_t{0}, half_rank) | ranges::to<perm_type>;
  // clang-format on
  //
  auto perm_vec = perm_type(half_rank);
  for (auto i = 0; i < half_rank; ++i) perm_vec[i] = i;

  auto add_half_rank = [half_rank](auto x) { return x + half_rank; };
  do {
    auto const total_perm =
        concat(perm_vec, perm_vec | transform(add_half_rank)) |
        ranges::to<perm_type>;
    std::forward<F>(symmetrizer)(total_perm);
  } while (std::next_permutation(std::begin(perm_vec), std::end(perm_vec)));
}
template <typename F, typename = std::enable_if_t<
                          std::is_invocable_v<F, double, perm_type const&>>>

void antisymmetric_permutation(size_t half_rank, F&& antisymmetrizer) noexcept {
  using phase_type = double;

  using ranges::views::concat;
  //  using ranges::views::iota;

  int bra_parity = 0;
  //
  // Caveat:
  // creating bra_perm_vec this way is not allowed in gcc-11,
  // sth to do with container::svector (boost::svector)
  // clang-format off
  // auto bra_perm_vec = iota(size_t{0}, half_rank) | ranges::to<perm_type>;
  // clang-format on
  //
  auto bra_perm_vec = perm_type(half_rank);
  for (auto i = 0; i < half_rank; ++i) bra_perm_vec[i] = i;

  do {
    int ket_parity = 0;
    // same problem as with bra_perm_vec
    // clang-format off
    // auto ket_perm_vec = iota(half_rank, 2 * half_rank) | ranges::to<perm_type>;
    // clang-format on
    auto ket_perm_vec = perm_type(half_rank);
    for (auto i = 0; i < half_rank; ++i) ket_perm_vec[i] = i + half_rank;
    do {
      phase_type const phase_factor =
          (bra_parity + ket_parity) % 2 == 0 ? 1 : -1;
      auto const perm_vec =
          concat(bra_perm_vec, ket_perm_vec) | ranges::to<perm_type>;
      std::forward<F>(antisymmetrizer)(phase_factor, perm_vec);
    } while (next_permutation_parity(ket_parity,
                                     std::begin(ket_perm_vec),  //
                                     std::end(ket_perm_vec)));
  } while (next_permutation_parity(bra_parity,
                                   std::begin(bra_perm_vec),  //
                                   std::end(bra_perm_vec)));
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
    // Which is a potentailly narrowing conversion leading to
    // integral overflow. Hence, the values in the returned
    // container are mixed negative and positive integers (long type)
    //
    return static_cast<long>(sequant::hash::value(idx.label()));
  });
}

template <typename... Args>
auto symmetrize_ta(TA::DistArray<Args...> const& arr) {
  auto result = TA::DistArray<Args...>{arr.world(), arr.trange()};
  result.fill(0);
  size_t rank = arr.trange().rank();
  auto const lannot = ords_to_annot(ranges::views::iota(size_t{0}, rank));

  auto symmetrizer = [&result, &lannot, &arr](auto const& permutation) {
    auto const rannot = ords_to_annot(permutation);
    result(lannot) += arr(rannot);
  };

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "[EVAL] symmetrizing rank-" << rank << " tensor" << std::endl;
#endif

  symmetric_permutation(rank / 2, symmetrizer);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

template <typename... Args>
auto antisymmetrize_ta(TA::DistArray<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.trange().rank();
  auto const lannot = ords_to_annot(iota(size_t{0}, rank));

  auto result = TA::DistArray<Args...>{arr.world(), arr.trange()};
  result.fill(0);

  auto antisymmetrizer = [&result, &lannot, &arr](auto phase,
                                                  auto const& permutation) {
    auto const rannot = ords_to_annot(permutation);
    result(lannot) += phase * arr(rannot);
  };

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "[EVAL] antisymmetrizing rank-" << rank << " tensor"
            << std::endl;
#endif

  antisymmetric_permutation(rank / 2, antisymmetrizer);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

template <typename... Args>
auto symmetrize_btas(btas::Tensor<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.rank();
  // Caveat:
  // clang-format off
  // auto const lannot = iota(size_t{0}, rank) | ranges::to<perm_type>;
  // clang-format on
  auto const lannot = [rank]() {
    auto p = perm_type(rank);
    for (auto i = 0; i < rank; ++i) p[i] = i;
    return p;
  }();

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  auto symmetrizer = [&result, &lannot, &arr](auto const& permutation) {
    auto const& rannot = permutation;
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, rannot);
    result += temp;
  };

  symmetric_permutation(rank / 2, symmetrizer);

  return result;
}

template <typename... Args>
auto antisymmetrize_btas(btas::Tensor<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.rank();
  // Caveat:
  // auto const lannot = iota(size_t{0}, rank) | ranges::to<perm_type>;
  //
  auto const lannot = [rank]() {
    auto p = perm_type(rank);
    for (auto i = 0; i < rank; ++i) p[i] = i;
    return p;
  }();

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  auto antisymmetrizer = [&result, &lannot, &arr](auto phase,
                                                  auto const& permutation) {
    auto const& rannot = permutation;
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, rannot);
    btas::scal(phase, temp);
    result += temp;
  };

  antisymmetric_permutation(rank / 2, antisymmetrizer);

  return result;
}

}  // namespace

struct EvalResult;

using ERPtr = std::shared_ptr<EvalResult>;

template <typename T, typename... Args>
ERPtr eval_result(Args&&... args) noexcept {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

class EvalResult {
 public:
  using id_t = size_t;

  virtual ~EvalResult() noexcept = default;

  template <typename T>
  [[nodiscard]] bool is() const noexcept {
    return this->type_id() == id_for_type<std::decay_t<T>>();
  }

  template <typename T>
  [[nodiscard]] T const& as() const {
    assert(this->is<std::decay_t<T>>());
    return static_cast<T const&>(*this);
  }

  [[nodiscard]] virtual ERPtr sum(EvalResult const&,
                                  std::array<std::any, 3> const&) const = 0;

  [[nodiscard]] virtual ERPtr prod(EvalResult const&,
                                   std::array<std::any, 3> const&) const = 0;

  [[nodiscard]] virtual ERPtr permute(std::array<std::any, 2> const&) const = 0;

  virtual void add_inplace(EvalResult const&) = 0;

  [[nodiscard]] virtual ERPtr symmetrize() const = 0;

  [[nodiscard]] virtual ERPtr antisymmetrize() const = 0;

  [[nodiscard]] bool has_value() const noexcept;

  template <typename T>
  [[nodiscard]] T& get() {
    assert(has_value());
    return *std::any_cast<T>(&value_);
  }

  template <typename T>
  [[nodiscard]] T const& get() const {
    return const_cast<EvalResult&>(*this).get<T>();
  }

 protected:
  template <typename T>
  [[nodiscard]] T& get_ref() {
    assert(has_value());
    return *std::any_cast<T>(&value_);
  }

  template <typename T>
  explicit EvalResult(T arg) noexcept
      : value_{std::make_any<T>(std::move(arg))} {}

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
/// T is numeric type such as double
///
template <typename T>
class EvalConstant final : public EvalResult {
 public:
  using EvalResult::id_t;

  explicit EvalConstant(T v) noexcept : EvalResult{std::move(v)} {}

  [[nodiscard]] T value() const noexcept { return get<T>(); }

  [[nodiscard]] ERPtr sum(EvalResult const& other,
                          std::array<std::any, 3> const&) const override {
    if (other.is<EvalConstant<T>>()) {
      auto const& o = other.as<EvalConstant<T>>();
      return eval_result<EvalConstant<T>>(value() + o.value());
    } else {
      throw invalid_operand();
    }
  }

  [[nodiscard]] ERPtr prod(
      EvalResult const& other,
      std::array<std::any, 3> const& maybe_empty) const override {
    if (other.is<EvalConstant<T>>()) {
      auto const& o = other.as<EvalConstant<T>>();
      return eval_result<EvalConstant<T>>(value() * o.value());
    } else {
      return other.prod(*this, maybe_empty);
    }
  }

  [[nodiscard]] ERPtr permute(std::array<std::any, 2> const&) const override {
    throw unimplemented_method("permute");
  }

  void add_inplace(EvalResult const& other) override {
    auto& val = get<T>();
    val += other.get<T>();
  }

  [[nodiscard]] ERPtr symmetrize() const override {
    throw unimplemented_method("symmetrize");
  }

  [[nodiscard]] ERPtr antisymmetrize() const override {
    throw unimplemented_method("antisymmetrize");
  }

 private:
  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<EvalConstant<T>>();
  }
};

template <typename T>
class EvalTensorTA final : public EvalResult {
 public:
  using EvalResult::id_t;
  using numeric_type = std::decay_t<typename T::numeric_type>;

  explicit EvalTensorTA(T arr) : EvalResult{std::move(arr)} {}

 private:
  using annot_wrap = Annot<std::string>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<EvalTensorTA<T>>();
  }

  [[nodiscard]] ERPtr sum(EvalResult const& other,
                          std::array<std::any, 3> const& annot) const override {
    assert(other.is<EvalTensorTA<T>>());
    auto const a = annot_wrap{annot};

    T result;
    result(a.this_annot) = get<T>()(a.lannot) + other.get<T>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<EvalTensorTA<T>>(std::move(result));
  }

  [[nodiscard]] ERPtr prod(
      EvalResult const& other,
      std::array<std::any, 3> const& annot) const override {
    auto const a = annot_wrap{annot};

    if (other.is<EvalConstant<numeric_type>>()) {
      auto result = get<T>();
      result(a.this_annot) = other.get<numeric_type>() * result(a.lannot);
      decltype(result)::wait_for_lazy_cleanup(result.world());
      return eval_result<EvalTensorTA<T>>(std::move(result));
    }

    assert(other.is<EvalTensorTA<T>>());

    if (a.this_annot.empty())  // DOT product
      return eval_result<EvalConstant<numeric_type>>(
          TA::dot(get<T>()(a.lannot), other.get<T>()(a.rannot)));

    T result;
    result(a.this_annot) = get<T>()(a.lannot) * other.get<T>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<EvalTensorTA<T>>(std::move(result));
  }

  [[nodiscard]] ERPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    T result;
    result(post_annot) = get<T>()(pre_annot);
    T::wait_for_lazy_cleanup(result.world());
    return eval_result<EvalTensorTA<T>>(std::move(result));
  }

  void add_inplace(EvalResult const& other) override {
    auto& t = get_ref<T>();
    auto const& o = other.get<T>();

    assert(t.trange() == o.trange());
    auto ann = TA::detail::dummy_annotation(t.trange().rank());
    t(ann) += o(ann);
    T::wait_for_lazy_cleanup(t.world());
  }

  [[nodiscard]] ERPtr symmetrize() const override {
    return eval_result<EvalTensorTA<T>>(symmetrize_ta(get<T>()));
  }

  [[nodiscard]] ERPtr antisymmetrize() const override {
    return eval_result<EvalTensorTA<T>>(antisymmetrize_ta(get<T>()));
  }
};

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

  [[nodiscard]] ERPtr prod(
      EvalResult const& other,
      std::array<std::any, 3> const& annot) const override {
    auto const a = annot_wrap{annot};

    if (other.is<EvalConstant<numeric_type>>()) {
      T result;
      btas::permute(get<T>(), a.lannot, result, a.this_annot);
      btas::scal(other.as<EvalConstant<numeric_type>>().value(), result);
      return eval_result<EvalTensorBTAS<T>>(std::move(result));
    }

    assert(other.is<EvalTensorBTAS<T>>());

    if (a.this_annot.empty()) {
      T rres;
      btas::permute(other.get<T>(), a.rannot, rres, a.lannot);
      return eval_result<EvalConstant<numeric_type>>(btas::dot(get<T>(), rres));
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
    auto& t = get_ref<T>();
    auto const& o = other.get<T>();
    assert(t.range() == o.range());
    t += o;
  }

  [[nodiscard]] ERPtr symmetrize() const override {
    return eval_result<EvalTensorBTAS<T>>(symmetrize_btas(get<T>()));
  }

  [[nodiscard]] ERPtr antisymmetrize() const override {
    return eval_result<EvalTensorBTAS<T>>(antisymmetrize_btas(get<T>()));
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_RESULT_HPP
