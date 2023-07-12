#ifndef SEQUANT_EVAL_RESULT_HPP
#define SEQUANT_EVAL_RESULT_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/hash.hpp>

#include <btas/btas.h>
#include <tiledarray.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>
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

auto valid_particle_range = [](auto const& tpl) -> bool {
  using std::distance;
  auto [b1, b2, l] = tpl;
  return distance(b1, b1 + l) == distance(b2, b2 + l);
};

auto iter_pairs = [](auto&& tpl) {
  using ranges::views::iota;
  using ranges::views::transform;
  using std::get;

  auto b1 = get<0>(tpl);
  auto b2 = get<1>(tpl);
  auto l = get<2>(tpl);

  return iota(size_t{0}, l) | transform([b1, b2](auto i) {
           return IterPair{b1 + i, b2 + i};
         }) |
         ranges::to_vector;
};

}  // namespace

namespace {

using perm_t = container::svector<size_t>;
using particle_range_t = std::array<size_t, 3>;

template <typename F,
          std::enable_if_t<std::is_invocable_v<F, int>, bool> = true>
void antisymmetric_permutation(
    container::svector<
        std::tuple<perm_t::iterator, perm_t::iterator, size_t>> const& groups,
    F const& call_back) {
  using ranges::views::transform;

  auto const n = groups.size();
  if (n == 0) return;

  assert(ranges::all_of(groups, valid_particle_range));

  call_back(0);

  for (int i = n - 1; i >= 0; --i) {
    auto [bra_beg, ket_beg, len] = groups[i];

    auto bra_end = bra_beg + len;
    auto ket_end = ket_beg + len;

    int bra_p = 0;
    auto outer = 0;

    for (auto bra_yn = true; bra_yn;
         bra_yn = next_permutation_parity(bra_p, bra_beg, bra_end), ++outer) {
      auto inner = 0;
      int ket_p = 0;

      for (auto ket_yn = true; ket_yn;
           ket_yn = next_permutation_parity(ket_p, ket_beg, ket_end), ++inner) {
        if (!(outer == 0 && inner == 0)) call_back((bra_p + ket_p) % 2);
      }
    }
  }
}

template <typename F, std::enable_if_t<std::is_invocable_v<F>, bool> = true>
void symmetric_permutation(
    container::svector<
        std::tuple<perm_t::iterator, perm_t::iterator, size_t>> const& groups,
    F const& call_back) {
  using ranges::views::transform;

  auto const n = groups.size();
  if (n == 0) return;

  assert(ranges::all_of(groups, valid_particle_range));

  auto groups_vec = groups | transform(iter_pairs) | ranges::to_vector;

  call_back();

  // using reverse iterator (instead of indices) not allowed for some reason
  // iter from the end group
  for (int i = n - 1; i >= 0; --i) {
    auto beg = groups_vec[i].begin();
    auto end = groups_vec[i].end();
    auto yn = std::next_permutation(beg, end);
    for (; yn; yn = std::next_permutation(beg, end)) call_back();
  }
}

template <
    typename F,
    std::enable_if_t<std::is_invocable_v<F, int, perm_t const&>, bool> = true>
void antisymmetrize_backend(size_t rank,
                            container::svector<particle_range_t> const& groups,
                            F const& call_back) {
  using ranges::views::iota;
  auto perm = iota(size_t{0}, rank) | ranges::to<perm_t>;

  auto groups_vec = container::svector<
      std::tuple<perm_t::iterator, perm_t::iterator, size_t>>{};
  groups_vec.reserve(groups.size());
  auto beg = perm.begin();
  for (auto&& g : groups) {
    groups_vec.emplace_back(beg + g[0], beg + g[1], g[2]);
  }
  antisymmetric_permutation(
      groups_vec,
      [&call_back, &perm = std::as_const(perm)](int p) { call_back(p, perm); });
}

template <typename F,
          std::enable_if_t<std::is_invocable_v<F, perm_t const&>, bool> = true>
void symmetrize_backend(size_t rank,
                        container::svector<particle_range_t> const& groups,
                        F const& call_back) {
  using ranges::views::iota;
  auto perm = iota(size_t{0}, rank) | ranges::to<perm_t>;
  auto groups_vec = container::svector<
      std::tuple<perm_t::iterator, perm_t::iterator, size_t>>{};
  groups_vec.reserve(groups.size());
  auto beg = perm.begin();
  for (auto&& g : groups) {
    groups_vec.emplace_back(beg + g[0], beg + g[1], g[2]);
  }
  symmetric_permutation(
      groups_vec,
      [&call_back, &perm = std::as_const(perm)]() { call_back(perm); });
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
auto symmetrize_ta(TA::DistArray<Args...> const& arr,
                   container::svector<particle_range_t> const& groups) {
  using ranges::views::iota;

  auto result = TA::DistArray<Args...>{arr.world(), arr.trange()};
  result.fill(0);

  size_t const rank = arr.trange().rank();

  auto const lannot = ords_to_annot(iota(size_t{0}, rank));

  auto call_back = [&result, &lannot, &arr](perm_t const& perm) {
    auto const rannot = ords_to_annot(perm);
    result(lannot) += arr(rannot);
  };

  symmetrize_backend(rank, groups, call_back);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

template <typename... Args>
auto antisymmetrize_ta(
    TA::DistArray<Args...> const& arr,
    container::svector<particle_range_t> const& groups = {}) {
  using ranges::views::iota;
  using ranges::views::transform;
  using ranges::views::zip;

  size_t const rank = arr.trange().rank();

  auto result = TA::DistArray<Args...>(arr.world(), arr.trange());
  result.fill(0);

  auto const lannot = ords_to_annot(iota(size_t{0}, rank));

  auto call_back = [&lannot, &arr, &result](int p, perm_t const& perm) {
    typename decltype(result)::numeric_type p_ = p == 0 ? 1 : -1;
    result(lannot) += p_ * arr(ords_to_annot(perm));
  };

  antisymmetrize_backend(rank, groups, call_back);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

template <typename... Args>
auto symmetrize_btas(btas::Tensor<Args...> const& arr,
                     container::svector<particle_range_t> const& groups) {
  using ranges::views::iota;

  size_t const rank = arr.rank();
  // Caveat:
  // clang-format off
  // auto const lannot = iota(size_t{0}, rank) | ranges::to<perm_t>;
  // clang-format on
  auto const lannot = [rank]() {
    auto p = perm_t(rank);
    for (auto i = 0; i < rank; ++i) p[i] = i;
    return p;
  }();

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  auto call_back = [&result, &lannot, &arr](auto const& permutation) {
    auto const& rannot = permutation;
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, rannot);
    result += temp;
  };

  symmetrize_backend(rank, groups, call_back);
  return result;
}

template <typename... Args>
auto antisymmetrize_btas(
    btas::Tensor<Args...> const& arr,
    container::svector<particle_range_t> const& groups = {}) {
  using ranges::views::iota;

  size_t const rank = arr.rank();
  // Caveat:
  // auto const lannot = iota(size_t{0}, rank) | ranges::to<perm_t>;
  //
  auto const lannot = [rank]() {
    auto p = perm_t(rank);
    for (auto i = 0; i < rank; ++i) p[i] = i;
    return p;
  }();

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  auto call_back = [&result, &lannot, &arr](int p, perm_t const& perm) {
    typename decltype(result)::numeric_type p_ = p == 0 ? 1 : -1;
    auto const& rannot = perm;
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, rannot);
    btas::scal(p_, temp);
    result += temp;
  };

  antisymmetrize_backend(rank, groups, call_back);

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

  ///
  /// @note In std::array<std::any, 3> is expected to be [l,r,res] where
  ///       the elements are the annotations for left, right and result
  ///       respectively.
  ///
  [[nodiscard]] virtual ERPtr sum(EvalResult const&,
                                  std::array<std::any, 3> const&) const = 0;

  ///
  /// @note In std::array<std::any, 3> is expected to be [l,r,res] where
  ///       the elements are the annotations for left, right and result
  ///       respectively.
  ///
  [[nodiscard]] virtual ERPtr prod(EvalResult const&,
                                   std::array<std::any, 3> const&) const = 0;

  ///
  /// @note In std::array<std::any, 2> is expected to be [pre,post] where
  ///       the elements are the annotations for the eval result before
  ///       permutation and after permutation respectively.
  ///
  [[nodiscard]] virtual ERPtr permute(std::array<std::any, 2> const&) const = 0;

  virtual void add_inplace(EvalResult const&) = 0;

  ///
  /// @note vector<array<size_t,3>> represents list of particle
  ///       symmetry index groups. array<size_t, 3> is expected to be
  ///         [b1, b2, len]
  ///       where b1, and b2 are the zero-based positions of the tensor indices.
  ///       [b1, b1+len) and [b2, b2+len) are two ranges that will be permuted
  ///       simultaneously and at the equivalent positions.
  ///
  [[nodiscard]] virtual ERPtr symmetrize(
      container::svector<std::array<size_t, 3>> const&) const = 0;

  ///
  /// @note vector<array<size_t,3>> represents list of antisymmetric index
  ///       groups. array<size_t,3> is expected to be
  ///         [b1,b2,len]
  ///       where b1, and b2 are the zero-based positions of the tensor indices.
  ///       [b1, b1+len) and [b2, b2+len) are two ranges that will be permuted
  ///       by keeping track of the parity (even/odd)-ness of the total
  ///       permutation.
  ///
  [[nodiscard]] virtual ERPtr antisymmetrize(
      container::svector<std::array<size_t, 3>> const&) const = 0;

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

  [[nodiscard]] ERPtr symmetrize(
      container::svector<std::array<size_t, 3>> const&) const override {
    throw unimplemented_method("symmetrize");
  }

  [[nodiscard]] ERPtr antisymmetrize(
      container::svector<std::array<size_t, 3>> const&) const override {
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

  [[nodiscard]] ERPtr symmetrize(
      container::svector<std::array<size_t, 3>> const& groups) const override {
    return eval_result<EvalTensorTA<T>>(symmetrize_ta(get<T>(), groups));
  }

  [[nodiscard]] ERPtr antisymmetrize(
      container::svector<std::array<size_t, 3>> const& groups) const override {
    return eval_result<EvalTensorTA<T>>(antisymmetrize_ta(get<T>(), groups));
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

  [[nodiscard]] ERPtr symmetrize(
      container::svector<std::array<size_t, 3>> const& groups) const override {
    return eval_result<EvalTensorBTAS<T>>(symmetrize_btas(get<T>(), groups));
  }

  [[nodiscard]] ERPtr antisymmetrize(
      container::svector<std::array<size_t, 3>> const& groups) const override {
    return eval_result<EvalTensorBTAS<T>>(
        antisymmetrize_btas(get<T>(), groups));
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_RESULT_HPP
