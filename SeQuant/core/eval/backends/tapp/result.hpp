#ifndef SEQUANT_EVAL_BACKENDS_TAPP_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_TAPP_RESULT_HPP

#ifdef SEQUANT_HAS_TAPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/tapp/ops.hpp>
#include <SeQuant/core/eval/backends/tapp/tensor.hpp>
#include <SeQuant/core/eval/elementwise.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <complex>
#include <cstddef>
#include <cstdint>

#include <range/v3/view/concat.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

namespace sequant {

// implementation details of the TAPP result backend; prefer sequant::detail
// over an unnamed namespace in a header (see CppCoreGuidelines SF.21 / "Use
// unnamed namespaces in headers ... no" guidance)
namespace detail {

///
/// \brief Symmetrize a TAPPTensor by summing over all symmetric particle
///        permutations.
///
/// \param arr The tensor to be symmetrized.
/// \pre The rank of the tensor must be even.
/// \return The symmetrized TAPPTensor.
///
template <typename T, typename Alloc>
auto column_symmetrize_tapp(TAPPTensor<T, Alloc> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.rank();

  if (rank % 2 != 0)
    throw Exception("This function only supports even-ranked tensors");

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;

  // lannot is identity ordering
  auto const lannot = perm | ranges::views::transform([](size_t v) {
                        return static_cast<int64_t>(v);
                      }) |
                      ranges::to<container::svector<int64_t>>;

  auto result = TAPPTensor<T, Alloc>(arr.extents());
  result.fill(T{0});

  auto call_back = [&result, &lannot, &arr, &perm = std::as_const(perm)]() {
    auto perm_annot = perm | ranges::views::transform([](size_t v) {
                        return static_cast<int64_t>(v);
                      }) |
                      ranges::to<container::svector<int64_t>>;
    TAPPTensor<T, Alloc> temp;
    tapp_ops::permute(arr, lannot, temp, perm_annot);
    result += temp;
  };

  auto const nparticles = rank / 2;
  symmetric_permutation(SymmetricParticleRange{perm.begin(),               //
                                               perm.begin() + nparticles,  //
                                               nparticles},
                        call_back);
  auto const nf = static_cast<double>(rational{1, factorial(nparticles)});
  tapp_ops::scal(static_cast<T>(nf), result);

  return result;
}

///
/// \brief Antisymmetrize a TAPPTensor over bra and ket index groups.
///
/// \param arr The tensor to be antisymmetrized.
/// \param bra_rank The rank of the bra indices.
/// \return The antisymmetrized TAPPTensor.
///
template <typename T, typename Alloc>
auto particle_antisymmetrize_tapp(TAPPTensor<T, Alloc> const& arr,
                                  size_t bra_rank) {
  using ranges::views::concat;
  using ranges::views::iota;
  size_t const rank = arr.rank();
  SEQUANT_ASSERT(bra_rank <= rank);
  size_t const ket_rank = rank - bra_rank;

  perm_t bra_perm = iota(size_t{0}, bra_rank) | ranges::to<perm_t>;
  perm_t ket_perm = iota(bra_rank, rank) | ranges::to<perm_t>;
  const auto lannot_perm = iota(size_t{0}, rank) | ranges::to<perm_t>;

  // Convert perm_t (svector<size_t>) to annot (svector<int64_t>)
  auto to_annot = [](perm_t const& p) {
    return p | ranges::views::transform([](size_t v) {
             return static_cast<int64_t>(v);
           }) |
           ranges::to<container::svector<int64_t>>;
  };

  auto const lannot = to_annot(lannot_perm);

  auto process_permutations = [&lannot, &to_annot](
                                  const TAPPTensor<T, Alloc>& input_arr,
                                  size_t range_rank, perm_t range_perm,
                                  const perm_t& other_perm, bool is_bra) {
    if (range_rank <= 1) return input_arr;
    TAPPTensor<T, Alloc> result{input_arr.extents()};
    result.fill(T{0});

    auto callback = [&](int parity) {
      const auto annot_perm =
          is_bra ? concat(range_perm, other_perm) | ranges::to<perm_t>()
                 : concat(other_perm, range_perm) | ranges::to<perm_t>();

      auto const annot = to_annot(annot_perm);

      T p_ = parity == 0 ? T{1} : T{-1};
      TAPPTensor<T, Alloc> temp;
      tapp_ops::permute(input_arr, lannot, temp, annot);
      tapp_ops::scal(p_, temp);
      result += temp;
    };

    antisymmetric_permutation(ParticleRange{range_perm.begin(), range_rank},
                              callback);
    return result;
  };

  // Process bra permutations first
  const auto ket_annot = ket_rank == 0 ? perm_t{} : ket_perm;
  auto result = process_permutations(arr, bra_rank, bra_perm, ket_annot, true);

  // Process ket permutations if needed
  const auto bra_annot = bra_rank == 0 ? perm_t{} : bra_perm;
  result = process_permutations(result, ket_rank, ket_perm, bra_annot, false);

  auto const nf = static_cast<double>(
      rational{1, factorial(bra_rank) * factorial(ket_rank)});
  tapp_ops::scal(static_cast<T>(nf), result);

  return result;
}

template <typename... Args>
inline void log_tapp(Args const&... args) noexcept {
  log_result("[TAPP] ", args...);
}

/// \return \p arr expanded over the modes of \p post that \p pre lacks (see
///         Result::broadcast), or, for a null \p arr, a tensor every element
///         of which is \p scalar
/// \param new_mode_extents the extent of each added mode, in the result's
///        mode order; a mode \p arr carries keeps \p arr's extent
template <typename T, typename Annot>
[[nodiscard]] T tapp_broadcast(
    T const* arr, typename T::numeric_type const& scalar, Annot const& pre,
    Annot const& post,
    container::svector<std::int64_t> const& new_mode_extents) {
  auto const pos = broadcast_operand_modes(pre, post);
  auto const rank = pos.size();
  SEQUANT_ASSERT(arr || pre.empty());
  container::svector<std::int64_t> extents(rank);
  for (std::size_t p = 0, k = 0; p < rank; ++p) {
    if (pos[p]) {
      extents[p] = arr->extents()[*pos[p]];
      continue;
    }
    SEQUANT_ASSERT(k < new_mode_extents.size());
    extents[p] = new_mode_extents[k++];
  }
  T result{extents};
  if (!arr) {
    for (auto& x : result) x = scalar;
    return result;
  }
  container::svector<std::size_t> uextents(rank), out_strides(rank),
      in_strides(rank);
  for (std::size_t p = 0; p < rank; ++p) {
    uextents[p] = static_cast<std::size_t>(extents[p]);
    out_strides[p] = static_cast<std::size_t>(result.strides()[p]);
    in_strides[p] =
        pos[p] ? static_cast<std::size_t>(arr->strides()[*pos[p]]) : 0;
  }
  broadcast_copy(result.data(), uextents, out_strides, arr->data(), in_strides);
  return result;
}

}  // namespace detail

///
/// \brief Result for a tensor value of TAPPTensor type.
/// \tparam T TAPPTensor type. Must be a specialization of TAPPTensor.
///
template <typename T>
class ResultTensorTAPP final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename T::numeric_type;

  explicit ResultTensorTAPP(T arr) : Result{std::move(arr)} {}

 private:
  using annot_t = container::svector<int64_t>;
  using annot_wrap = Annot<annot_t>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<ResultTensorTAPP<T>>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<ResultTensorTAPP<T>>());
    auto const a = annot_wrap{annot};

    detail::log_tapp(detail::ords_to_labels(a.lannot), " + ",
                     detail::ords_to_labels(a.rannot), " = ",
                     detail::ords_to_labels(a.this_annot), "\n");

    T lres, rres;
    tapp_ops::permute(get<T>(), a.lannot, lres, a.this_annot);
    tapp_ops::permute(other.get<T>(), a.rannot, rres, a.this_annot);
    lres += rres;
    return eval_result<ResultTensorTAPP<T>>(std::move(lres));
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      auto const scalar = other.as<ResultScalar<numeric_type>>().value();
      detail::log_tapp(detail::ords_to_labels(a.lannot), " * ", scalar, " = ",
                       detail::ords_to_labels(a.this_annot), "\n");

      T result;
      tapp_ops::permute(get<T>(), a.lannot, result, a.this_annot);
      tapp_ops::scal(scalar, result);
      return eval_result<ResultTensorTAPP<T>>(std::move(result));
    }

    SEQUANT_ASSERT(other.is<ResultTensorTAPP<T>>());

    if (a.this_annot.empty()) {
      // Full contraction -> scalar result (dot product)
      T rres;
      tapp_ops::permute(other.get<T>(), a.rannot, rres, a.lannot);
      numeric_type const d = tapp_ops::dot(get<T>(), rres);
      detail::log_tapp(detail::ords_to_labels(a.lannot), " * ",
                       detail::ords_to_labels(a.rannot), " = ", d, "\n");
      return eval_result<ResultScalar<numeric_type>>(d);
    }

    detail::log_tapp(detail::ords_to_labels(a.lannot), " * ",
                     detail::ords_to_labels(a.rannot), " = ",
                     detail::ords_to_labels(a.this_annot), "\n");

    T result;
    tapp_ops::contract(numeric_type{1},           //
                       get<T>(), a.lannot,        //
                       other.get<T>(), a.rannot,  //
                       numeric_type{0},           //
                       result, a.this_annot);
    return eval_result<ResultTensorTAPP<T>>(std::move(result));
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto pre = get<T>();
    tapp_ops::scal(numeric_type(factor), pre);
    return eval_result<ResultTensorTAPP<T>>(std::move(pre));
  }

  /// Deep copy: the backing tensor type owns its elements, so its copy
  /// constructor already produces an independently owned buffer.
  [[nodiscard]] ResultPtr clone() const override {
    return eval_result<ResultTensorTAPP<T>>(get<T>());
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<annot_t>(ann[0]);
    auto const post_annot = std::any_cast<annot_t>(ann[1]);

    detail::log_tapp(detail::ords_to_labels(pre_annot), " = ",
                     detail::ords_to_labels(post_annot), "\n");

    T result;
    tapp_ops::permute(get<T>(), pre_annot, result, post_annot);
    return eval_result<ResultTensorTAPP<T>>(std::move(result));
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    // T†{a;i} = conj(T{i;a}): bra/ket-swapped layout (operand annot in ann[0],
    // adjoint annot in ann[1]) plus elementwise conjugation. conj is a no-op
    // for real numeric_type, so it is elided there. (TAPP also exposes a native
    // TAPP_CONJUGATE element op; the explicit loop keeps this self-contained
    // and backend-version-independent.)
    auto const pre_annot = std::any_cast<annot_t>(ann[0]);
    auto const post_annot = std::any_cast<annot_t>(ann[1]);

    detail::log_tapp(detail::ords_to_labels(post_annot), " = adjoint(",
                     detail::ords_to_labels(pre_annot), ")\n");

    T result;
    tapp_ops::permute(get<T>(), pre_annot, result, post_annot);
    if constexpr (meta::is_complex_v<numeric_type>) {
      for (auto& x : result) x = std::conj(x);
    }
    return eval_result<ResultTensorTAPP<T>>(std::move(result));
  }

  [[nodiscard]] ResultPtr pow(rational const& exponent) const override {
    T result = get<T>();
    detail::ElementwisePow const op{exponent};
    for (auto& x : result) x = op(x);
    return eval_result<ResultTensorTAPP<T>>(std::move(result));
  }

  [[nodiscard]] ResultPtr broadcast(
      std::array<std::any, 2> const& ann,
      container::svector<ResultPtr> const& new_mode_ranges) const override {
    auto const pre_annot = std::any_cast<annot_t>(ann[0]);
    auto const post_annot = std::any_cast<annot_t>(ann[1]);
    detail::log_tapp(detail::ords_to_labels(post_annot), " = broadcast(",
                     detail::ords_to_labels(pre_annot), ")\n");
    return eval_result<ResultTensorTAPP<T>>(
        detail::tapp_broadcast(&get<T>(), numeric_type{}, pre_annot, post_annot,
                               new_mode_extents(new_mode_ranges)));
  }

  [[nodiscard]] ResultPtr broadcast_scalar(
      Result const& scalar, std::array<std::any, 2> const& ann,
      container::svector<ResultPtr> const& new_mode_ranges) const override {
    if (!scalar.is<ResultScalar<numeric_type>>())
      throw Exception("broadcast: the scalar's type is not this tensor's");
    auto const post_annot = std::any_cast<annot_t>(ann[1]);
    detail::log_tapp(detail::ords_to_labels(post_annot), " = broadcast(",
                     scalar.get<numeric_type>(), ")\n");
    return eval_result<ResultTensorTAPP<T>>(detail::tapp_broadcast<T>(
        nullptr, scalar.get<numeric_type>(), annot_t{}, post_annot,
        new_mode_extents(new_mode_ranges)));
  }

  void add_inplace(Result const& other) override {
    auto& t = get<T>();
    auto const& o = other.get<T>();
    SEQUANT_ASSERT(t.extents() == o.extents());

    auto const ann = detail::ords_to_annot(
        ranges::views::iota(size_t{0}, static_cast<size_t>(t.rank())));
    detail::log_tapp(ann, " += ", ann, "\n");

    t += o;
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultTensorTAPP<T>>(
        detail::column_symmetrize_tapp(get<T>()));
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<ResultTensorTAPP<T>>(
        detail::particle_antisymmetrize_tapp(get<T>(), bra_rank));
  }

 private:
  /// the extent of each rank-1 range Result::broadcast receives
  [[nodiscard]] static container::svector<std::int64_t> new_mode_extents(
      container::svector<ResultPtr> const& new_mode_ranges) {
    container::svector<std::int64_t> extents;
    extents.reserve(new_mode_ranges.size());
    for (auto const& r : new_mode_ranges) {
      if (!r->is<ResultTensorTAPP<T>>())
        throw Exception(
            "broadcast: a mode range is not of this tensor's type; the "
            "BackendArrayOps must build tensors of the evaluated type");
      SEQUANT_ASSERT(r->get<T>().rank() == 1);
      extents.push_back(r->get<T>().extents()[0]);
    }
    return extents;
  }

  [[nodiscard]] std::size_t size_in_bytes() const final {
    const auto& tensor = get<T>();
    return tensor.volume() * sizeof(typename T::value_type);
  }
};

}  // namespace sequant

#endif  // SEQUANT_HAS_TAPP

#endif  // SEQUANT_EVAL_BACKENDS_TAPP_RESULT_HPP
