#ifndef SEQUANT_EVAL_BACKENDS_TAPP_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_TAPP_RESULT_HPP

#ifdef SEQUANT_HAS_TAPP

#include <SeQuant/core/eval/backends/tapp/ops.hpp>
#include <SeQuant/core/eval/backends/tapp/tensor.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/utility/exception.hpp>

namespace sequant {

namespace {

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

}  // namespace

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
      T result;
      tapp_ops::permute(get<T>(), a.lannot, result, a.this_annot);
      tapp_ops::scal(other.as<ResultScalar<numeric_type>>().value(), result);
      return eval_result<ResultTensorTAPP<T>>(std::move(result));
    }

    SEQUANT_ASSERT(other.is<ResultTensorTAPP<T>>());

    if (a.this_annot.empty()) {
      // Full contraction -> scalar result (dot product)
      T rres;
      tapp_ops::permute(other.get<T>(), a.rannot, rres, a.lannot);
      return eval_result<ResultScalar<numeric_type>>(
          tapp_ops::dot(get<T>(), rres));
    }

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

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<annot_t>(ann[0]);
    auto const post_annot = std::any_cast<annot_t>(ann[1]);
    T result;
    tapp_ops::permute(get<T>(), pre_annot, result, post_annot);
    return eval_result<ResultTensorTAPP<T>>(std::move(result));
  }

  void add_inplace(Result const& other) override {
    auto& t = get<T>();
    auto const& o = other.get<T>();
    SEQUANT_ASSERT(t.extents() == o.extents());
    t += o;
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultTensorTAPP<T>>(column_symmetrize_tapp(get<T>()));
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<ResultTensorTAPP<T>>(
        particle_antisymmetrize_tapp(get<T>(), bra_rank));
  }

 private:
  [[nodiscard]] std::size_t size_in_bytes() const final {
    const auto& tensor = get<T>();
    return tensor.volume() * sizeof(typename T::value_type);
  }
};

}  // namespace sequant

#endif  // SEQUANT_HAS_TAPP

#endif  // SEQUANT_EVAL_BACKENDS_TAPP_RESULT_HPP
