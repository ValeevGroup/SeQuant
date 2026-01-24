#ifndef SEQUANT_EVAL_BACKENDS_BTAS_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_BTAS_RESULT_HPP

#ifdef SEQUANT_HAS_BTAS

#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/math.hpp>

#include <btas/btas.h>

namespace sequant {

namespace {

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
auto column_symmetrize_btas(btas::Tensor<Args...> const& arr) {
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
  auto const nf = static_cast<double>(rational{1, factorial(nparticles)});
  btas::scal(nf, result);

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
  using ranges::views::concat;
  using ranges::views::iota;
  size_t const rank = arr.rank();
  SEQUANT_ASSERT(bra_rank <= rank);
  size_t const ket_rank = rank - bra_rank;

  perm_t bra_perm = iota(size_t{0}, bra_rank) | ranges::to<perm_t>;
  perm_t ket_perm = iota(bra_rank, rank) | ranges::to<perm_t>;
  const auto lannot = iota(size_t{0}, rank) | ranges::to<perm_t>;

  auto process_permutations = [&lannot](const btas::Tensor<Args...>& input_arr,
                                        size_t range_rank, perm_t range_perm,
                                        const perm_t& other_perm, bool is_bra) {
    if (range_rank <= 1) return input_arr;
    btas::Tensor<Args...> result{input_arr.range()};

    auto callback = [&](int parity) {
      const auto annot =
          is_bra ? concat(range_perm, other_perm) | ranges::to<perm_t>()
                 : concat(other_perm, range_perm) | ranges::to<perm_t>();

      typename decltype(result)::numeric_type p_ = parity == 0 ? 1 : -1;
      btas::Tensor<Args...> temp;
      btas::permute(input_arr, lannot, temp, annot);
      btas::scal(p_, temp);
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
  btas::scal(nf, result);

  return result;
}

}  // namespace

///
/// \brief Result for a tensor value of btas::Tensor type.
/// \tparam T btas::Tensor type. Must be a specialization of btas::Tensor.
///
template <typename T>
class ResultTensorBTAS final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename T::numeric_type;

  explicit ResultTensorBTAS(T arr) : Result{std::move(arr)} {}

 private:
  // TODO make it same as that used by EvalExprBTAS class from eval.hpp file
  using annot_t = container::svector<long>;
  using annot_wrap = Annot<annot_t>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<ResultTensorBTAS<T>>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<ResultTensorBTAS<T>>());
    auto const a = annot_wrap{annot};

    T lres, rres;
    btas::permute(get<T>(), a.lannot, lres, a.this_annot);
    btas::permute(other.get<T>(), a.rannot, rres, a.this_annot);
    return eval_result<ResultTensorBTAS<T>>(lres + rres);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      T result;
      btas::permute(get<T>(), a.lannot, result, a.this_annot);
      btas::scal(other.as<ResultScalar<numeric_type>>().value(), result);
      return eval_result<ResultTensorBTAS<T>>(std::move(result));
    }

    SEQUANT_ASSERT(other.is<ResultTensorBTAS<T>>());

    if (a.this_annot.empty()) {
      T rres;
      btas::permute(other.get<T>(), a.rannot, rres, a.lannot);
      return eval_result<ResultScalar<numeric_type>>(btas::dot(get<T>(), rres));
    }

    T result;
    btas::contract(numeric_type{1},           //
                   get<T>(), a.lannot,        //
                   other.get<T>(), a.rannot,  //
                   numeric_type{0},           //
                   result, a.this_annot);
    return eval_result<ResultTensorBTAS<T>>(std::move(result));
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto pre = get<T>();
    btas::scal(numeric_type(factor), pre);
    return eval_result<ResultTensorBTAS<T>>(std::move(pre));
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<annot_t>(ann[0]);
    auto const post_annot = std::any_cast<annot_t>(ann[1]);
    T result;
    btas::permute(get<T>(), pre_annot, result, post_annot);
    return eval_result<ResultTensorBTAS<T>>(std::move(result));
  }

  void add_inplace(Result const& other) override {
    auto& t = get<T>();
    auto const& o = other.get<T>();
    SEQUANT_ASSERT(t.range() == o.range());
    t += o;
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultTensorBTAS<T>>(column_symmetrize_btas(get<T>()));
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<ResultTensorBTAS<T>>(
        particle_antisymmetrize_btas(get<T>(), bra_rank));
  }

 private:
  [[nodiscard]] std::size_t size_in_bytes() const final {
    static_assert(std::is_arithmetic_v<typename T::value_type>);
    const auto& tensor = get<T>();
    // only count data
    return tensor.range().volume() * sizeof(T);
  }
};

}  // namespace sequant

#endif  // SEQUANT_HAS_BTAS

#endif  // SEQUANT_EVAL_BACKENDS_BTAS_RESULT_HPP
