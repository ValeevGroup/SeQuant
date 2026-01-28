#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/math.hpp>

#include <TiledArray/einsum/tiledarray.h>
#include <tiledarray.h>

namespace sequant {

namespace {

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
auto column_symmetrize_ta(TA::DistArray<Args...> const& arr) {
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

  auto const nf = static_cast<double>(rational{1, factorial(nparticles)});
  result(lannot) = nf * result(lannot);

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
  SEQUANT_ASSERT(bra_rank <= rank);
  size_t const ket_rank = rank - bra_rank;

  if (bra_rank <= 1 && ket_rank <= 1) {
    // nothing to do
    return arr;
  }

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;
  perm_t bra_perm = iota(size_t{0}, bra_rank) | ranges::to<perm_t>;
  perm_t ket_perm = iota(bra_rank, rank) | ranges::to<perm_t>;

  const auto lannot = ords_to_annot(perm);

  auto process_permutations = [&lannot](const TA::DistArray<Args...>& input_arr,
                                        size_t range_rank, perm_t range_perm,
                                        const std::string& other_annot,
                                        bool is_bra) -> TA::DistArray<Args...> {
    if (range_rank <= 1) return input_arr;
    TA::DistArray<Args...> result;

    auto callback = [&](int parity) {
      const auto range_annot = ords_to_annot(range_perm);
      const auto annot = other_annot.empty()
                             ? range_annot
                             : (is_bra ? range_annot + "," + other_annot
                                       : other_annot + "," + range_annot);

      typename decltype(result)::numeric_type p_ = parity == 0 ? 1 : -1;
      if (result.is_initialized()) {
        result(lannot) += p_ * input_arr(annot);
      } else {
        result(lannot) = p_ * input_arr(annot);
      }
    };
    antisymmetric_permutation(ParticleRange{range_perm.begin(), range_rank},
                              callback);
    return result;
  };

  // Process bra permutations first
  const auto ket_annot = ket_rank == 0 ? "" : ords_to_annot(ket_perm);
  auto result = process_permutations(arr, bra_rank, bra_perm, ket_annot, true);

  // Process ket permutations
  const auto bra_annot = bra_rank == 0 ? "" : ords_to_annot(bra_perm);
  result = process_permutations(result, ket_rank, ket_perm, bra_annot, false);

  auto const nf = static_cast<double>(
      rational{1, factorial(bra_rank) * factorial(ket_rank)});
  result(lannot) = nf * result(lannot);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());
  return result;
}

template <typename... Args>
inline void log_ta(Args const&... args) noexcept {
  log_result("[TA] ", args...);
}

/// Convert sequant::DeNest to TA::DeNest
inline constexpr TA::DeNest to_ta_denest(DeNest d) noexcept {
  return d == DeNest::True ? TA::DeNest::True : TA::DeNest::False;
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

///
/// \brief Result for a tensor value of TA::DistArray type.
/// \tparam ArrayT TA::DistArray type. Tile type of ArrayT is regular tensor of
///                scalars (not a tensor of tensors)
///
template <typename ArrayT, typename = std::enable_if_t<TA::detail::is_tensor_v<
                               typename ArrayT::value_type>>>
class ResultTensorTA final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename ArrayT::numeric_type;

  explicit ResultTensorTA(ArrayT arr) : Result{std::move(arr)} {}

 private:
  using this_type = ResultTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const a = annot_wrap{annot};

    log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    result(a.this_annot) =
        get<ArrayT>()(a.lannot) + other.get<ArrayT>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      auto result = get<ArrayT>();
      auto scalar = other.get<numeric_type>();

      log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

      result(a.this_annot) = scalar * result(a.lannot);

      decltype(result)::wait_for_lazy_cleanup(result.world());
      return eval_result<this_type>(std::move(result));
    }

    if (a.this_annot.empty()) {
      // DOT product
      SEQUANT_ASSERT(other.is<this_type>());
      numeric_type d =
          TA::dot(get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot));
      ArrayT::wait_for_lazy_cleanup(get<ArrayT>().world());
      ArrayT::wait_for_lazy_cleanup(other.get<ArrayT>().world());

      log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      return eval_result<ResultScalar<numeric_type>>(d);
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

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto pre = get<ArrayT>();
    TA::scale(pre, numeric_type(factor));
    return eval_result<this_type>(std::move(pre));
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    log_ta(pre_annot, " = ", post_annot, "\n");

    ArrayT result;
    result(post_annot) = get<ArrayT>()(pre_annot);
    ArrayT::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<this_type>());

    auto& t = get<ArrayT>();
    auto const& o = other.get<ArrayT>();

    SEQUANT_ASSERT(t.trange() == o.trange());
    auto ann = TA::detail::dummy_annotation(t.trange().rank());

    log_ta(ann, " += ", ann, "\n");

    t(ann) += o(ann);
    ArrayT::wait_for_lazy_cleanup(t.world());
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<this_type>(column_symmetrize_ta(get<ArrayT>()));
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t bra_rank) const override {
    return eval_result<this_type>(
        particle_antisymmetrize_ta(get<ArrayT>(), bra_rank));
  }

 private:
  [[nodiscard]] std::size_t size_in_bytes() const final {
    auto& v = get<ArrayT>();
    auto local_size = TA::size_of<TA::MemorySpace::Host>(v);
    v.world().gop.sum(local_size);
    return local_size;
  }
};

template <typename ArrayT,
          typename = std::enable_if_t<
              TA::detail::is_tensor_of_tensor_v<typename ArrayT::value_type>>>
class ResultTensorOfTensorTA final : public Result {
 public:
  using Result::id_t;
  using numeric_type = typename ArrayT::numeric_type;

  explicit ResultTensorOfTensorTA(ArrayT arr) : Result{std::move(arr)} {}

 private:
  using this_type = ResultTensorOfTensorTA<ArrayT>;
  using annot_wrap = Annot<std::string>;

  using _inner_tensor_type = typename ArrayT::value_type::value_type;

  using compatible_regular_distarray_type =
      TA::DistArray<_inner_tensor_type, typename ArrayT::policy_type>;

  // Only @c that_type type is allowed for ToT * T computation
  using that_type = ResultTensorTA<compatible_regular_distarray_type>;

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<this_type>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    SEQUANT_ASSERT(other.is<this_type>());
    auto const a = annot_wrap{annot};

    log_ta(a.lannot, " + ", a.rannot, " = ", a.this_annot, "\n");

    ArrayT result;
    result(a.this_annot) =
        get<ArrayT>()(a.lannot) + other.get<ArrayT>()(a.rannot);
    decltype(result)::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest DeNestFlag) const override {
    auto const a = annot_wrap{annot};

    if (other.is<ResultScalar<numeric_type>>()) {
      auto result = get<ArrayT>();
      auto scalar = other.get<numeric_type>();

      log_ta(a.lannot, " * ", scalar, " = ", a.this_annot, "\n");

      result(a.this_annot) = scalar * result(a.lannot);

      decltype(result)::wait_for_lazy_cleanup(result.world());
      return eval_result<this_type>(std::move(result));
    } else if (a.this_annot.empty()) {
      // DOT product
      SEQUANT_ASSERT(other.is<this_type>());
      numeric_type d =
          TA::dot(get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot));
      ArrayT::wait_for_lazy_cleanup(get<ArrayT>().world());
      ArrayT::wait_for_lazy_cleanup(other.get<ArrayT>().world());

      log_ta(a.lannot, " * ", a.rannot, " = ", d, "\n");

      return eval_result<ResultScalar<numeric_type>>(d);
    }

    log_ta(a.lannot, " * ", a.rannot, " = ", a.this_annot, "\n");

    if (other.is<that_type>()) {
      // ToT * T -> ToT
      auto result =
          TA::einsum(get<ArrayT>()(a.lannot),
                     other.get<compatible_regular_distarray_type>()(a.rannot),
                     a.this_annot);
      return eval_result<this_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == DeNest::True) {
      // ToT * ToT -> T
      auto result = TA::einsum<TA::DeNest::True>(
          get<ArrayT>()(a.lannot), other.get<ArrayT>()(a.rannot), a.this_annot);
      return eval_result<that_type>(std::move(result));

    } else if (other.is<this_type>() && DeNestFlag == DeNest::False) {
      // ToT * ToT -> ToT
      auto result = TA::einsum(get<ArrayT>()(a.lannot),
                               other.get<ArrayT>()(a.rannot), a.this_annot);
      return eval_result<this_type>(std::move(result));
    } else {
      throw invalid_operand();
    }
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t factor) const override {
    auto pre = get<ArrayT>();
    TA::scale(pre, numeric_type(factor));
    return eval_result<this_type>(std::move(pre));
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    auto const pre_annot = std::any_cast<std::string>(ann[0]);
    auto const post_annot = std::any_cast<std::string>(ann[1]);

    log_ta(pre_annot, " = ", post_annot, "\n");

    ArrayT result;
    result(post_annot) = get<ArrayT>()(pre_annot);
    ArrayT::wait_for_lazy_cleanup(result.world());
    return eval_result<this_type>(std::move(result));
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<this_type>());

    auto& t = get<ArrayT>();
    auto const& o = other.get<ArrayT>();

    SEQUANT_ASSERT(t.trange() == o.trange());
    auto ann = TA::detail::dummy_annotation(t.trange().rank());

    log_ta(ann, " += ", ann, "\n");

    t(ann) += o(ann);
    ArrayT::wait_for_lazy_cleanup(t.world());
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    // not implemented yet
    return nullptr;
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    // not implemented yet
    return nullptr;
  }

 private:
  [[nodiscard]] std::size_t size_in_bytes() const final {
    auto& v = get<ArrayT>();
    auto local_size = TA::size_of<TA::MemorySpace::Host>(v);
    v.world().gop.sum(local_size);
    return local_size;
  }
};

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_RESULT_HPP
