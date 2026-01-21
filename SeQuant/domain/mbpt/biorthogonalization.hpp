#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#if defined(SEQUANT_HAS_TILEDARRAY) || defined(SEQUANT_HAS_BTAS)
#include <SeQuant/core/eval/eval.hpp>
#endif

#include <condition_variable>

namespace sequant {

namespace {
static constexpr double default_biorth_threshold = 1e-12;
}  // namespace

// clang-format off
/// \brief Provides the first row of the biorthogonal coefficients matrix,
/// hardcoded from Mathematica to avoid numerical precision loss.
///
/// The Myrvold-Ruskey unrank1 algorithm (doi.org/10.1016/S0020-0190(01)00141-7)
/// is used to order permutations, then the permutational overlap matrix M is
/// constructed with elements (-2)^{c} × (-1)^{n_particles}, where c is the
/// number of cycles in the relative permutation.
///
/// The biorthogonal coefficients are obtained from the normalized pseudoinverse
/// of M: first compute M_pinv (the pseudoinverse), then normalize it by the
/// factor ((n_particles)!/rank(M)).
/// Finally, biorthogonal coefficients = normalized_M_pinv · e_1,
/// where e_1 is the first unit vector.
/// See [DOI 10.48550/ARXIV.1805.00565](https://doi.org/10.48550/ARXIV.1805.00565)
/// for more details.
///
/// \param n_particles The rank of external index pairs
///
/// \return Vector of rational coefficients representing the first row
///
/// \throw std::runtime_error if n_particles is not in the range [1,5]
// clang-format on
std::vector<sequant::rational> hardcoded_biorth_coeffs_first_row(
    std::size_t n_particles);

[[nodiscard]] ResultExpr biorthogonal_transform_copy(
    const ResultExpr& expr, double threshold = default_biorth_threshold);

[[nodiscard]] container::svector<ResultExpr> biorthogonal_transform_copy(
    const container::svector<ResultExpr>& exprs,
    double threshold = default_biorth_threshold);

void biorthogonal_transform(ResultExpr& expr,
                            double threshold = default_biorth_threshold);

void biorthogonal_transform(container::svector<ResultExpr>& exprs,
                            double threshold = default_biorth_threshold);

/// performs symbolic biorthogonal transform of CC-like equation using
///(for rank-3 and higher
/// Wang-Knizia biorthogonalization (https://arxiv.org/abs/1805.00565) is used
/// uses hardcoded coefficients for ranks 1-6, computed coefficients for higher
/// ranks
[[nodiscard]] ExprPtr biorthogonal_transform(
    const ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups = {},
    double threshold = default_biorth_threshold);

/// \brief Computes the non-null space (NNS) projection coefficients
///
/// \param n_particles The rank of external index pairs
/// \param threshold The threshold to compute the pseudoinverse matrix
///        (set to default_biorth_threshold)
///
/// \return Vector of computed NNS projection coefficients
[[nodiscard]] std::vector<double> compute_nns_p_coeffs(
    std::size_t n_particles, double threshold = default_biorth_threshold);

/// \brief Provides permuted indices using libperm unrank function
///
/// \param indices The indices to permute
/// \param perm_rank The rank of the permutation
/// \param n_particles The rank of external index pairs
///
/// \return The permuted indices
container::svector<size_t> compute_permuted_indices(
    const container::svector<size_t>& indices, size_t perm_rank,
    size_t n_particles);

/// \brief Provides one row of the NNS projector matrix,
/// hardcoded from Mathematica to avoid numerical precision loss.
///
/// The NNS projector weights are obtained from the normalized pseudoinverse
/// of M: first compute M_pinv (the pseudoinverse), then normalize it by the
/// factor ((n_particles)!/rank(M)).
//// Finally, NNS projector = normalized_M_pinv · M.
///
/// \param n_particles The rank of external index pairs
///
/// \return Vector of NNS projector weights representing the last row
///
/// \throw std::runtime_error if n_particles is not in the range [1,5]
template <typename T>
  requires(std::floating_point<T> || meta::is_complex_v<T>)
std::optional<std::vector<T>> hardcoded_nns_projector(std::size_t n_particles) {
  switch (n_particles) {
    case 1:
      return std::vector<T>{T(1) / T(1)};

    case 2:
      return std::vector<T>{T(0) / T(1), T(1) / T(1)};

    case 3:
      return std::vector<T>{T(-1) / T(5), T(-1) / T(5), T(-1) / T(5),
                            T(-1) / T(5), T(-1) / T(5), T(1) / T(1)};

    case 4:
      return std::vector<T>{
          T(1) / T(7),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(14),
          T(1) / T(7),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(1) / T(7),   T(-2) / T(7),
          T(-1) / T(14), T(1) / T(7),   T(-1) / T(14), T(-2) / T(7),
          T(1) / T(7),   T(-1) / T(14), T(-1) / T(14), T(-2) / T(7),
          T(-2) / T(7),  T(-2) / T(7),  T(-2) / T(7),  T(1) / T(1)};

    case 5:
      return std::vector<T>{
          T(-1) / T(14), T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(2) / T(21),  T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(2) / T(21),  T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(0) / T(1),
          T(-1) / T(14), T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(2) / T(21),  T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(2) / T(21),  T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(0) / T(1),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(2) / T(21),
          T(0) / T(1),   T(2) / T(21),  T(2) / T(21),  T(-1) / T(21),
          T(2) / T(21),  T(0) / T(1),   T(-1) / T(21), T(-1) / T(21),
          T(-1) / T(21), T(-1) / T(21), T(1) / T(7),   T(0) / T(1),
          T(0) / T(1),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(3),
          T(2) / T(21),  T(-1) / T(21), T(2) / T(21),  T(2) / T(21),
          T(0) / T(1),   T(-1) / T(21), T(-1) / T(21), T(-1) / T(21),
          T(-1) / T(21), T(1) / T(7),   T(2) / T(21),  T(-1) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(0) / T(1),   T(0) / T(1),
          T(1) / T(7),   T(0) / T(1),   T(1) / T(7),   T(-1) / T(3),
          T(-1) / T(21), T(-1) / T(21), T(-1) / T(21), T(-1) / T(21),
          T(1) / T(7),   T(-1) / T(21), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(0) / T(1),   T(-1) / T(21), T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(0) / T(1),   T(1) / T(7),
          T(0) / T(1),   T(0) / T(1),   T(1) / T(7),   T(-1) / T(3),
          T(0) / T(1),   T(1) / T(7),   T(1) / T(7),   T(0) / T(1),
          T(-1) / T(3),  T(1) / T(7),   T(0) / T(1),   T(1) / T(7),
          T(0) / T(1),   T(-1) / T(3),  T(1) / T(7),   T(1) / T(7),
          T(0) / T(1),   T(0) / T(1),   T(-1) / T(3),  T(-1) / T(3),
          T(-1) / T(3),  T(-1) / T(3),  T(-1) / T(3),  T(1) / T(1)};

    default:
      return std::nullopt;
  }
}

/// \brief Provides NNS projection weights for a given rank
///
/// \tparam T The numeric type (must be floating point or complex)
/// \param n_particles The rank of external index pairs
/// \param threshold The threshold to compute the pseudoinverse matrix
///        (set to default_biorth_threshold)
///
/// \return (memoized) Vector of hrdcoded/computed NNS projection weights
template <typename T>
  requires(std::floating_point<T> || meta::is_complex_v<T>)
[[nodiscard]] const std::vector<T>& nns_projection_weights(
    std::size_t n_particles, double threshold = default_biorth_threshold) {
  static const std::vector<T> empty_vec{};

  if (n_particles < 3) {
    return empty_vec;
  }

  using CacheKey = std::pair<std::size_t, double>;
  using CacheValue = std::optional<std::vector<T>>;

  static std::mutex cache_mutex;
  static std::condition_variable cache_cv;
  static container::map<CacheKey, CacheValue> cache;

  CacheKey key{n_particles, threshold};

  {
    std::unique_lock<std::mutex> lock(cache_mutex);
    auto [it, inserted] = cache.try_emplace(key, std::nullopt);
    if (!inserted) {
      cache_cv.wait(lock, [&] { return it->second.has_value(); });
      return it->second.value();
    }
  }

  std::vector<T> nns_p_coeffs;

  constexpr std::size_t max_rank_hardcoded_nns_projector = 5;
  if (n_particles <= max_rank_hardcoded_nns_projector) {
    auto hardcoded_coeffs = hardcoded_nns_projector<T>(n_particles);
    if (hardcoded_coeffs) {
      nns_p_coeffs = std::move(hardcoded_coeffs.value());
    }
  } else {
    auto coeffs = compute_nns_p_coeffs(n_particles, threshold);
    nns_p_coeffs.reserve(coeffs.size());
    for (const auto& c : coeffs) {
      nns_p_coeffs.push_back(static_cast<T>(c));
    }
  }

  {
    std::lock_guard<std::mutex> lock(cache_mutex);
    cache[key] = std::move(nns_p_coeffs);
    cache_cv.notify_all();
    return cache[key].value();
  }
}

#if defined(SEQUANT_HAS_TILEDARRAY) || defined(SEQUANT_HAS_BTAS)

/// \brief This function is used to implement
/// ResultPtr::biorthogonal_nns_project for TA::DistArray
///
/// \param arr The array to be "cleaned up"
/// \param bra_rank The rank of the bra indices
///
/// \return The cleaned TA::DistArray.
template <typename... Args>
auto biorthogonal_nns_project_ta(TA::DistArray<Args...> const& arr,
                                 size_t bra_rank) {
  using ranges::views::iota;
  size_t const rank = arr.trange().rank();
  SEQUANT_ASSERT(bra_rank <= rank);
  size_t const ket_rank = rank - bra_rank;

  // Residuals of rank 4 or less have no redundancy and don't require NNS
  // projection
  if (rank <= 4) return arr;

  using numeric_type = typename TA::DistArray<Args...>::numeric_type;

  const auto& nns_p_coeffs = nns_projection_weights<numeric_type>(ket_rank);

  TA::DistArray<Args...> result;

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;
  perm_t bra_perm = iota(size_t{0}, bra_rank) | ranges::to<perm_t>;
  perm_t ket_perm = iota(bra_rank, rank) | ranges::to<perm_t>;

  const auto lannot = ords_to_annot(perm);

  if (ket_rank > 2 && !nns_p_coeffs.empty()) {
    const auto bra_annot = bra_rank == 0 ? "" : ords_to_annot(bra_perm);

    size_t num_perms = nns_p_coeffs.size();
    for (size_t perm_rank = 0; perm_rank < num_perms; ++perm_rank) {
      perm_t permuted_ket =
          compute_permuted_indices(ket_perm, perm_rank, ket_rank);

      numeric_type coeff = nns_p_coeffs[perm_rank];

      const auto ket_annot = ords_to_annot(permuted_ket);
      const auto annot =
          bra_annot.empty() ? ket_annot : bra_annot + "," + ket_annot;

      if (result.is_initialized()) {
        result(lannot) += coeff * arr(annot);
      } else {
        result(lannot) = coeff * arr(annot);
      }
    }
  } else {
    result(lannot) = arr(lannot);
  }

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());
  return result;
}

/// \brief This function is used to implement
/// ResultPtr::biorthogonal_nns_project for btas::Tensor
///
/// \param arr The array to be "cleaned up"
/// \param bra_rank The rank of the bra indices
///
/// \return The cleaned btas::Tensor.
template <typename... Args>
auto biorthogonal_nns_project_btas(btas::Tensor<Args...> const& arr,
                                   size_t bra_rank) {
  using ranges::views::iota;
  size_t const rank = arr.rank();
  SEQUANT_ASSERT(bra_rank <= rank);
  size_t const ket_rank = rank - bra_rank;

  // Residuals of rank 4 or less have no redundancy and don't require NNS
  // projection
  if (rank <= 4) return arr;

  using numeric_type = typename btas::Tensor<Args...>::numeric_type;

  const auto& nns_p_coeffs = nns_projection_weights<numeric_type>(ket_rank);

  btas::Tensor<Args...> result;

  perm_t perm = iota(size_t{0}, rank) | ranges::to<perm_t>;
  perm_t bra_perm = iota(size_t{0}, bra_rank) | ranges::to<perm_t>;
  perm_t ket_perm = iota(bra_rank, rank) | ranges::to<perm_t>;

  if (ket_rank > 2 && !nns_p_coeffs.empty()) {
    bool result_initialized = false;

    size_t num_perms = nns_p_coeffs.size();
    for (size_t perm_rank = 0; perm_rank < num_perms; ++perm_rank) {
      perm_t permuted_ket =
          compute_permuted_indices(ket_perm, perm_rank, ket_rank);

      numeric_type coeff = nns_p_coeffs[perm_rank];

      perm_t annot = bra_perm;
      annot.insert(annot.end(), permuted_ket.begin(), permuted_ket.end());

      btas::Tensor<Args...> temp;
      btas::permute(arr, annot, temp, perm);
      btas::scal(coeff, temp);

      if (result_initialized) {
        result += temp;
      } else {
        result = temp;
        result_initialized = true;
      }
    }

  } else {
    result = arr;
  }

  return result;
}

#endif

}  // namespace sequant

#endif
