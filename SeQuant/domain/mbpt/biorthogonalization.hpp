#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization_hardcoded.hpp>

#include <condition_variable>

namespace sequant {

namespace {
static constexpr double default_biorth_threshold = 1e-12;
}  // namespace

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
}  // namespace sequant

#endif
