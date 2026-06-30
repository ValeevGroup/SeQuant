#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/slotted_index.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/memoize.hpp>

#if defined(SEQUANT_HAS_TILEDARRAY)
#include <SeQuant/core/eval/backends/tiledarray/eval_expr.hpp>
#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#endif
#if defined(SEQUANT_HAS_BTAS)
#include <SeQuant/core/eval/backends/btas/eval_expr.hpp>
#include <SeQuant/core/eval/backends/btas/result.hpp>
#endif
#if defined(SEQUANT_HAS_TAPP)
#include <SeQuant/core/eval/backends/tapp/ops.hpp>
#include <SeQuant/core/eval/backends/tapp/tensor.hpp>
#endif

#include <range/v3/view/iota.hpp>

#include <concepts>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <optional>
#include <vector>

namespace sequant::mbpt {

static constexpr double default_biorthogonalizer_pseudoinverse_threshold =
    1e-12;

void biorthogonal_transform(
    ResultExpr& expr, double pseudoinverse_threshold =
                          default_biorthogonalizer_pseudoinverse_threshold);

void biorthogonal_transform(
    container::svector<ResultExpr>& exprs,
    double pseudoinverse_threshold =
        default_biorthogonalizer_pseudoinverse_threshold);

/// performs symbolic biorthogonal transform of CC-like equation using
///(for rank-3 and higher
/// [Wang-Knizia biorthogonalization](https://arxiv.org/abs/1805.00565).
///
/// @note uses hardcoded coefficients for ranks 1-5,
///  for higher ranks computes coefficients (if Eigen3 is available, else throws
///  an exception)
[[nodiscard]] ExprPtr biorthogonal_transform(
    const ExprPtr& expr,
    const container::svector<container::svector<sequant::SlottedIndex>>&
        ext_index_groups = {},
    double pseudoinverse_threshold =
        default_biorthogonalizer_pseudoinverse_threshold);
[[nodiscard]] ExprPtr biorthogonal_transform(
    const ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups = {},
    double pseudoinverse_threshold =
        default_biorthogonalizer_pseudoinverse_threshold);

/// @brief filters out the nonunique terms in Wang-Knizia biorthogonalization
/// WK biorthogonalization rewrites biorthogonal expressions as a projector
/// onto non-null-space (NNS)
/// applied to the biorthogonal expressions where out of each
/// group of terms related by permutation of external indices
/// those with the largest coefficients are selected.
/// This function performs the selection by forming groups of terms that
/// are equivalent modulo external index permutation (all terms in a group
/// have identical graph hashes).
/// @details This function processes a sum expression, grouping product terms by
/// hash of their canonicalized tensor network forms. For each group, it
/// retains only the terms with the largest absolute scalar coefficient.
/// @param expr The input expression, expected to be a `Sum` of `Product` terms.
/// @param ext_idxs A vector of external index groups. The function will not
/// apply the filtering logic if `ext_idxs.size()` is 2 or less.
/// @return A new `ExprPtr` representing the filtered and compacted expression.
ExprPtr WK_biorthogonalization_filter(
    ExprPtr expr,
    const container::svector<container::svector<SlottedIndex>>& ext_idxs);
ExprPtr WK_biorthogonalization_filter(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs);

/// @brief Four index layouts for triplet R2 metric primitives (id, bra, ket,
/// pair swaps on external index pairs).
struct TripletDoublesSwapLayouts {
  std::string orig;
  std::string bra_swap;
  std::string ket_swap;
  std::string pair_swap;
};

/// @brief Build bra/ket/pair swap layouts from a rank-4 comma-separated
/// annotation (a_1,a_2,i_1,i_2).
[[nodiscard]] TripletDoublesSwapLayouts triplet_doubles_swap_layouts(
    std::string const& orig_layout);

/// @brief Compact closed-shell triplet R2 by applying metric NNS weights
/// symbolically within each tensor-network hash group (id/bra/ket/pair swaps).
/// Produces one combined term per group (135 for 2h2p). Evaluation does not
/// need triplet_doubles_nns_project. Only applies when @p ext_idxs.size() == 2.
[[nodiscard]] ExprPtr triplet_doubles_symbolic_nns_compact(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs);

/// @brief Hash filter for closed-shell triplet R2: keep one term per tensor-
/// network hash, preferring the identity external-index layout (id swap) so
/// numerical triplet_doubles_nns_project can recover the full metric residual.
/// Falls back to largest |coefficient| within a group. Only applies when
/// @p ext_idxs.size() == 2.
[[nodiscard]] ExprPtr triplet_doubles_hash_filter(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs = {});

/// @brief Lossless compaction of the closed-shell triplet R2 residual: keep one
/// term per tensor-network hash, choosing the largest-|coefficient| member of
/// each group. Each group has the structure {c, c, c, -3c} (three external
/// Klein-four images with weight c and one odd-one-out with -3c), so the kept
/// term is always the -3c representative. This is the representative required
/// by triplet_doubles_symbolic_reconstruct to rebuild the dropped terms; the
/// layout-preferring triplet_doubles_hash_filter is NOT suitable for that since
/// it may keep a +c term. Only applies when @p ext_idxs.size() == 2.
[[nodiscard]] ExprPtr triplet_doubles_maxcoeff_compact(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs);

/// @brief Symbolic inverse of triplet_doubles_maxcoeff_compact: from each kept
/// -3c representative T (coefficient w = -3c), rebuild its external Klein-four
/// orbit by applying the bra-swap (i1<->i2), ket-swap (a1<->a2) and pair-swap
/// (both) to the external indices, each scaled by -1/3 (= c/w), and keeping T
/// itself. Reproduces the full 4-term group {w*T, c*bra(T), c*ket(T),
/// c*pair(T)}. Requires the -3c representative (use
/// triplet_doubles_maxcoeff_compact); reconstruction from a +c term is
/// ambiguous. Only applies when @p ext_idxs.size() == 2.
[[nodiscard]] ExprPtr triplet_doubles_symbolic_reconstruct(
    ExprPtr compact_expr,
    const container::svector<container::svector<Index>>& ext_idxs);

/// @brief Performs biorthogonal transformation with factored out NNS projector
/// @details Applies biorthogonal transformation. When factor_out_nns_projector
/// is true (default), factors out the NNS projector by applying additional
/// steps (S_maps and WK_biorthogonalization_filter) to produce compact
/// biorthogonal equations, necessitating a subsequent numerical NNS-projection
/// evaluation. When false, the NNS projector is not factored out, so no need to
/// apply numerical NNS-projection evaluation.
/// @param expr The input expression.
/// @param ext_idxs A vector of external index groups.
/// @param factor_out_nns_projector If true (default), factored out NNS
/// projector. If false, NNS projector is not factored out.
/// @return Expression pointer to the biorthogonalized result with leading S
/// operator.
ExprPtr biorthogonal_transform_pre_nnsproject(
    ExprPtr& expr,
    const container::svector<container::svector<SlottedIndex>>& ext_idxs,
    bool factor_out_nns_projector = true);
ExprPtr biorthogonal_transform_pre_nnsproject(
    ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_idxs,
    bool factor_out_nns_projector = true);

namespace detail {

/// \brief Computes the non-null space (NNS) projection coefficients
///
/// \param n_particles The rank of external index pairs
/// \param threshold The threshold to compute the pseudoinverse matrix
///        (set to default_biorth_threshold)
///
/// \return Vector of computed NNS projection coefficients
[[nodiscard]] std::vector<double> compute_nns_p_coeffs(
    std::size_t n_particles,
    double pseudoinverse_threshold =
        default_biorthogonalizer_pseudoinverse_threshold);

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
/// Finally, NNS projector = normalized_M_pinv · M.
///
/// \param n_particles The rank of external index pairs
///
/// \return Optional vector of NNS projector weights representing the last row,
///         std::nullopt if n_particles is outside the range [1,5].
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
///        (set to default_biorthogonalizer_pseudoinverse_threshold)
///
/// \return (memoized) Vector of hrdcoded/computed NNS projection weights
template <typename T>
  requires(std::floating_point<T> || meta::is_complex_v<T>)
[[nodiscard]] const std::vector<T>& nns_projection_weights(
    std::size_t n_particles,
    double pseudoinverse_threshold =
        default_biorthogonalizer_pseudoinverse_threshold) {
  static const std::vector<T> empty_vec{};

  if (n_particles < 3) {
    return empty_vec;
  }

  using CacheKey = std::pair<std::size_t, double>;

  static std::mutex cache_mutex;
  static std::condition_variable cache_cv;
  static container::map<CacheKey, std::optional<std::vector<T>>> cache;

  CacheKey key{n_particles, pseudoinverse_threshold};

  return sequant::detail::memoize(
      cache, cache_mutex, cache_cv, key, [&]() -> std::vector<T> {
        constexpr std::size_t max_rank_hardcoded_nns_projector = 5;
        if (n_particles <= max_rank_hardcoded_nns_projector) {
          if (auto hardcoded_coeffs = hardcoded_nns_projector<T>(n_particles)) {
            return std::move(hardcoded_coeffs.value());
          }
        }
        auto coeffs =
            detail::compute_nns_p_coeffs(n_particles, pseudoinverse_threshold);
        std::vector<T> nns_p_coeffs;
        nns_p_coeffs.reserve(coeffs.size());
        for (const auto& c : coeffs) {
          nns_p_coeffs.push_back(static_cast<T>(c));
        }
        return nns_p_coeffs;
      });
}

}  // namespace detail

#if defined(SEQUANT_HAS_TILEDARRAY)

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

  const auto& nns_p_coeffs =
      detail::nns_projection_weights<numeric_type>(ket_rank);

  TA::DistArray<Args...> result;

  sequant::detail::perm_t perm =
      iota(size_t{0}, rank) | ranges::to<sequant::detail::perm_t>;
  sequant::detail::perm_t bra_perm =
      iota(size_t{0}, bra_rank) | ranges::to<sequant::detail::perm_t>;
  sequant::detail::perm_t ket_perm =
      iota(bra_rank, rank) | ranges::to<sequant::detail::perm_t>;

  const auto lannot = sequant::detail::ords_to_annot(perm);

  if (ket_rank > 2 && !nns_p_coeffs.empty()) {
    const auto bra_annot =
        bra_rank == 0 ? "" : sequant::detail::ords_to_annot(bra_perm);

    size_t num_perms = nns_p_coeffs.size();
    for (size_t perm_rank = 0; perm_rank < num_perms; ++perm_rank) {
      sequant::detail::perm_t permuted_ket =
          detail::compute_permuted_indices(ket_perm, perm_rank, ket_rank);

      numeric_type coeff = nns_p_coeffs[perm_rank];

      const auto ket_annot = sequant::detail::ords_to_annot(permuted_ket);
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

template <typename... Args>
auto biorthogonal_nns_project(TA::DistArray<Args...> const& arr,
                              size_t bra_rank) {
  return biorthogonal_nns_project_ta(arr, bra_rank);
}

#endif  // defined(SEQUANT_HAS_TILEDARRAY)

#if defined(SEQUANT_HAS_BTAS)

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

  const auto& nns_p_coeffs =
      detail::nns_projection_weights<numeric_type>(ket_rank);

  btas::Tensor<Args...> result;

  sequant::detail::perm_t perm =
      iota(size_t{0}, rank) | ranges::to<sequant::detail::perm_t>;
  sequant::detail::perm_t bra_perm =
      iota(size_t{0}, bra_rank) | ranges::to<sequant::detail::perm_t>;
  sequant::detail::perm_t ket_perm =
      iota(bra_rank, rank) | ranges::to<sequant::detail::perm_t>;

  if (ket_rank > 2 && !nns_p_coeffs.empty()) {
    bool result_initialized = false;

    size_t num_perms = nns_p_coeffs.size();
    for (size_t perm_rank = 0; perm_rank < num_perms; ++perm_rank) {
      sequant::detail::perm_t permuted_ket =
          detail::compute_permuted_indices(ket_perm, perm_rank, ket_rank);

      numeric_type coeff = nns_p_coeffs[perm_rank];

      sequant::detail::perm_t annot = bra_perm;
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

template <typename... Args>
auto biorthogonal_nns_project(btas::Tensor<Args...> const& arr,
                              size_t bra_rank) {
  return biorthogonal_nns_project_btas(arr, bra_rank);
}

#endif  // defined(SEQUANT_HAS_BTAS)

#if defined(SEQUANT_HAS_TAPP)

/// \brief This function is used to implement
/// ResultPtr::biorthogonal_nns_project for TAPPTensor
///
/// \param arr The tensor to be "cleaned up"
/// \param bra_rank The rank of the bra indices
///
/// \return The cleaned TAPPTensor.
template <typename T, typename Alloc>
auto biorthogonal_nns_project_tapp(TAPPTensor<T, Alloc> const& arr,
                                   size_t bra_rank) {
  using ranges::views::iota;
  size_t const rank = arr.rank();
  SEQUANT_ASSERT(bra_rank <= rank);
  size_t const ket_rank = rank - bra_rank;

  // Residuals of rank 4 or less have no redundancy and don't require NNS
  // projection
  if (rank <= 4) return arr;

  using numeric_type = T;

  const auto& nns_p_coeffs =
      detail::nns_projection_weights<numeric_type>(ket_rank);

  using perm_type = container::svector<size_t>;

  TAPPTensor<T, Alloc> result;

  perm_type perm = iota(size_t{0}, rank) | ranges::to<perm_type>;
  perm_type bra_perm = iota(size_t{0}, bra_rank) | ranges::to<perm_type>;
  perm_type ket_perm = iota(bra_rank, rank) | ranges::to<perm_type>;

  if (ket_rank > 2 && !nns_p_coeffs.empty()) {
    bool result_initialized = false;

    size_t num_perms = nns_p_coeffs.size();
    for (size_t perm_rank = 0; perm_rank < num_perms; ++perm_rank) {
      perm_type permuted_ket =
          detail::compute_permuted_indices(ket_perm, perm_rank, ket_rank);

      numeric_type coeff = nns_p_coeffs[perm_rank];

      perm_type annot = bra_perm;
      annot.insert(annot.end(), permuted_ket.begin(), permuted_ket.end());

      container::svector<int64_t> annot_i64(annot.begin(), annot.end());
      container::svector<int64_t> perm_i64(perm.begin(), perm.end());

      TAPPTensor<T, Alloc> temp;
      tapp_ops::permute(arr, annot_i64, temp, perm_i64);
      tapp_ops::scal(coeff, temp);

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

template <typename T, typename Alloc>
auto biorthogonal_nns_project(TAPPTensor<T, Alloc> const& arr,
                              size_t bra_rank) {
  return biorthogonal_nns_project_tapp(arr, bra_rank);
}

#endif  // defined(SEQUANT_HAS_TAPP)

#if defined(SEQUANT_HAS_TILEDARRAY)

namespace detail {

/// @brief Weighted average over the triplet-doubles Klein-four orbit:
///   out(orig) = w_self * arr(orig)
///             + w_swap * (arr(bra_swap) + arr(ket_swap) + arr(pair_swap)).
/// All three external swap images share a single weight by the Klein-four
/// symmetry, so the metric NNS reconstruction (w_self = 1, w_swap = -1/3) and
/// the idempotent null-space projector (w_self = 3/4, w_swap = -1/4) are both
/// instances of this combine; the two differ only by the overall scale 4/3.
template <typename... Args>
auto triplet_doubles_orbit_combine_ta(
    TA::DistArray<Args...> const& arr, TripletDoublesSwapLayouts const& layouts,
    typename TA::DistArray<Args...>::numeric_type w_self,
    typename TA::DistArray<Args...>::numeric_type w_swap) {
  TA::DistArray<Args...> result;
  result(layouts.orig) =
      w_self * arr(layouts.orig) +
      w_swap * (arr(layouts.bra_swap) + arr(layouts.ket_swap) +
                arr(layouts.pair_swap));
  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());
  result.truncate();
  return result;
}

}  // namespace detail

/// @brief Metric NNS reconstruction for compact triplet R2 residuals: rebuilds
/// the full residual from its compact (-3c) representative via the Klein-four
/// orbit with the n=2 triplet Gram pseudo-inverse weights {1, -1/3, -1/3,
/// -1/3}, i.e. out(orig) = arr(orig) - (1/3)(bra + ket + pair). Numerical
/// analogue of triplet_doubles_symbolic_reconstruct; apply to the H*R residual
/// when the compact equations were evaluated.
template <typename... Args>
auto triplet_doubles_nns_project_ta(TA::DistArray<Args...> const& arr,
                                    std::string const& orig_layout) {
  using numeric_type = typename TA::DistArray<Args...>::numeric_type;
  if (arr.trange().rank() != 4) return arr;
  const auto layouts = triplet_doubles_swap_layouts(orig_layout);
  return detail::triplet_doubles_orbit_combine_ta<Args...>(
      arr, layouts, numeric_type(1), -numeric_type(1) / numeric_type(3));
}

template <typename... Args>
auto triplet_doubles_nns_project(TA::DistArray<Args...> const& arr,
                                 std::string const& orig_layout) {
  return triplet_doubles_nns_project_ta(arr, orig_layout);
}

/// @brief Idempotent null-space projector for triplet R2: removes the
/// metric-null fully symmetric component, P = I - (1/4)Sigma, i.e. the
/// Klein-four orbit with weights {3/4, -1/4, -1/4, -1/4}. Equals
/// (3/4) * triplet_doubles_nns_project. Apply to the Davidson trial R2 each
/// iteration to keep it in the physical (non-null) subspace.
template <typename... Args>
auto triplet_doubles_nullspace_project_ta(TA::DistArray<Args...> const& arr,
                                          std::string const& orig_layout) {
  using numeric_type = typename TA::DistArray<Args...>::numeric_type;
  if (arr.trange().rank() != 4) return arr;
  const auto layouts = triplet_doubles_swap_layouts(orig_layout);
  return detail::triplet_doubles_orbit_combine_ta<Args...>(
      arr, layouts, numeric_type(3) / numeric_type(4),
      -numeric_type(1) / numeric_type(4));
}

template <typename... Args>
auto triplet_doubles_nullspace_project(TA::DistArray<Args...> const& arr,
                                       std::string const& orig_layout) {
  return triplet_doubles_nullspace_project_ta(arr, orig_layout);
}

/// @brief TE-only -> full-metric reconstruction for triplet R2 (PI conjecture
/// experiment). The bare-TE residual te_a = TE/4 sits in the same Klein-four
/// orbits as the production residual Omega = (3TE - TE_ps)/16, and the exact
/// identity Omega = te_a + (1/4)(bra_swap(te_a) + ket_swap(te_a)) recovers the
/// production residual via a partial (bra + ket, NO pair) Klein-four
/// symmetrization. Apply to the H*R residual when the te_only equations were
/// evaluated (Knob A: TE_ps dropped, R+R_swap amplitude kept).
template <typename... Args>
auto triplet_doubles_te_reconstruct_ta(TA::DistArray<Args...> const& arr,
                                       std::string const& orig_layout) {
  using numeric_type = typename TA::DistArray<Args...>::numeric_type;
  if (arr.trange().rank() != 4) return arr;
  const auto layouts = triplet_doubles_swap_layouts(orig_layout);
  TA::DistArray<Args...> result;
  result(layouts.orig) =
      arr(layouts.orig) + (numeric_type(1) / numeric_type(4)) *
                              (arr(layouts.bra_swap) + arr(layouts.ket_swap));
  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());
  result.truncate();
  return result;
}

template <typename... Args>
auto triplet_doubles_te_reconstruct(TA::DistArray<Args...> const& arr,
                                    std::string const& orig_layout) {
  return triplet_doubles_te_reconstruct_ta(arr, orig_layout);
}

#endif  // defined(SEQUANT_HAS_TILEDARRAY)

#if defined(SEQUANT_HAS_BTAS)

namespace detail {

/// @brief BTAS analogue of triplet_doubles_orbit_combine_ta (no truncation).
template <typename... Args>
auto triplet_doubles_orbit_combine_btas(
    btas::Tensor<Args...> const& arr, TripletDoublesSwapLayouts const& layouts,
    typename btas::Tensor<Args...>::value_type w_self,
    typename btas::Tensor<Args...>::value_type w_swap) {
  btas::Tensor<Args...> result;
  result(layouts.orig) =
      w_self * arr(layouts.orig) +
      w_swap * (arr(layouts.bra_swap) + arr(layouts.ket_swap) +
                arr(layouts.pair_swap));
  return result;
}

}  // namespace detail

template <typename... Args>
auto triplet_doubles_nns_project_btas(btas::Tensor<Args...> const& arr,
                                      std::string const& orig_layout) {
  using numeric_type = typename btas::Tensor<Args...>::value_type;
  if (arr.rank() != 4) return arr;
  const auto layouts = triplet_doubles_swap_layouts(orig_layout);
  return detail::triplet_doubles_orbit_combine_btas<Args...>(
      arr, layouts, numeric_type(1), -numeric_type(1) / numeric_type(3));
}

template <typename... Args>
auto triplet_doubles_nns_project(btas::Tensor<Args...> const& arr,
                                 std::string const& orig_layout) {
  return triplet_doubles_nns_project_btas(arr, orig_layout);
}

template <typename... Args>
auto triplet_doubles_nullspace_project_btas(btas::Tensor<Args...> const& arr,
                                            std::string const& orig_layout) {
  using numeric_type = typename btas::Tensor<Args...>::value_type;
  if (arr.rank() != 4) return arr;
  const auto layouts = triplet_doubles_swap_layouts(orig_layout);
  return detail::triplet_doubles_orbit_combine_btas<Args...>(
      arr, layouts, numeric_type(3) / numeric_type(4),
      -numeric_type(1) / numeric_type(4));
}

template <typename... Args>
auto triplet_doubles_nullspace_project(btas::Tensor<Args...> const& arr,
                                       std::string const& orig_layout) {
  return triplet_doubles_nullspace_project_btas(arr, orig_layout);
}

#endif  // defined(SEQUANT_HAS_BTAS)

}  // namespace sequant::mbpt

#endif
