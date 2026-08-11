#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/rational.hpp>
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

/// TODO separate them to two categories:  biorthogonalization and NNS
enum class TripletWeightKind {
  /// idempotent null-space projector row
  NullspaceProjector,
  /// identity-normalized NullspaceProjector row (metric NNS reconstruction)
  NnsReconstruction,
  /// bare-TE (te_only) undo-compact row {1, -1/2, -1/2, 0}; n = 2 only
  TeNnsReconstruction,
  /// bare-TE -> paper-metric reconstruction row {1, 1/4, 1/4, 0}; n = 2 only
  TeReconstruction,
  /// paper-combined residual assembly weights (see triplet_combined_residual)
  CombinedResidual,
  /// bare-TE (te_only) residual assembly weights {1/4, 0, 0, 0}; n = 2 only
  TeCombinedResidual
};

/// \brief Compaction of the closed-shell triplet residual: keeps one
/// representative slot permutation per tensor-network hash group.
/// \param expr The residual; returned unchanged unless it is a Sum
/// \param ext_idxs A vector of external index groups
/// \param kind The weight row that generated the residual
/// \return The compacted expression
/// \throw Exception if a hash group does not match the single-representative
[[nodiscard]] ExprPtr triplet_maxcoeff_compact(
    ExprPtr expr, const container::svector<container::svector<Index>>& ext_idxs,
    TripletWeightKind kind = TripletWeightKind::NnsReconstruction);

// clang-format off
/// \brief Assembles the closed-shell triplet residual from the
/// sector-summed primitive V, for any supported rank
/// (rank 2 based on Kohn's paper)
///
///   n = 2: Omega = (3 V - V_ps)/16 (Faber's paper)
///   te_only (n = 2, EFV experiment): Omega = V/4, i.e. the bare TE primitive
///          with the external pair swap dropped. The dropped part is restored
///          by the postprocessing Omega = V/4 + (V_bs + V_ks)/16
///   n = 3: Omega = (6 V - V_ps01 - V_ps02 + 2 V_ks12)/160
///          the 18 TEE ops (6 pairings x 3 T positions) have rank 9;
///
/// Supporting a new rank means adding weights to the CombinedResidual case
/// of the triplet weight table
///
/// \param V The sector-summed triplet primitive
/// \param ext_idxs A vector of external index groups (must have 2 or 3 groups)
/// \param te_only Assemble the bare-TE weights instead (n = 2 only)
/// \note I might want to have te_only for triples if EFV prefer it (tee_only)
/// \return The combined residual, simplified
/// \throw Exception for unsupported \p ext_idxs sizes
// clang-format on
[[nodiscard]] ExprPtr triplet_combined_residual(
    const ExprPtr& V,
    const container::svector<container::svector<Index>>& ext_idxs,
    bool te_only = false);

/// \brief Symbolic inverse of triplet_maxcoeff_compact: rebuilds the full
/// residual as the weighted sum of the kept terms over the (n!)^2 external
/// slot permutations
///
/// \param compact_expr The compact residual; returned unchanged unless a Sum
/// \param ext_idxs A vector of external index groups (must have 2 or 3 groups)
/// \param kind The weight row to apply
/// \return The reconstructed expression
[[nodiscard]] ExprPtr triplet_symbolic_reconstruct(
    ExprPtr compact_expr,
    const container::svector<container::svector<Index>>& ext_idxs,
    TripletWeightKind kind = TripletWeightKind::NnsReconstruction);

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

/// \brief Provides bra and ket permuted indices
/// for the S_n x S_n external slot permutations, by applying
/// compute_permuted_indices to the bra and the ket slots separately
///
/// \param perm_index The flat perm index, in [0, (n!)^2)
/// \param n_particles The rank of external index pairs
///
/// \return The 2*n_particles slot orders: slot s displays slot ords[s], the
///         first n_particles slots being the bra and the rest the ket ones
[[nodiscard]] container::svector<size_t> compute_bra_ket_permuted_indices(
    size_t perm_index, size_t n_particles);

/// \brief Provides the annotations of the S_n x S_n
/// external slot permutations of a rank-2n annotation
///
/// \param orig_annot A comma-separated rank-2n annotation, e.g.
///        "a_1,a_2,i_1,i_2" (bra slots first, then ket slots)
/// \param n_particles The rank of external index pairs
///
/// \return The (n!)^2 permuted annotations in flat perm order
[[nodiscard]] container::svector<std::string> slot_perm_annots(
    std::string const& orig_annot, size_t n_particles);

/// \brief Computes the weight row over the S_n x S_n
/// external slot permutations for the closed-shell triplet primitives
///
/// \param n_particles The rank of external index pairs (2 or 3; 2 only for
///        the bare-TE kinds)
/// \param kind The weight row to compute
///
/// \return Vector of weights in flat perm order
///
/// \throw Exception for unsupported \p n_particles / \p kind combinations
[[nodiscard]] std::vector<double> compute_triplet_weights(
    std::size_t n_particles, TripletWeightKind kind);

/// \brief Provides the numeric weight row over the S_n x S_n
/// external slot permutations for the closed-shell triplet primitives
///
/// \tparam T The numeric type (must be floating point or complex)
/// \param n_particles The rank of external index pairs
/// \param kind The weight row to provide
///
/// \return (memoized) Vector of weights in flat perm order; empty unless
///         \p n_particles is 2 or 3 (2 only for the bare-TE kinds)
template <typename T>
  requires(std::floating_point<T> || meta::is_complex_v<T>)
[[nodiscard]] const std::vector<T>& triplet_weights(std::size_t n_particles,
                                                    TripletWeightKind kind) {
  static const std::vector<T> empty_vec{};

  const bool te_kind = kind == TripletWeightKind::TeNnsReconstruction ||
                       kind == TripletWeightKind::TeReconstruction ||
                       kind == TripletWeightKind::TeCombinedResidual;
  if ((n_particles != 2 && n_particles != 3) || (te_kind && n_particles != 2)) {
    return empty_vec;
  }

  using CacheKey = std::pair<std::size_t, TripletWeightKind>;
  // the rows are cached behind a pointer since the cache holds several keys
  // and flat_map insertions would invalidate references into it
  using CachedRow = std::unique_ptr<const std::vector<T>>;

  static std::mutex cache_mutex;
  static std::condition_variable cache_cv;
  static container::map<CacheKey, std::optional<CachedRow>> cache;

  CacheKey key{n_particles, kind};

  return *sequant::detail::memoize(
      cache, cache_mutex, cache_cv, key, [&]() -> CachedRow {
        auto coeffs = compute_triplet_weights(n_particles, kind);
        std::vector<T> weights;
        weights.reserve(coeffs.size());
        for (const auto& c : coeffs) {
          weights.push_back(static_cast<T>(c));
        }
        return std::make_unique<const std::vector<T>>(std::move(weights));
      });
}

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

/// \brief Weighted sum over the S_n x S_n external-slot permutations of a
/// rank-2n TA::DistArray:
///   out(orig) = sum_p weights[p] * arr(annot_p),
/// with the annotations generated by slot_perm_annots and the weight row
/// selected by the caller (see TripletWeightKind). Zero-weight
/// permutations are skipped.
template <typename... Args>
auto triplet_perm_combine_ta(
    TA::DistArray<Args...> const& arr, std::string const& orig_layout,
    std::size_t n_particles,
    std::vector<typename TA::DistArray<Args...>::numeric_type> const& weights) {
  using numeric_type = typename TA::DistArray<Args...>::numeric_type;
  const auto annots = slot_perm_annots(orig_layout, n_particles);
  SEQUANT_ASSERT(annots.size() == weights.size());

  TA::DistArray<Args...> result;
  bool result_initialized = false;
  for (std::size_t p = 0; p != weights.size(); ++p) {
    if (weights[p] == numeric_type(0)) continue;
    if (result_initialized) {
      result(orig_layout) += weights[p] * arr(annots[p]);
    } else {
      result(orig_layout) = weights[p] * arr(annots[p]);
      result_initialized = true;
    }
  }
  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());
  result.truncate();
  return result;
}

/// \brief Applies the \p kind weight row over the triplet slot permutations of
/// \p arr; no-op unless the array rank is 4 or 6 (n_particles = rank/2) and
/// the row is available for that rank
template <typename... Args>
auto triplet_perm_project_ta(TA::DistArray<Args...> const& arr,
                             std::string const& orig_layout,
                             TripletWeightKind kind) {
  using numeric_type = typename TA::DistArray<Args...>::numeric_type;
  const std::size_t rank = arr.trange().rank();
  if (rank != 4 && rank != 6) return arr;
  const auto& weights = triplet_weights<numeric_type>(rank / 2, kind);
  if (weights.empty()) return arr;
  return triplet_perm_combine_ta<Args...>(arr, orig_layout, rank / 2, weights);
}

}  // namespace detail

/// \brief Idempotent null-space projector for the closed-shell triplet R:
/// removes the metric-null component of the array. Apply to the Davidson
/// trial vector each iteration; no-op unless the array rank is 4 or 6.
template <typename... Args>
auto triplet_nullspace_project_ta(TA::DistArray<Args...> const& arr,
                                  std::string const& orig_layout) {
  return detail::triplet_perm_project_ta(arr, orig_layout,
                                         TripletWeightKind::NullspaceProjector);
}

template <typename... Args>
auto triplet_nullspace_project(TA::DistArray<Args...> const& arr,
                               std::string const& orig_layout) {
  return triplet_nullspace_project_ta(arr, orig_layout);
}

/// \brief Metric NNS reconstruction for compact closed-shell triplet
/// residuals: rebuilds the full residual from the representatives kept by
/// triplet_maxcoeff_compact (numerical analogue of
/// triplet_symbolic_reconstruct). Apply to the H*R residual when the compact
/// equations were evaluated; no-op unless the array rank is 4 or 6.
template <typename... Args>
auto triplet_nns_project_ta(TA::DistArray<Args...> const& arr,
                            std::string const& orig_layout) {
  return detail::triplet_perm_project_ta(arr, orig_layout,
                                         TripletWeightKind::NnsReconstruction);
}

template <typename... Args>
auto triplet_nns_project(TA::DistArray<Args...> const& arr,
                         std::string const& orig_layout) {
  return triplet_nns_project_ta(arr, orig_layout);
}

/// \brief Bare-TE undo-compact for compact triplet R2 residuals
/// (TeNnsReconstruction row). Apply to the H*R residual when the compact
/// te_only equations were evaluated, before triplet_te_reconstruct; no-op
/// unless the array rank is 4.
template <typename... Args>
auto triplet_te_nns_project_ta(TA::DistArray<Args...> const& arr,
                               std::string const& orig_layout) {
  return detail::triplet_perm_project_ta(
      arr, orig_layout, TripletWeightKind::TeNnsReconstruction);
}

template <typename... Args>
auto triplet_te_nns_project(TA::DistArray<Args...> const& arr,
                            std::string const& orig_layout) {
  return triplet_te_nns_project_ta(arr, orig_layout);
}

/// \brief TE-only -> full-metric reconstruction for triplet R2 (EFV
/// experiment), via the exact identity
/// Omega = te_a + (1/4)(bra_swap(te_a) + ket_swap(te_a)).
template <typename... Args>
auto triplet_te_reconstruct_ta(TA::DistArray<Args...> const& arr,
                               std::string const& orig_layout) {
  return detail::triplet_perm_project_ta(arr, orig_layout,
                                         TripletWeightKind::TeReconstruction);
}

template <typename... Args>
auto triplet_te_reconstruct(TA::DistArray<Args...> const& arr,
                            std::string const& orig_layout) {
  return triplet_te_reconstruct_ta(arr, orig_layout);
}

#endif  // defined(SEQUANT_HAS_TILEDARRAY)

#if defined(SEQUANT_HAS_BTAS)

namespace detail {

/// \brief BTAS analogue of triplet_perm_combine_ta
template <typename... Args>
auto triplet_perm_combine_btas(
    btas::Tensor<Args...> const& arr, std::size_t n_particles,
    std::vector<typename btas::Tensor<Args...>::numeric_type> const& weights) {
  using ranges::views::iota;
  using numeric_type = typename btas::Tensor<Args...>::numeric_type;
  const std::size_t rank = 2 * n_particles;

  sequant::detail::perm_t perm =
      iota(size_t{0}, rank) | ranges::to<sequant::detail::perm_t>;

  btas::Tensor<Args...> result;
  bool result_initialized = false;
  for (std::size_t p = 0; p != weights.size(); ++p) {
    if (weights[p] == numeric_type(0)) continue;
    const sequant::detail::perm_t annot =
        compute_bra_ket_permuted_indices(p, rank / 2);

    btas::Tensor<Args...> temp;
    btas::permute(arr, annot, temp, perm);
    btas::scal(weights[p], temp);

    if (result_initialized) {
      result += temp;
    } else {
      result = temp;
      result_initialized = true;
    }
  }
  return result;
}

/// \brief BTAS analogue of triplet_perm_project_ta; \p orig_layout is
/// accepted for signature parity with the TA overloads but not consulted
/// (slots are positional)
template <typename... Args>
auto triplet_perm_project_btas(btas::Tensor<Args...> const& arr,
                               [[maybe_unused]] std::string const& orig_layout,
                               TripletWeightKind kind) {
  using numeric_type = typename btas::Tensor<Args...>::numeric_type;
  const std::size_t rank = arr.rank();
  if (rank != 4 && rank != 6) return arr;
  const auto& weights = triplet_weights<numeric_type>(rank / 2, kind);
  if (weights.empty()) return arr;
  return triplet_perm_combine_btas<Args...>(arr, rank / 2, weights);
}

}  // namespace detail

/// @brief BTAS analogue of the TA triplet_nns_project
template <typename... Args>
auto triplet_nns_project_btas(btas::Tensor<Args...> const& arr,
                              std::string const& orig_layout) {
  return detail::triplet_perm_project_btas(
      arr, orig_layout, TripletWeightKind::NnsReconstruction);
}

template <typename... Args>
auto triplet_nns_project(btas::Tensor<Args...> const& arr,
                         std::string const& orig_layout) {
  return triplet_nns_project_btas(arr, orig_layout);
}

/// @brief BTAS analogue of the TA triplet_nullspace_project
template <typename... Args>
auto triplet_nullspace_project_btas(btas::Tensor<Args...> const& arr,
                                    std::string const& orig_layout) {
  return detail::triplet_perm_project_btas(
      arr, orig_layout, TripletWeightKind::NullspaceProjector);
}

template <typename... Args>
auto triplet_nullspace_project(btas::Tensor<Args...> const& arr,
                               std::string const& orig_layout) {
  return triplet_nullspace_project_btas(arr, orig_layout);
}

#endif  // defined(SEQUANT_HAS_BTAS)

}  // namespace sequant::mbpt

#endif
