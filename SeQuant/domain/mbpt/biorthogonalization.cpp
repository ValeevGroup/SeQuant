#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/detail/concepts.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/permutation.hpp>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <libperm/Permutation.hpp>
#include <libperm/Rank.hpp>
#include <libperm/Utils.hpp>

#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

namespace sequant::mbpt {

template <typename T>
struct compare_first_less {
  bool operator()(const T& lhs, const T& rhs) const {
    return lhs.first < rhs.first;
  }
};

using IndexPair = std::pair<Index, Index>;
using ParticlePairings = container::svector<IndexPair>;

// clang-format off
/// \brief Provides the first row of the biorthogonal coefficients matrix,
/// hardcoded from Mathematica to avoid numerical precision loss.
///
/// The Myrvold-Ruskey unrank1 algorithm (doi.org/10.1016/S0020-0190(01)00141-7)
/// is used to order permutations, then the permutational overlap matrix M is
/// constructed with elements (-2)^{c} × (-1)^{n_particles}, where c is the
/// number of cycles in the relative permutation.
///
/// The biorthogonal coefficients are obtained from the unnormalized
/// pseudoinverse of M.
/// biorthogonal coefficients = M_pinv · e_1, where e_1 is the first unit vector.
/// See [DOI 10.48550/ARXIV.1805.00565](https://doi.org/10.48550/ARXIV.1805.00565)
/// for more details (Note that these biorthogonal coefficients are unnormalized,
/// unlike those in the reference).
///
/// \param n_particles The rank of external index pairs
///
/// \return Vector of rational coefficients representing the first row
///
/// \throw Exception if n_particles is not in the range [1,5]
// clang-format on
std::vector<sequant::rational> hardcoded_biorthogonalizer_row(
    std::size_t n_particles) {
  switch (n_particles) {
    case 1:
      return std::vector<sequant::rational>{ratio(1, 2)};

    case 2:
      return std::vector<sequant::rational>{ratio(1, 3), ratio(1, 6)};

    case 3:
      return std::vector<sequant::rational>{ratio(17, 144), ratio(-7, 144),
                                            ratio(-1, 144), ratio(-1, 144),
                                            ratio(-1, 144), ratio(-7, 144)};

    case 4:
      return std::vector<sequant::rational>{
          ratio(43, 1440), ratio(-19, 2880), ratio(-19, 2880),
          ratio(-1, 180),  ratio(-19, 2880), ratio(-19, 2880),
          ratio(13, 1440), ratio(7, 1440),   ratio(-1, 180),
          ratio(7, 1440),  ratio(-1, 180),   ratio(-19, 2880),
          ratio(-1, 180),  ratio(7, 1440),   ratio(7, 1440),
          ratio(13, 1440), ratio(-1, 180),   ratio(-1, 180),
          ratio(7, 1440),  ratio(-19, 2880), ratio(-19, 2880),
          ratio(13, 1440), ratio(-19, 2880), ratio(7, 1440)};

    case 5:
      return std::vector<sequant::rational>{
          ratio(59, 10800),  ratio(-1, 1728),   ratio(-1, 1728),
          ratio(-1, 1728),   ratio(-31, 21600), ratio(-1, 1728),
          ratio(-1, 1728),   ratio(-23, 86400), ratio(19, 21600),
          ratio(37, 43200),  ratio(-1, 1728),   ratio(-23, 86400),
          ratio(-1, 1728),   ratio(19, 21600),  ratio(37, 43200),
          ratio(-31, 21600), ratio(37, 43200),  ratio(37, 43200),
          ratio(-31, 21600), ratio(-1, 1728),   ratio(-1, 1728),
          ratio(-23, 86400), ratio(-23, 86400), ratio(-23, 86400),
          ratio(-13, 21600), ratio(-1, 1728),   ratio(-1, 1728),
          ratio(19, 21600),  ratio(-23, 86400), ratio(37, 43200),
          ratio(19, 21600),  ratio(-23, 86400), ratio(19, 21600),
          ratio(-23, 86400), ratio(-13, 21600), ratio(37, 43200),
          ratio(-13, 21600), ratio(-13, 21600), ratio(37, 43200),
          ratio(-23, 86400), ratio(-31, 21600), ratio(-13, 21600),
          ratio(37, 43200),  ratio(37, 43200),  ratio(19, 21600),
          ratio(37, 43200),  ratio(37, 43200),  ratio(-13, 21600),
          ratio(-13, 21600), ratio(-23, 86400), ratio(-31, 21600),
          ratio(37, 43200),  ratio(-31, 21600), ratio(37, 43200),
          ratio(-1, 1728),   ratio(-1, 1728),   ratio(-23, 86400),
          ratio(19, 21600),  ratio(-1, 1728),   ratio(37, 43200),
          ratio(-31, 21600), ratio(37, 43200),  ratio(37, 43200),
          ratio(-13, 21600), ratio(19, 21600),  ratio(37, 43200),
          ratio(37, 43200),  ratio(-13, 21600), ratio(-13, 21600),
          ratio(-23, 86400), ratio(37, 43200),  ratio(-13, 21600),
          ratio(37, 43200),  ratio(-13, 21600), ratio(-23, 86400),
          ratio(19, 21600),  ratio(-23, 86400), ratio(-23, 86400),
          ratio(19, 21600),  ratio(-13, 21600), ratio(-31, 21600),
          ratio(37, 43200),  ratio(-13, 21600), ratio(37, 43200),
          ratio(19, 21600),  ratio(-31, 21600), ratio(-31, 21600),
          ratio(37, 43200),  ratio(37, 43200),  ratio(-1, 1728),
          ratio(37, 43200),  ratio(-13, 21600), ratio(37, 43200),
          ratio(-13, 21600), ratio(-23, 86400), ratio(-1, 1728),
          ratio(19, 21600),  ratio(-23, 86400), ratio(-1, 1728),
          ratio(37, 43200),  ratio(-1, 1728),   ratio(-23, 86400),
          ratio(-23, 86400), ratio(-23, 86400), ratio(-13, 21600),
          ratio(19, 21600),  ratio(19, 21600),  ratio(-23, 86400),
          ratio(-23, 86400), ratio(-13, 21600), ratio(-1, 1728),
          ratio(19, 21600),  ratio(-1, 1728),   ratio(-23, 86400),
          ratio(37, 43200),  ratio(37, 43200),  ratio(-13, 21600),
          ratio(-13, 21600), ratio(37, 43200),  ratio(-23, 86400)};

    default:
      throw Exception(
          "hardcoded biorthogonal coefficients only available for ranks 1-5, "
          "requested rank is : " +
          std::to_string(n_particles));
  }
}

Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
make_hardcoded_biorthogonalizer_matrix(
    const std::vector<sequant::rational>& first_row, std::size_t n_particles) {
  const auto n = first_row.size();
  Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic> M(n, n);

  // entry (row, col) depends only on the permutation relating the permutations
  // of rank row and col (formed as in permutational_overlap_matrix), and
  // first_row tabulates it for row 0, so look up the column whose relation to
  // row 0 is the same
  const perm::Permutation first = perm::unrank(0, n_particles);
  for (std::size_t row = 0; row < n; ++row) {
    perm::Permutation ref = perm::unrank(row, n_particles);
    ref->invert();

    for (std::size_t col = 0; col < n; ++col) {
      perm::Permutation current = perm::unrank(col, n_particles);
      current->postMultiply(ref);
      current->postMultiply(first);

      std::size_t source_idx = perm::rank(current, n_particles);
      M(row, col) = first_row[source_idx];
    }
  }
  return M;
}

Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
hardcoded_biorthogonalizer_matrix(std::size_t n_particles) {
  auto first_row = hardcoded_biorthogonalizer_row(n_particles);
  return make_hardcoded_biorthogonalizer_matrix(first_row, n_particles);
}

ResultExpr biorthogonal_transform_copy(
    const ResultExpr& expr,
    double threshold = default_biorthogonalizer_pseudoinverse_threshold) {
  container::svector<ResultExpr> wrapper = {expr.clone()};

  biorthogonal_transform(wrapper, threshold);

  return wrapper.front();
}

container::svector<ResultExpr> biorthogonal_transform_copy(
    const container::svector<ResultExpr>& exprs,
    double threshold = default_biorthogonalizer_pseudoinverse_threshold) {
  container::svector<ResultExpr> copy;
  copy.reserve(exprs.size());

  std::transform(exprs.begin(), exprs.end(), std::back_inserter(copy),
                 [](const ResultExpr& expr) { return expr.clone(); });

  biorthogonal_transform(copy, threshold);

  return copy;
}

void biorthogonal_transform(ResultExpr& expr, double threshold) {
  // TODO: avoid copy
  expr = biorthogonal_transform_copy(expr, threshold);
}

Eigen::MatrixXd permutational_overlap_matrix(std::size_t n_particles) {
  const auto n = boost::numeric_cast<Eigen::Index>(factorial(n_particles));

  // The matrix only contains integer entries but all operations we want to do
  // with the matrix will (in general) require non-integer scalars which in
  // Eigen only works if you start from a non-integer matrix.
  Eigen::MatrixXd M(n, n);
  M.setZero();

  // TODO: Can we fill the entire matrix only by knowing the entries of one
  // row/column? For n_particles < 4, every consecutive col/row is only rotated
  // by one compared to the one before
  for (std::size_t row = 0; row < n; ++row) {
    perm::Permutation ref = perm::unrank(row, n_particles);
    ref->invert();

    // The identity permutation always has as many disjoint cycles as the number
    // of elements it acts on
    M(row, row) = std::pow(-2, n_particles);

    for (std::size_t col = row + 1; col < n; ++col) {
      // Get permutation that transforms the permutation of rank1 into the one
      // of current rank i
      perm::Permutation current = perm::unrank(col, n_particles);
      current->postMultiply(ref);

      auto cycles = current->toDisjointCycles(n_particles);
      std::size_t n_cycles = std::distance(cycles.begin(), cycles.end());

      auto entry = std::pow(-2, n_cycles);

      M(row, col) = entry;
      M(col, row) = entry;
    }
  }

  if (n_particles % 2 != 0) {
    M *= -1;
  }

  SEQUANT_ASSERT(M.isApprox(M.transpose()));

  return M;
}

/// rank of the permutational overlap matrix of @p n_particles particles = the
/// number of linearly independent singlet couplings = Catalan number C_n
sequant::intmax_t catalan(std::size_t n) {
  return factorial(2 * n) / (factorial(n + 1) * factorial(n));
}

Eigen::MatrixXd compute_biorthogonalizer_matrix(std::size_t n_particles,
                                                double threshold) {
  auto perm_ovlp_mat = permutational_overlap_matrix(n_particles);
  SEQUANT_ASSERT(perm_ovlp_mat.rows() == perm_ovlp_mat.cols());
  SEQUANT_ASSERT(perm_ovlp_mat.isApprox(perm_ovlp_mat.transpose()));

  // Find Pseudo Inverse
  auto decomp =
      Eigen::CompleteOrthogonalDecomposition<decltype(perm_ovlp_mat)>();
  decomp.setThreshold(threshold);
  decomp.compute(perm_ovlp_mat);
  // WK_biorthogonalization_filter_impl assumes this rank
  SEQUANT_ASSERT(decomp.rank() == catalan(n_particles));

  Eigen::MatrixXd pinv = decomp.pseudoInverse();
  // The pseudo inverse of a symmetric matrix should also be symmetric
  SEQUANT_ASSERT(pinv.isApprox(pinv.transpose()));

  return pinv;
}

void sort_pairings(ParticlePairings& pairing) {
  std::stable_sort(pairing.begin(), pairing.end(),
                   compare_first_less<IndexPair>{});
}

std::size_t rank_transformation_perms(const ParticlePairings& reference,
                                      const ParticlePairings& current) {
  SEQUANT_ASSERT(reference.size() == current.size());
  SEQUANT_ASSERT(std::is_sorted(reference.begin(), reference.end(),
                                compare_first_less<IndexPair>{}));
  SEQUANT_ASSERT(std::is_sorted(current.begin(), current.end(),
                                compare_first_less<IndexPair>{}));

  perm::Permutation perm = perm::computeTransformationPermutation(
      reference, current, [](const IndexPair& lhs, const IndexPair& rhs) {
        return lhs.second < rhs.second;
      });

  return perm::rank(perm, reference.size());
}

ExprPtr create_expr_for(const ParticlePairings& ref_pairing,
                        const perm::Permutation& perm,
                        const container::svector<ParticlePairings>& pairings,
                        const container::svector<ExprPtr>& base_exprs) {
  // Note: perm only applies to the p->second for every pair p in ref_pairing

  // assert that all pairings are sorted w.r.t. first
  SEQUANT_ASSERT(std::all_of(
      pairings.begin(), pairings.end(), [](const ParticlePairings& pairing) {
        return std::is_sorted(pairing.begin(), pairing.end(),
                              compare_first_less<IndexPair>{});
      }));
  SEQUANT_ASSERT(std::is_sorted(ref_pairing.begin(), ref_pairing.end(),
                                compare_first_less<IndexPair>{}));

  container::set<std::pair<IndexSpace, IndexSpace>> ref_space_pairing;
  ref_space_pairing.reserve(ref_pairing.size());
  for (std::size_t i = 0; i < ref_pairing.size(); ++i) {
    ref_space_pairing.emplace(ref_pairing[i].first.space(),
                              ref_pairing[perm->image(i)].second.space());
  }

  // Look for a ParticlePairings object that pairs indices belonging to index
  // spaces compatible with ref_space_pairing
  auto it = std::find_if(
      pairings.begin(), pairings.end(), [&](const ParticlePairings& p) {
        SEQUANT_ASSERT(p.size() == ref_pairing.size());

        for (const IndexPair& pair : p) {
          if (ref_space_pairing.find(
                  std::make_pair(pair.first.space(), pair.second.space())) ==
              ref_space_pairing.end()) {
            return false;
          }
        }

        return true;
      });

  if (it == pairings.end()) {
    throw Exception(
        "Missing explicit expression for a required index pairing in "
        "biorthogonalization");
  }

  auto idx = std::distance(pairings.begin(), it);
  const ParticlePairings& base = *it;

  SEQUANT_ASSERT(base.size() == ref_pairing.size());

  container::map<Index, Index> replacements;
  for (std::size_t i = 0; i < base.size(); ++i) {
    std::size_t ref_idx = perm->image(i);

    // Remember that all index pairings are sorted w.r.t. first and hence we are
    // only looking for permutations in second
    SEQUANT_ASSERT(base[i].first == ref_pairing[i].first);
    const bool differs_in_second =
        base[i].second != ref_pairing[ref_idx].second;

    if (!differs_in_second) {
      // This particle pairing is identical
      continue;
    }

    SEQUANT_ASSERT(differs_in_second);

    // Note: we may only permute indices belonging to the same space
    // (otherwise, we would produce non-sensical expressions)
    if (base[i].second.space() == ref_pairing[ref_idx].second.space()) {
      // base and ref_pairing differ in the second index of the current
      // pairing and their index space matches -> can just permute them
      replacements.emplace(base[i].second, ref_pairing[ref_idx].second);
    } else {
      // Index spaces of the differing index (second) in the pairings are
      // different as well. Since the tensors are assumed to be
      // particle-symmetric, we can instead permute the first indices in the
      // pairings, which are of the same space (that's guaranteed by the way we
      // chose base).
      SEQUANT_ASSERT(base[i].first.space() ==
                     ref_pairing[ref_idx].first.space());
      replacements.emplace(base[i].first, ref_pairing[ref_idx].first);
    }
  }

  ExprPtr expr = base_exprs.at(idx)->clone();

  if (!replacements.empty()) {
    if constexpr (assert_enabled()) {
      for ([[maybe_unused]] const auto& [first, second] : replacements) {
        SEQUANT_ASSERT(first.space() == second.space());
      }
    }
    expr = transform_expr(expr, replacements);
  }

  return expr;
}

void biorthogonal_transform(container::svector<ResultExpr>& result_exprs,
                            double threshold) {
  if (result_exprs.empty()) {
    return;
  }

  // We expect all ResultExpr objects to be equal except for the permutation of
  // indices
  // Also, we are assuming that all given ResultExpr objects are
  // particle-symmetric
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.has_label() == result_exprs.front().has_label() &&
               (!expr.has_label() ||
                expr.label() == result_exprs.front().label());
      }));
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.symmetry() == result_exprs.front().symmetry();
      }));
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.braket_symmetry() == result_exprs.front().braket_symmetry();
      }));
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.column_symmetry() == result_exprs.front().column_symmetry();
      }));
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.bra().size() == result_exprs.front().bra().size() &&
               std::is_permutation(expr.bra().begin(), expr.bra().end(),
                                   result_exprs.front().bra().begin());
      }));
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.ket().size() == result_exprs.front().ket().size() &&
               std::is_permutation(expr.ket().begin(), expr.ket().end(),
                                   result_exprs.front().ket().begin());
      }));
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.aux().size() == result_exprs.front().aux().size() &&
               std::is_permutation(expr.aux().begin(), expr.aux().end(),
                                   result_exprs.front().aux().begin());
      }));
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [](const ResultExpr& res) {
        return res.column_symmetry() == ColumnSymmetry::Symm;
      }));

  // Furthermore, we expect that there is no symmetrization operator present in
  // the expressions as that would imply transforming also the symmetrization
  // operator, which is incorrect. This is because the idea during
  // biorthogonalization is that we project onto e.g.
  // \tilde{E}^{IJ}_{AB} = c_1 E^{IJ}_{AB} + c_2 E^{JI}_{AB}
  // instead of E^{IJ}_{AB} directly. In either case though, the result looks
  // like R^{IJ}_{AB} and the index pairing of the result is what determines
  // the required symmetrization. Hence, the symmetrization operator must not
  // be changed when transforming from one representation into the other.
  SEQUANT_ASSERT(std::all_of(
      result_exprs.begin(), result_exprs.end(), [](const ResultExpr& res) {
        bool found = false;
        res.expression()->visit(
            [&](const ExprPtr& expr) {
              if (expr->is<Tensor>() &&
                  (expr->as<Tensor>().label() == reserved::symm_label() ||
                   expr->as<Tensor>().label() == reserved::antisymm_label())) {
                found = true;
              };
            },
            true);
        return !found;
      }));

  auto externals = result_exprs |
                   ranges::views::transform([](const ResultExpr& expr) {
                     return expr.index_particle_grouping<IndexPair>();
                   }) |
                   ranges::to<container::svector<ParticlePairings>>();
  ranges::for_each(externals, sort_pairings);

  auto ranks = externals | ranges::views::transform([&](const auto& p) {
                 return rank_transformation_perms(externals.front(), p);
               }) |
               ranges::to<container::svector<std::size_t>>();

  const std::size_t n_particles = externals.front().size();
  auto num_perms = factorial(n_particles);

  auto original_exprs = result_exprs |
                        ranges::views::transform([](const ResultExpr& res) {
                          return res.expression();
                        }) |
                        ranges::to<container::svector<ExprPtr>>();

  using HardcodedMatrix =
      Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>;
  using ComputedMatrix = Eigen::MatrixXd;
  using CacheKey = std::pair<std::size_t, double>;

  static std::mutex cache_mutex;
  static std::condition_variable cache_cv;
  static container::map<CacheKey, std::optional<HardcodedMatrix>>
      hardcoded_cache;
  static container::map<CacheKey, std::optional<ComputedMatrix>> computed_cache;

  constexpr std::size_t max_rank_hardcoded_biorthogonalizer_matrix = 5;
  CacheKey key{n_particles, threshold};

  const HardcodedMatrix* hardcoded_coefficients = nullptr;
  const ComputedMatrix* computed_coefficients = nullptr;

  if (n_particles <= max_rank_hardcoded_biorthogonalizer_matrix) {
    hardcoded_coefficients = &sequant::detail::memoize(
        hardcoded_cache, cache_mutex, cache_cv, key,
        [&] { return hardcoded_biorthogonalizer_matrix(n_particles); });
  } else {
    computed_coefficients = &sequant::detail::memoize(
        computed_cache, cache_mutex, cache_cv, key, [&] {
          return compute_biorthogonalizer_matrix(n_particles, threshold);
        });
    SEQUANT_ASSERT(num_perms == computed_coefficients->rows());
    SEQUANT_ASSERT(num_perms == computed_coefficients->cols());
  }

  const auto perm_ranks =
      ranges::views::iota(std::size_t{0}, static_cast<std::size_t>(num_perms));
  for (std::size_t i = 0; i < result_exprs.size(); ++i) {
    perm::Permutation reference = perm::unrank(ranks.at(i), n_particles);
    reference->invert();

    // Merge canonical terms before retaining the next permuted expression.
    result_exprs.at(i).expression() = transform_sum_expr(
        perm_ranks,
        [&](std::size_t rank) {
          perm::Permutation perm = perm::unrank(rank, n_particles);
          perm->postMultiply(reference);

          sequant::rational coeff =
              (n_particles <= max_rank_hardcoded_biorthogonalizer_matrix)
                  ? (*hardcoded_coefficients)(ranks.at(i), rank)
                  : to_rational((*computed_coefficients)(ranks.at(i), rank),
                                threshold);

          auto contribution =
              ex<Constant>(coeff) *
              create_expr_for(externals.at(i), perm, externals, original_exprs);
          simplify(contribution);
          return contribution;
        },
        {.canonicalize = false});
    simplify(result_exprs.at(i).expression());
  }
}

template <detail::index_group_range IdxGroups>
ExprPtr biorthogonal_transform_impl(const sequant::ExprPtr& expr,
                                    IdxGroups&& ext_index_groups,
                                    const double threshold) {
  ResultExpr res(
      bra(ext_index_groups | ranges::views::transform([](const auto& pair) {
            return get_ket_idx(pair);
          }) |
          ranges::to<container::svector<Index>>()),
      ket(ext_index_groups | ranges::views::transform([](const auto& pair) {
            return get_bra_idx(pair);
          }) |
          ranges::to<container::svector<Index>>()),
      aux(IndexList{}), Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
      ColumnSymmetry::Symm, {}, expr);

  biorthogonal_transform(res, threshold);

  return res.expression();
}

ExprPtr biorthogonal_transform(
    const sequant::ExprPtr& expr,
    const container::svector<container::svector<sequant::SlottedIndex>>&
        ext_index_groups,
    const double threshold) {
  return biorthogonal_transform_impl(
      expr, as_view_of_index_groups(ext_index_groups), threshold);
}

ExprPtr biorthogonal_transform(
    const sequant::ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups,
    const double threshold) {
  return biorthogonal_transform_impl(expr, ext_index_groups, threshold);
}

namespace {

// groups the terms of \p expr by tensor-network hash, keeping in each group
// the terms with the largest |coefficient|
container::map<std::size_t, container::vector<ExprPtr>>
largest_coeff_terms_by_hash(const ExprPtr& expr) {
  container::map<std::size_t, container::vector<ExprPtr>> largest_coeff_terms;

  const auto ctx = get_default_context_snapshot();
  const auto& cardinal_tensor_labels = ctx.cardinal_tensor_labels();
  for (const auto& term : *expr) {
    if (!term->is<Product>()) continue;

    auto product = term.as_shared_ptr<Product>();
    auto scalar = product->scalar();

    sequant::TensorNetwork tn(*product);
    auto hash = tn.canonicalize_slots(cardinal_tensor_labels).hash_value();

    auto it = largest_coeff_terms.find(hash);
    if (it == largest_coeff_terms.end()) {
      largest_coeff_terms[hash] = {term};
    } else {
      if (!it->second.empty()) {
        auto existing_scalar = it->second[0]->as<Product>().scalar();
        auto existing_abs = abs(existing_scalar);
        auto current_abs = abs(scalar);

        if (current_abs > existing_abs) {
          it->second.clear();
          it->second.push_back(term);
        } else if (current_abs == existing_abs) {
          it->second.push_back(term);
        }
      }
    }
  }

  return largest_coeff_terms;
}

}  // namespace

template <detail::index_group_range IdxGroups>
ExprPtr WK_biorthogonalization_filter_impl(ExprPtr expr, IdxGroups&& ext_idxs) {
  if (!expr->is<Sum>()) return expr;
  if (ext_idxs.size() <= 2) return expr;  // always skip R1 and R2

  // the diagonal of the idempotent NNS projector is rank(M)/n!, so each kept
  // term carries the reciprocal for the projector to rebuild its group at full
  // weight; all terms tied for the largest coefficient are kept and scaled
  // alike
  const auto n = ext_idxs.size();
  const auto rank = catalan(n);
  const rational normalization(factorial(n), rank);
  Sum filtered;
  for (const auto& [_, terms] : largest_coeff_terms_by_hash(expr)) {
    for (const auto& t : terms) {
      auto kept = ex<Product>(t->as<Product>());
      kept->as<Product>().scale(normalization);
      filtered.append(std::move(kept));
    }
  }
  auto result = ex<Sum>(filtered);

  return result;
}

ExprPtr WK_biorthogonalization_filter(
    ExprPtr expr,
    const container::svector<container::svector<SlottedIndex>>& ext_idxs) {
  return WK_biorthogonalization_filter_impl(expr,
                                            as_view_of_index_groups(ext_idxs));
}
ExprPtr WK_biorthogonalization_filter(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs) {
  return WK_biorthogonalization_filter_impl(expr, ext_idxs);
}

template <detail::index_group_range IdxGroups>
ExprPtr biorthogonal_transform_pre_nnsproject_impl(
    ExprPtr& expr, IdxGroups&& ext_idxs, bool factor_out_nns_projector) {
  using ranges::views::transform;

  // Remove leading S operator if present
  for (auto& term : *expr) {
    if (term->is<Product>())
      term =
          remove_tensor(term.as_shared_ptr<Product>(), reserved::symm_label());
  }

  auto bt = biorthogonal_transform_impl(
      expr, ext_idxs, default_biorthogonalizer_pseudoinverse_threshold);

  auto bixs = ext_idxs | transform([](auto&& vec) { return get_bra_idx(vec); });
  auto kixs = ext_idxs | transform([](auto&& vec) { return get_ket_idx(vec); });
  ExprPtr S_tensor = make_symmetrizer(bra(kixs), ket(bixs));

  if (factor_out_nns_projector) {
    if (ext_idxs.size() > 1) {
      bt = S_tensor * bt;
    }
    simplify(bt);

    bt = S_maps(bt);
    canonicalize(bt);
    bt = WK_biorthogonalization_filter_impl(bt, ext_idxs);
  }

  bt = S_tensor * bt;
  simplify(bt);

  return bt;
}

ExprPtr biorthogonal_transform_pre_nnsproject(
    ExprPtr& expr,
    const container::svector<container::svector<SlottedIndex>>& ext_idxs,
    bool factor_out_nns_projector) {
  return biorthogonal_transform_pre_nnsproject_impl(
      expr, as_view_of_index_groups(ext_idxs), factor_out_nns_projector);
}

ExprPtr biorthogonal_transform_pre_nnsproject(
    ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_idxs,
    bool factor_out_nns_projector) {
  return biorthogonal_transform_pre_nnsproject_impl(expr, ext_idxs,
                                                    factor_out_nns_projector);
}

namespace detail {

std::vector<double> compute_nns_p_coeffs(std::size_t n_particles,
                                         double threshold) {
  auto perm_ovlp_mat = permutational_overlap_matrix(n_particles);
  auto pinv = compute_biorthogonalizer_matrix(n_particles, threshold);
  Eigen::MatrixXd nns_matrix = perm_ovlp_mat * pinv;
  SEQUANT_ASSERT((nns_matrix * nns_matrix).isApprox(nns_matrix));

  auto num_perms = nns_matrix.rows();
  std::vector<double> coeffs;
  coeffs.reserve(num_perms);
  for (std::size_t i = 0; i < num_perms; ++i) {
    coeffs.push_back(nns_matrix(num_perms - 1, i));
  }
  return coeffs;
}

template <typename Matrix>
auto flatten_row_major(const Matrix& matrix) {
  std::vector<typename Matrix::Scalar> coeffs;
  coeffs.reserve(matrix.size());
  for (Eigen::Index row = 0; row < matrix.rows(); ++row)
    for (Eigen::Index col = 0; col < matrix.cols(); ++col)
      coeffs.push_back(matrix(row, col));
  return coeffs;
}

std::vector<double> compute_biorthogonalizer_coeffs(std::size_t n_particles,
                                                    double threshold) {
  return flatten_row_major(
      compute_biorthogonalizer_matrix(n_particles, threshold));
}

std::vector<rational> hardcoded_biorthogonalizer_coeffs(
    std::size_t n_particles) {
  return flatten_row_major(hardcoded_biorthogonalizer_matrix(n_particles));
}

container::svector<size_t> compute_permuted_indices(
    const container::svector<size_t>& indices, size_t perm_rank,
    size_t n_particles) {
  perm::Permutation perm_obj = perm::unrank(perm_rank, n_particles);

  container::svector<size_t> permuted_indices(n_particles);
  for (size_t i = 0; i < n_particles; ++i) {
    permuted_indices[i] = indices[perm_obj[i]];
  }
  return permuted_indices;
}

}  // namespace detail

// clang-format off
/// \brief Provides the closed-shell triplet residual assembly row over the
/// (n!)^2 external slot permutations
///
/// A combination of T and E operators with their permutations
/// [Kohn's triplet paper](http://dx.doi.org/10.1063/1.1457434).
///
/// \param n_particles The rank of external index pairs
///
/// \return Vector of rational weights in flat perm order
///
/// \throw Exception if the row is not available for \p n_particles
// clang-format on
std::vector<rational> hardcoded_triplet_residual_row(std::size_t n_particles) {
  switch (n_particles) {
    case 1:
      return {ratio(1, 2)};

    case 2:
      return {ratio(-1, 16), ratio(0, 1), ratio(0, 1), ratio(3, 16)};

    case 3:
      return {ratio(0, 1), ratio(0, 1),    ratio(0, 1),    ratio(0, 1),
              ratio(0, 1), ratio(0, 1),    ratio(0, 1),    ratio(0, 1),
              ratio(0, 1), ratio(0, 1),    ratio(0, 1),    ratio(0, 1),
              ratio(0, 1), ratio(0, 1),    ratio(-1, 160), ratio(0, 1),
              ratio(0, 1), ratio(0, 1),    ratio(0, 1),    ratio(0, 1),
              ratio(0, 1), ratio(-1, 160), ratio(0, 1),    ratio(0, 1),
              ratio(0, 1), ratio(0, 1),    ratio(0, 1),    ratio(0, 1),
              ratio(0, 1), ratio(0, 1),    ratio(0, 1),    ratio(0, 1),
              ratio(0, 1), ratio(0, 1),    ratio(1, 80),   ratio(3, 80)};
  }
  throw Exception(
      "hardcoded triplet residual weights are not available for "
      "n_particles = " +
      std::to_string(n_particles));
}

// Returns weight * term with its external bra and ket indices relabeled by slot
// permutation p; bra and ket are permuted independently because triplet
// terms are not particle-symmetric
ExprPtr triplet_expr_for(const ExprPtr& term,
                         const container::svector<Index>& b,
                         const container::svector<Index>& k, std::size_t p,
                         const rational& w = 1) {
  const std::size_t n_particles = b.size();
  const auto ords = detail::compute_permuted_indices_bra_ket(p, n_particles);
  container::map<Index, Index> replacements;
  for (std::size_t i = 0; i != n_particles; ++i) {
    if (ords[i] != i) replacements.emplace(b[i], b[ords[i]]);
    if (ords[n_particles + i] != n_particles + i)
      replacements.emplace(k[i], k[ords[n_particles + i] - n_particles]);
  }
  if (replacements.empty())
    return w == 1 ? term->clone() : ex<Constant>(w) * term;
  return transform_expr(term, replacements, w);
}

/// Overlap matrix of the triplet state over the (n!)^2 slot permutations
Eigen::MatrixXd triplet_overlap_matrix(std::size_t n_particles) {
  const auto num_perms = static_cast<std::size_t>(factorial(n_particles));
  const auto dim = boost::numeric_cast<Eigen::Index>(num_perms * num_perms);

  // bra slot -> ket slot map of primitive idx
  auto pairing = [&](std::size_t idx) {
    perm::Permutation result = perm::unrank(idx / num_perms, n_particles);
    result->invert();
    result->postMultiply(perm::unrank(idx % num_perms, n_particles));
    return result;
  };
  auto triplet_bra = [&](std::size_t idx) {
    return perm::unrank(idx / num_perms, n_particles)->image(0);
  };

  Eigen::MatrixXd G(dim, dim);
  G.setZero();

  for (std::size_t row = 0; row < dim; ++row) {
    perm::Permutation ref = pairing(row);
    ref->invert();
    const auto row_triplet = triplet_bra(row);

    G(row, row) = std::pow(-2, n_particles);

    for (std::size_t col = row + 1; col < dim; ++col) {
      perm::Permutation current = pairing(col);
      current->postMultiply(ref);

      auto cycles = current->toDisjointCycles(n_particles);
      const auto col_triplet = triplet_bra(col);
      const bool same_cycle =
          std::any_of(cycles.begin(), cycles.end(), [&](const auto& cycle) {
            return std::find(cycle.begin(), cycle.end(), row_triplet) !=
                       cycle.end() &&
                   std::find(cycle.begin(), cycle.end(), col_triplet) !=
                       cycle.end();
          });
      if (!same_cycle) continue;

      std::size_t n_cycles = std::distance(cycles.begin(), cycles.end());
      auto entry = std::pow(-2, n_cycles);

      G(row, col) = entry;
      G(col, row) = entry;
    }
  }

  if (n_particles % 2 != 0) {
    G *= -1;
  }

  SEQUANT_ASSERT(G.isApprox(G.transpose()));

  return G;
}

/// rank of the triplet overlap matrix
sequant::intmax_t triplet_rank(std::size_t n) {
  return 3 * factorial(2 * n) / (factorial(n + 2) * factorial(n - 1));
}

/// last row of the pseudoinverse of the triplet overlap matrix; the hardcoded
/// triplet residual row differs from it by a null-space vector (for now)
std::vector<double> compute_triplet_biorthogonalizer_coeffs(
    std::size_t n_particles, double threshold) {
  auto overlap = triplet_overlap_matrix(n_particles);

  // Find Pseudo Inverse
  auto decomp = Eigen::CompleteOrthogonalDecomposition<decltype(overlap)>();
  decomp.setThreshold(threshold);
  decomp.compute(overlap);
  SEQUANT_ASSERT(decomp.rank() == triplet_rank(n_particles));

  Eigen::MatrixXd pinv = decomp.pseudoInverse();
  SEQUANT_ASSERT(pinv.isApprox(pinv.transpose()));

  return detail::flatten_row_major(pinv.row(pinv.rows() - 1));
}

ExprPtr triplet_combined_residual(
    const ExprPtr& V,
    const container::svector<container::svector<Index>>& ext_idxs) {
  const std::size_t n_particles = ext_idxs.size();
  const auto weights = hardcoded_triplet_residual_row(n_particles);

  // the n bra slots and the n ket slots are permuted independently
  SEQUANT_ASSERT(weights.size() ==
                 factorial(n_particles) * factorial(n_particles));

  const auto b = ext_idxs | ranges::views::transform([](const auto& group) {
                   return get_bra_idx(group);
                 }) |
                 ranges::to<container::svector<Index>>();
  const auto k = ext_idxs | ranges::views::transform([](const auto& group) {
                   return get_ket_idx(group);
                 }) |
                 ranges::to<container::svector<Index>>();
  const ExprPtr terms = V->is<Sum>() ? V : ex<Sum>(ExprPtrList{V});

  Sum out;
  for (const auto& term : *terms) {
    SEQUANT_ASSERT(term->is<Product>());
    if (!term->is<Product>()) {
      out.append(term);
      continue;
    }
    for (std::size_t p = 0; p != weights.size(); ++p) {
      const auto& w = weights[p];
      if (w == 0) continue;
      out.append(triplet_expr_for(term, b, k, p, w));
    }
  }

  auto result = ex<Sum>(out);
  simplify(result);
  return result;
}

ExprPtr triplet_maxcoeff_compact(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs) {
  const std::size_t n_particles = ext_idxs.size();
  if (n_particles <= 1) return expr;
  if (n_particles > 3)
    throw Exception(
        "triplet_maxcoeff_compact: only up to 3 external index groups are "
        "supported, got " +
        std::to_string(n_particles));
  if (!expr->is<Sum>()) return expr;

  auto work = expr->clone();
  for (auto& term : *work) {
    if (term->is<Product>())
      term =
          remove_tensor(term.as_shared_ptr<Product>(), reserved::symm_label());
  }
  canonicalize(work);
  simplify(work);

  // columns of the external groups other than group 0: V was signed by its
  // spin, so permutations moving it do not fix the TEE operator
  container::set<std::pair<Index, Index>> ext_columns;
  for (const auto& g : ext_idxs | ranges::views::drop(1))
    ext_columns.emplace(std::minmax(get_bra_idx(g), get_ket_idx(g)));

  // count the number of slot-permuted terms equal to this term
  auto stabilizer_order = [&](const Product& term) {
    sequant::intmax_t order = 1;
    for (const auto& f : term.factors()) {
      if (!f->is<Tensor>()) continue;
      const auto& t = f->as<Tensor>();
      if (t.column_symmetry() != ColumnSymmetry::Symm) continue;
      order *= factorial(static_cast<std::size_t>(ranges::count_if(
          ranges::views::zip(t.bra(), t.ket()), [&](const auto& col) {
            const auto& [x, y] = col;
            return ext_columns.contains(std::minmax(x, y));
          })));
    }
    return order;
  };

  // scale the compact term by rank/dim (doubles: rank 3 of 4, triples: 9 of 36)
  const rational dim_over_rank = n_particles == 2 ? rational(4, 3) : 4;

  Sum compact;
  for (const auto& [_, terms] : largest_coeff_terms_by_hash(work)) {
    // the max-|coeff| member is an identity-perm representative; its
    // coefficient also collects the images that map back onto it
    auto kept = terms.front()->clone();
    kept->as<Product>().scale(
        dim_over_rank * rational(1, stabilizer_order(kept->as<Product>())));
    compact.append(std::move(kept));
  }
  return ex<Sum>(compact);
}

namespace detail {

container::svector<size_t> compute_permuted_indices_bra_ket(
    size_t perm_index, size_t n_particles) {
  const auto num_perms = static_cast<size_t>(factorial(n_particles));
  SEQUANT_ASSERT(perm_index < num_perms * num_perms);

  container::svector<size_t> bra_slots(n_particles);
  container::svector<size_t> ket_slots(n_particles);
  for (size_t i = 0; i < n_particles; ++i) {
    bra_slots[i] = i;
    ket_slots[i] = n_particles + i;
  }

  container::svector<size_t> ords =
      compute_permuted_indices(bra_slots, perm_index / num_perms, n_particles);
  const auto permuted_ket =
      compute_permuted_indices(ket_slots, perm_index % num_perms, n_particles);
  ords.insert(ords.end(), permuted_ket.begin(), permuted_ket.end());
  return ords;
}

std::vector<double> compute_triplet_nns_p_coeffs(std::size_t n_particles,
                                                 double threshold) {
  const auto overlap = triplet_overlap_matrix(n_particles);
  const auto pinv_row =
      compute_triplet_biorthogonalizer_coeffs(n_particles, threshold);
  // G and pinv(G) are symmetric, so the last row of G pinv(G) is G times the
  // last row of pinv(G)
  return flatten_row_major(
      overlap *
      Eigen::Map<const Eigen::VectorXd>(
          pinv_row.data(), boost::numeric_cast<Eigen::Index>(pinv_row.size())));
}

std::vector<rational> compute_triplet_sparse_biorthogonalizer_coeffs(
    std::size_t n_particles, double threshold) {
  const auto overlap = triplet_overlap_matrix(n_particles);
  const auto dense_row =
      compute_triplet_biorthogonalizer_coeffs(n_particles, threshold);
  const auto dim = dense_row.size();
  const Eigen::Map<const Eigen::VectorXd> dense(
      dense_row.data(), boost::numeric_cast<Eigen::Index>(dim));
  // every valid row x satisfies overlap.x == target
  const Eigen::VectorXd target = overlap * dense;
  constexpr double rel_tol = 1e-8;
  const double tol = rel_tol * dense.cwiseAbs().maxCoeff();

  // column is the map key, so swapping two Es in TEE will be the same key.
  container::map<std::vector<double>, std::size_t> representative;
  for (std::size_t idx = 0; idx < dim; ++idx) {
    std::vector<double> column(dim);
    for (std::size_t r = 0; r < dim; ++r)
      column[r] =
          overlap(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(idx));
    representative[std::move(column)] =
        idx;  // 18 col for rank 3 (from 36 to 18 for rank 3)
  }

  // group the unique ones by equal dense coefficient (group by values)
  std::vector<std::size_t> by_value;
  for (const auto& [_, idx] : representative) by_value.push_back(idx);
  std::ranges::stable_sort(by_value, {},
                           [&](std::size_t i) { return dense_row[i]; });
  std::vector<std::vector<std::size_t>> groups;
  for (std::size_t k = 0; k < by_value.size(); ++k) {
    if (k == 0 || dense_row[by_value[k]] - dense_row[by_value[k - 1]] > tol)
      groups.emplace_back();
    groups.back().push_back(by_value[k]);
  }
  for (auto& g : groups) std::ranges::sort(g);
  std::ranges::sort(groups, {}, [](const auto& g) { return g.front(); });

  // 2^g = all, zero every subset of groups; keep the solvable one with the
  // higher sparsity
  std::vector<std::size_t> best_support;
  std::vector<double> best_values;
  for (std::uint64_t zeroed = 0; zeroed < (std::uint64_t{1} << groups.size());
       ++zeroed) {
    std::vector<std::size_t> cols;
    for (std::size_t g = 0; g < groups.size(); ++g)
      if (!((zeroed >> g) & 1u))
        cols.insert(cols.end(), groups[g].begin(), groups[g].end());
    if (cols.empty()) continue;
    std::ranges::sort(cols);

    Eigen::MatrixXd A(overlap.rows(), static_cast<Eigen::Index>(cols.size()));
    for (std::size_t j = 0; j < cols.size(); ++j)
      A.col(static_cast<Eigen::Index>(j)) =
          overlap.col(static_cast<Eigen::Index>(cols[j]));
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomp;
    decomp.setThreshold(threshold);
    decomp.compute(A);
    const Eigen::VectorXd x = decomp.solve(target);
    if ((A * x - target).norm() > rel_tol * target.norm())
      continue;  // if Gd=Gs

    std::vector<std::size_t> support;
    std::vector<double> values;
    for (std::size_t j = 0; j < cols.size(); ++j)
      if (const double v = x(static_cast<Eigen::Index>(j)); std::abs(v) > tol) {
        support.push_back(cols[j]);
        values.push_back(v);
      }
    if (best_support.empty() || support.size() < best_support.size() ||
        (support.size() == best_support.size() &&
         std::ranges::lexicographical_compare(best_support, support))) {
      best_support = std::move(support);
      best_values = std::move(values);
    }
  }
  SEQUANT_ASSERT(!best_support.empty());

  std::vector<rational> result(dim, rational{0});
  for (std::size_t j = 0; j < best_support.size(); ++j)
    result[best_support[j]] = to_rational(best_values[j], threshold);

  return result;
}

}  // namespace detail

}  // namespace sequant::mbpt
