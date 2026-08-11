#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/detail/concepts.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/container.hpp>
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

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <iostream>
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
      return std::vector<sequant::rational>{ratio(17, 120), ratio(-7, 120),
                                            ratio(-1, 120), ratio(-1, 120),
                                            ratio(-1, 120), ratio(-7, 120)};

    case 4:
      return std::vector<sequant::rational>{
          ratio(43, 840), ratio(-19, 1680), ratio(-19, 1680),
          ratio(-1, 105), ratio(-19, 1680), ratio(-19, 1680),
          ratio(13, 840), ratio(1, 120),    ratio(-1, 105),
          ratio(1, 120),  ratio(-1, 105),   ratio(-19, 1680),
          ratio(-1, 105), ratio(1, 120),    ratio(1, 120),
          ratio(13, 840), ratio(-1, 105),   ratio(-1, 105),
          ratio(1, 120),  ratio(-19, 1680), ratio(-19, 1680),
          ratio(13, 840), ratio(-19, 1680), ratio(1, 120)};

    case 5:
      return std::vector<sequant::rational>{
          ratio(59, 3780),   ratio(-5, 3024),   ratio(-5, 3024),
          ratio(-5, 3024),   ratio(-31, 7560),  ratio(-5, 3024),
          ratio(-5, 3024),   ratio(-23, 30240), ratio(19, 7560),
          ratio(37, 15120),  ratio(-5, 3024),   ratio(-23, 30240),
          ratio(-5, 3024),   ratio(19, 7560),   ratio(37, 15120),
          ratio(-31, 7560),  ratio(37, 15120),  ratio(37, 15120),
          ratio(-31, 7560),  ratio(-5, 3024),   ratio(-5, 3024),
          ratio(-23, 30240), ratio(-23, 30240), ratio(-23, 30240),
          ratio(-13, 7560),  ratio(-5, 3024),   ratio(-5, 3024),
          ratio(19, 7560),   ratio(-23, 30240), ratio(37, 15120),
          ratio(19, 7560),   ratio(-23, 30240), ratio(19, 7560),
          ratio(-23, 30240), ratio(-13, 7560),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(-23, 30240), ratio(-31, 7560),  ratio(-13, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(19, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-13, 7560),
          ratio(-13, 7560),  ratio(-23, 30240), ratio(-31, 7560),
          ratio(37, 15120),  ratio(-31, 7560),  ratio(37, 15120),
          ratio(-5, 3024),   ratio(-5, 3024),   ratio(-23, 30240),
          ratio(19, 7560),   ratio(-5, 3024),   ratio(37, 15120),
          ratio(-31, 7560),  ratio(37, 15120),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(19, 7560),   ratio(37, 15120),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(-13, 7560),
          ratio(-23, 30240), ratio(37, 15120),  ratio(-13, 7560),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(-23, 30240),
          ratio(19, 7560),   ratio(-23, 30240), ratio(-23, 30240),
          ratio(19, 7560),   ratio(-13, 7560),  ratio(-31, 7560),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(19, 7560),   ratio(-31, 7560),  ratio(-31, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-5, 3024),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(-23, 30240), ratio(-5, 3024),
          ratio(19, 7560),   ratio(-23, 30240), ratio(-5, 3024),
          ratio(37, 15120),  ratio(-5, 3024),   ratio(-23, 30240),
          ratio(-23, 30240), ratio(-23, 30240), ratio(-13, 7560),
          ratio(19, 7560),   ratio(19, 7560),   ratio(-23, 30240),
          ratio(-23, 30240), ratio(-13, 7560),  ratio(-5, 3024),
          ratio(19, 7560),   ratio(-5, 3024),   ratio(-23, 30240),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-13, 7560),
          ratio(-13, 7560),  ratio(37, 15120),  ratio(-23, 30240)};

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

  for (std::size_t row = 0; row < n; ++row) {
    for (std::size_t col = 0; col < n; ++col) {
      perm::Permutation row_perm = perm::unrank(n - 1 - row, n_particles);
      perm::Permutation col_perm = perm::unrank(col, n_particles);

      col_perm->preMultiply(row_perm);

      std::size_t source_idx = perm::rank(col_perm, n_particles);
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

  Eigen::MatrixXd pinv = decomp.pseudoInverse();
  // The pseudo inverse of a symmetric matrix should also be symmetric
  SEQUANT_ASSERT(pinv.isApprox(pinv.transpose()));

  // We need to normalize to the amount of non-zero eigenvalues via
  // normalization = #eigenvalues / #non-zero eigenvalues
  // Since perm_ovlp_mat is symmetric, it is diagonalizable and for every
  // diagonalizable matrix, its rank equals the amount of non-zero eigenvalues.
  double normalization =
      static_cast<double>(perm_ovlp_mat.rows()) / decomp.rank();

  pinv *= normalization;

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

  for (std::size_t i = 0; i < result_exprs.size(); ++i) {
    result_exprs.at(i).expression() = ex<Constant>(0);
    perm::Permutation reference = perm::unrank(ranks.at(i), n_particles);
    reference->invert();

    for (std::size_t rank = 0; rank < num_perms; ++rank) {
      perm::Permutation perm = perm::unrank(rank, n_particles);
      perm->postMultiply(reference);

      sequant::rational coeff =
          (n_particles <= max_rank_hardcoded_biorthogonalizer_matrix)
              ? (*hardcoded_coefficients)(ranks.at(i), rank)
              : to_rational((*computed_coefficients)(ranks.at(i), rank),
                            threshold);

      result_exprs.at(i).expression() +=
          ex<Constant>(coeff) *
          create_expr_for(externals.at(i), perm, externals, original_exprs);
    }

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

template <detail::index_group_range IdxGroups>
ExprPtr WK_biorthogonalization_filter_impl(ExprPtr expr, IdxGroups&& ext_idxs) {
  if (!expr->is<Sum>()) return expr;
  if (ext_idxs.size() <= 2) return expr;  // always skip R1 and R2

  // hash filtering logic for R > 2
  container::map<std::size_t, container::vector<ExprPtr>> largest_coeff_terms;

  for (const auto& term : *expr) {
    if (!term->is<Product>()) continue;

    auto product = term.as_shared_ptr<Product>();
    auto scalar = product->scalar();

    sequant::TensorNetwork tn(*product);
    auto hash =
        tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels())
            .hash_value();

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

  Sum filtered;
  for (const auto& [_, terms] : largest_coeff_terms) {
    for (const auto& t : terms) {
      filtered.append(t);
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

namespace {

std::size_t product_network_hash(const ExprPtr& term) {
  SEQUANT_ASSERT(term->is<Product>());
  auto product = term.as_shared_ptr<Product>();
  sequant::TensorNetwork tn(*product);
  return tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels())
      .hash_value();
}

container::svector<std::string> split_annotation(std::string const& annot) {
  container::svector<std::string> parts;
  std::string part;
  for (char c : annot) {
    if (c == ',') {
      parts.push_back(part);
      part.clear();
    } else if (c != ' ') {
      part.push_back(c);
    }
  }
  if (!part.empty()) parts.push_back(part);
  return parts;
}

std::string join_annotation(container::svector<std::string> const& parts) {
  std::string out;
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (i) out.push_back(',');
    out += parts[i];
  }
  return out;
}

/// The triplet n = 3 non-null-space weights
rational triplet_triples_slot_weight(
    const container::svector<std::size_t>& ords) {
  // w10 is the first row of the non-null-sapce projector: G . pinv(G)
  // of the 18x18 TEE overlap matrix
  constexpr std::array<int, 18> w10{5, -1, -1, -1, 1, 1, -2, -2, 1,
                                    0, 1,  0,  0,  0, 1, -2, 1,  -2};

  SEQUANT_ASSERT(ords.size() == 6);
  constexpr std::size_t n = 3;

  std::array<std::size_t, n> ph{}, pv{}, pv_inv{};
  for (std::size_t i = 0; i != n; ++i) {
    ph[i] = ords[i];
    pv[i] = ords[n + i] - n;
  }
  for (std::size_t i = 0; i != n; ++i) pv_inv[pv[i]] = i;

  std::array<std::size_t, n> sigma{};
  for (std::size_t i = 0; i != n; ++i) sigma[i] = ph[pv_inv[i]];

  // lexicographic rank of sigma among the n! permutations of {0, .., n-1}
  std::array<std::size_t, n> perm{};
  std::iota(perm.begin(), perm.end(), std::size_t{0});
  std::size_t sigma_rank = 0;
  while (perm != sigma) {
    ++sigma_rank;
    [[maybe_unused]] const bool has_next =
        std::next_permutation(perm.begin(), perm.end());
    SEQUANT_ASSERT(has_next);
  }

  return ratio(w10[n * sigma_rank + pv[0]], 20);
}

// (n!)^2: the n bra slots and the n ket slots are permuted
// independently
std::size_t slot_perm_count(std::size_t n_particles) {
  return static_cast<std::size_t>(factorial(n_particles) *
                                  factorial(n_particles));
}

// flat index of a slot permutation
std::size_t slot_perm_index(const container::svector<std::size_t>& ords) {
  SEQUANT_ASSERT(ords.size() % 2 == 0);
  const std::size_t n_particles = ords.size() / 2;
  const auto num_perms = slot_perm_count(n_particles);
  for (std::size_t p = 0; p != num_perms; ++p) {
    if (detail::compute_bra_ket_permuted_indices(p, n_particles) == ords)
      return p;
  }
  throw Exception("slot_perm_index: not a valid S_n x S_n slot permutation");
}

// slot perms: {3/4, -1/4, -1/4, -1/4} for n = 2;
// for n = 3 the per-representative weights w10[m]/20 distributed over the 36
// slot perms by 18x2
std::vector<rational> make_triplet_nullspace_weights(std::size_t n_particles) {
  switch (n_particles) {
    case 2: {
      std::vector<rational> weights(slot_perm_count(2), ratio(-1, 4));
      container::svector<std::size_t> identity_ords(4);
      std::iota(identity_ords.begin(), identity_ords.end(), std::size_t{0});
      weights[slot_perm_index(identity_ords)] = ratio(3, 4);
      return weights;
    }

    case 3: {
      std::vector<rational> weights(slot_perm_count(3), 0);
      for (std::size_t p = 0; p != weights.size(); ++p) {
        weights[p] = triplet_triples_slot_weight(
            detail::compute_bra_ket_permuted_indices(p, 3));
      }
      return weights;
    }

    default:
      throw Exception(
          "triplet slot-perm weights only available for n_particles = "
          "2, 3, requested rank is : " +
          std::to_string(n_particles));
  }
}

// the paper-combined residual assembly weights; see triplet_combined_residual
// for the formulas
std::vector<rational> make_triplet_combined_residual_weights(
    std::size_t n_particles, bool te_only) {
  std::vector<rational> weights(slot_perm_count(n_particles), 0);

  auto set = [&](container::svector<std::size_t> bra,
                 const container::svector<std::size_t>& ket, rational coeff) {
    SEQUANT_ASSERT(bra.size() == n_particles && ket.size() == n_particles);
    for (const auto s : ket) bra.push_back(n_particles + s);
    weights[slot_perm_index(bra)] = coeff;
  };

  if (te_only) {
    if (n_particles != 2)
      throw Exception(
          "the bare-TE triplet residual weights are only defined for "
          "n_particles = 2");
    set({0, 1}, {0, 1}, ratio(1, 4));
    return weights;
  }

  switch (n_particles) {
    case 2:
      set({0, 1}, {0, 1}, ratio(3, 16));   // identity
      set({1, 0}, {1, 0}, ratio(-1, 16));  // whole-pair swap
      break;

    case 3:
      set({0, 1, 2}, {0, 1, 2}, ratio(6, 160));   // identity
      set({1, 0, 2}, {1, 0, 2}, ratio(-1, 160));  // pair swap 0 <-> 1
      set({2, 1, 0}, {2, 1, 0}, ratio(-1, 160));  // pair swap 0 <-> 2
      set({0, 1, 2}, {0, 2, 1}, ratio(2, 160));   // ket swap 1 <-> 2
      break;

    default:
      throw Exception(
          "triplet paper-combined residual weights are only available for "
          "n_particles = 2, 3, requested rank is : " +
          std::to_string(n_particles));
  }
  return weights;
}

// (memoized) rational weight rows for the symbolic triplet primitives
const std::vector<rational>& triplet_weights_rational(std::size_t n_particles,
                                                      TripletWeightKind kind) {
  using CacheKey = std::pair<std::size_t, TripletWeightKind>;
  // the rows are cached behind a pointer since the cache holds several keys
  // and flat_map insertions would invalidate references into it
  using CachedRow = std::unique_ptr<const std::vector<rational>>;

  static std::mutex cache_mutex;
  static std::condition_variable cache_cv;
  static container::map<CacheKey, std::optional<CachedRow>> cache;

  CacheKey key{n_particles, kind};

  return *sequant::detail::memoize(
      cache, cache_mutex, cache_cv, key, [&]() -> CachedRow {
        auto make_row = [&]() -> std::vector<rational> {
          switch (kind) {
            case TripletWeightKind::NullspaceProjector:
              return make_triplet_nullspace_weights(n_particles);

            case TripletWeightKind::NnsReconstruction: {
              auto weights = make_triplet_nullspace_weights(n_particles);
              container::svector<std::size_t> identity_ords(2 * n_particles);
              std::iota(identity_ords.begin(), identity_ords.end(),
                        std::size_t{0});
              const auto identity_weight =
                  weights.at(slot_perm_index(identity_ords));
              SEQUANT_ASSERT(identity_weight != 0);
              for (auto& w : weights) w /= identity_weight;
              return weights;
            }

            case TripletWeightKind::CombinedResidual:
              return make_triplet_combined_residual_weights(n_particles,
                                                            /*te_only=*/false);

            case TripletWeightKind::TeCombinedResidual:
              return make_triplet_combined_residual_weights(n_particles,
                                                            /*te_only=*/true);

            case TripletWeightKind::TeNnsReconstruction:
            case TripletWeightKind::TeReconstruction: {
              if (n_particles != 2)
                throw Exception(
                    "bare-TE triplet weights are only defined for "
                    "n_particles = 2");
              const auto swap_weight =
                  kind == TripletWeightKind::TeNnsReconstruction ? ratio(-1, 2)
                                                                 : ratio(1, 4);
              std::vector<rational> weights(slot_perm_count(2), 0);
              for (std::size_t p = 0; p != weights.size(); ++p) {
                const auto ords =
                    detail::compute_bra_ket_permuted_indices(p, 2);
                const bool bra_swapped = ords[0] != 0;
                const bool ket_swapped = ords[2] != 2;
                if (!bra_swapped && !ket_swapped)
                  weights[p] = 1;
                else if (bra_swapped != ket_swapped)
                  weights[p] = swap_weight;
                // the pair swap carries weight 0 in the bare-TE rows
              }
              return weights;
            }
          }
          SEQUANT_UNREACHABLE;
        };
        return std::make_unique<const std::vector<rational>>(make_row());
      });
}

// index replacement map realizing slot permutation p on the external bra/ket
// index n-tuples: b[i] -> b[ords[i]], k[i] -> k[ords[n + i] - n]
container::map<Index, Index> slot_perm_replacements(
    const container::svector<Index>& b, const container::svector<Index>& k,
    std::size_t p) {
  const std::size_t n_particles = b.size();
  const auto ords = detail::compute_bra_ket_permuted_indices(p, n_particles);
  container::map<Index, Index> m;
  for (std::size_t i = 0; i != n_particles; ++i) {
    if (ords[i] != i) m.emplace(b[i], b[ords[i]]);
    if (ords[n_particles + i] != n_particles + i)
      m.emplace(k[i], k[ords[n_particles + i] - n_particles]);
  }
  return m;
}

std::pair<container::svector<Index>, container::svector<Index>>
external_bra_ket(
    const container::svector<container::svector<Index>>& ext_idxs) {
  container::svector<Index> b, k;
  b.reserve(ext_idxs.size());
  k.reserve(ext_idxs.size());
  for (const auto& group : ext_idxs) {
    b.push_back(get_bra_idx(group));
    k.push_back(get_ket_idx(group));
  }
  return {std::move(b), std::move(k)};
}

// the sign the canonicalization produced is returned separately so callers
// can track coefficients
std::pair<ExprPtr, Product::scalar_type> canonical_unit_product(
    const ExprPtr& t) {
  auto out = t->clone();
  out->as<Product>().scale(Product::scalar_type{1} /
                           out->as<Product>().scalar());
  canonicalize(out);
  simplify(out);
  const auto sign = out->as<Product>().scalar();
  out->as<Product>().scale(Product::scalar_type{1} / sign);
  return {std::move(out), sign};
}

ExprPtr triplet_weighted_perm_sum(
    const ExprPtr& compact_expr,
    const container::svector<container::svector<Index>>& ext_idxs,
    const std::vector<rational>& weights) {
  const auto [b, k] = external_bra_ket(ext_idxs);

  Sum out;
  for (const auto& term : *compact_expr) {
    if (!term->is<Product>()) {
      out.append(term);
      continue;
    }
    for (std::size_t p = 0; p != weights.size(); ++p) {
      const auto& w = weights[p];
      if (w == 0) continue;
      const auto map = slot_perm_replacements(b, k, p);
      if (map.empty())
        out.append(w == 1 ? term : ex<Constant>(w) * term);
      else
        out.append(transform_expr(term, map, w));
    }
  }

  auto result = ex<Sum>(out);
  simplify(result);
  return result;
}

}  // namespace

ExprPtr triplet_combined_residual(
    const ExprPtr& V,
    const container::svector<container::svector<Index>>& ext_idxs,
    bool te_only) {
  const std::size_t n_particles = ext_idxs.size();
  const auto& weights = triplet_weights_rational(
      n_particles, te_only ? TripletWeightKind::TeCombinedResidual
                           : TripletWeightKind::CombinedResidual);
  SEQUANT_ASSERT(weights.size() == slot_perm_count(n_particles));

  if (V->is<Sum>()) return triplet_weighted_perm_sum(V, ext_idxs, weights);

  const auto [b, k] = external_bra_ket(ext_idxs);
  Sum out;
  for (std::size_t p = 0; p != weights.size(); ++p) {
    const auto& w = weights[p];
    if (w == 0) continue;
    const auto map = slot_perm_replacements(b, k, p);
    out.append(map.empty() ? ex<Constant>(w) * V : transform_expr(V, map, w));
  }
  auto result = ex<Sum>(out);
  simplify(result);
  return result;
}

ExprPtr triplet_maxcoeff_compact(
    ExprPtr expr, const container::svector<container::svector<Index>>& ext_idxs,
    TripletWeightKind kind) {
  if (!expr->is<Sum>()) return expr;
  const std::size_t n_particles = ext_idxs.size();
  if (n_particles != 2 && n_particles != 3) return expr;

  auto work = expr->clone();
  for (auto& term : *work) {
    if (term->is<Product>())
      term =
          remove_tensor(term.as_shared_ptr<Product>(), reserved::symm_label());
  }
  canonicalize(work);
  simplify(work);

  const auto& weights = triplet_weights_rational(n_particles, kind);
  const auto [b, k] = external_bra_ket(ext_idxs);

  container::map<std::size_t, container::vector<ExprPtr>> groups;
  for (const auto& term : *work) {
    if (!term->is<Product>()) continue;
    groups[product_network_hash(term)].push_back(term);
  }

  Sum compact;
  for (const auto& [_, terms] : groups) {
    // the max-|coeff| member is an identity-perm representative
    const ExprPtr* rep = &terms.front();
    for (const auto& t : terms)
      if (abs(t->as<Product>().scalar()) > abs((*rep)->as<Product>().scalar()))
        rep = &t;
    const auto [rep_unit, rep_sign] = canonical_unit_product(*rep);
    const auto rep_coeff = (*rep)->as<Product>().scalar() * rep_sign;

    container::map<std::size_t, Product::scalar_type> pred_weight;
    for (std::size_t p = 0; p != weights.size(); ++p) {
      const auto map = slot_perm_replacements(b, k, p);
      ExprPtr t =
          map.empty() ? rep_unit->clone() : transform_expr(rep_unit, map);
      auto [u, s] = canonical_unit_product(t);
      const auto w = Product::scalar_type(weights[p]) * s;
      if (auto it = pred_weight.find(u->hash_value()); it != pred_weight.end())
        it->second += w;
      else
        pred_weight.emplace(u->hash_value(), w);
    }

    const auto rep_weight = pred_weight.at(rep_unit->hash_value());
    SEQUANT_ASSERT(rep_weight != Product::scalar_type(0));
    const auto kept_coeff = rep_coeff / rep_weight;

    std::size_t n_predicted_members = 0;
    for (const auto& [h, w] : pred_weight)
      if (w != Product::scalar_type(0)) ++n_predicted_members;
    if (n_predicted_members != terms.size())
      throw Exception(
          "triplet_maxcoeff_compact: hash group does not match the "
          "single-representative pattern; the residual cannot be "
          "compacted losslessly");
    for (const auto& t : terms) {
      const auto [u, s] = canonical_unit_product(t);
      const auto it = pred_weight.find(u->hash_value());
      if (it == pred_weight.end() ||
          t->as<Product>().scalar() * s != kept_coeff * it->second)
        throw Exception(
            "triplet_maxcoeff_compact: hash group does not match the "
            "single-representative pattern; the residual cannot be "
            "compacted losslessly");
    }

    auto kept = rep_unit->clone();
    kept->as<Product>().scale(kept_coeff);
    compact.append(std::move(kept));
  }
  return ex<Sum>(compact);
}

ExprPtr triplet_symbolic_reconstruct(
    ExprPtr compact_expr,
    const container::svector<container::svector<Index>>& ext_idxs,
    TripletWeightKind kind) {
  if (!compact_expr->is<Sum>()) return compact_expr;
  const std::size_t n_particles = ext_idxs.size();
  if (n_particles != 2 && n_particles != 3) return compact_expr;
  return triplet_weighted_perm_sum(compact_expr, ext_idxs,
                                   triplet_weights_rational(n_particles, kind));
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
  ExprPtr S_tensor =
      ex<Tensor>(Tensor{reserved::symm_label(), bra(kixs), ket(bixs)});

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
  auto normalized_pinv =
      compute_biorthogonalizer_matrix(n_particles, threshold);
  Eigen::MatrixXd nns_matrix = perm_ovlp_mat * normalized_pinv;

  auto num_perms = nns_matrix.rows();
  std::vector<double> coeffs;
  coeffs.reserve(num_perms);
  for (std::size_t i = 0; i < num_perms; ++i) {
    coeffs.push_back(nns_matrix(num_perms - 1, i));
  }
  return coeffs;
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

container::svector<size_t> compute_bra_ket_permuted_indices(
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

container::svector<std::string> slot_perm_annots(std::string const& orig_annot,
                                                 size_t n_particles) {
  const auto parts = split_annotation(orig_annot);
  SEQUANT_ASSERT(parts.size() == 2 * n_particles &&
                 "slot_perm_annots: annotation rank must equal "
                 "2 * n_particles");

  const auto num_perms = slot_perm_count(n_particles);
  container::svector<std::string> annots;
  annots.reserve(num_perms);
  for (size_t p = 0; p != num_perms; ++p) {
    const auto ords = compute_bra_ket_permuted_indices(p, n_particles);
    container::svector<std::string> permuted(parts.size());
    for (size_t s = 0; s != parts.size(); ++s) permuted[s] = parts[ords[s]];
    annots.push_back(join_annotation(permuted));
  }
  return annots;
}

std::vector<double> compute_triplet_weights(std::size_t n_particles,
                                            TripletWeightKind kind) {
  const auto& weights = triplet_weights_rational(n_particles, kind);
  std::vector<double> coeffs;
  coeffs.reserve(weights.size());
  for (const auto& w : weights) {
    coeffs.push_back(static_cast<double>(w));
  }
  return coeffs;
}

}  // namespace detail

}  // namespace sequant::mbpt
