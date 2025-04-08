#include <SeQuant/core/biorthogonalization.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/permutation.hpp>

#include <Eigen/Eigenvalues>

#include <libperm/Permutation.hpp>
#include <libperm/Rank.hpp>
#include <libperm/Utils.hpp>

#include <algorithm>

namespace sequant {

template <typename T>
struct compare_first_less {
  bool operator()(const T& lhs, const T& rhs) const {
    return lhs.first < rhs.first;
  }
};

using IndexPair = std::pair<Index, Index>;
using ParticlePairings = container::svector<IndexPair>;

ResultExpr biorthogonal_transform_copy(const ResultExpr& expr,
                                       double threshold) {
  container::svector<ResultExpr> wrapper = {expr.clone()};

  biorthogonal_transform(wrapper, threshold);

  return wrapper.front();
}

container::svector<ResultExpr> biorthogonal_transform_copy(
    const container::svector<ResultExpr>& exprs, double threshold) {
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

  // TODO: Can we fill the entire matrix only by knowing the entires of one
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

  assert(M.isApprox(M.transpose()));

  return M;
}

Eigen::MatrixXd compute_biorth_coeffs(std::size_t n_particles,
                                      double threshold) {
  auto perm_ovlp_mat = permutational_overlap_matrix(n_particles);
  assert(perm_ovlp_mat.rows() == perm_ovlp_mat.cols());
  assert(perm_ovlp_mat.isApprox(perm_ovlp_mat.transpose()));

  // Find Pseudo Inverse
  auto decomp =
      Eigen::CompleteOrthogonalDecomposition<decltype(perm_ovlp_mat)>();
  decomp.setThreshold(threshold);
  decomp.compute(perm_ovlp_mat);

  Eigen::MatrixXd pinv = decomp.pseudoInverse();
  // The pseudo inverse of a symmetric matrix should also be symmetric
  assert(pinv.isApprox(pinv.transpose()));

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
  assert(reference.size() == current.size());
  assert(std::is_sorted(reference.begin(), reference.end(),
                        compare_first_less<IndexPair>{}));
  assert(std::is_sorted(current.begin(), current.end(),
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
  assert(std::all_of(pairings.begin(), pairings.end(),
                     [](const ParticlePairings& pairing) {
                       return std::is_sorted(pairing.begin(), pairing.end(),
                                             compare_first_less<IndexPair>{});
                     }));

  container::set<std::pair<IndexSpace, IndexSpace>> ref_space_pairing;
  ref_space_pairing.reserve(ref_pairing.size());
  for (std::size_t i = 0; i < ref_pairing.size(); ++i) {
    ref_space_pairing.insert(
        std::make_pair(ref_pairing[i].first.space(),
                       ref_pairing[perm->image(i)].second.space()));
  }

  auto it = std::find_if(
      pairings.begin(), pairings.end(), [&](const ParticlePairings& p) {
        assert(p.size() == ref_pairing.size());

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
    throw std::runtime_error(
        "Missing explicit expression for a required index pairing in "
        "biorthogonalization");
  }

  auto idx = std::distance(pairings.begin(), it);
  const ParticlePairings& base = *it;

  assert(base.size() == ref_pairing.size());

  container::map<Index, Index> replacements;
  for (std::size_t i = 0; i < base.size(); ++i) {
    std::size_t ref_idx = perm->image(i);

    const bool differs_in_first = base[i].first != ref_pairing[i].first;
    const bool differs_in_second =
        base[i].second != ref_pairing[ref_idx].second;

    if (!differs_in_first && !differs_in_second) {
      // This particle pairing is identical
      continue;
    }

    assert(differs_in_first || differs_in_second);

    if (differs_in_first) {
      if (base[i].first.space() == ref_pairing[i].first.space()) {
        replacements.insert(
            std::make_pair(base[i].first, ref_pairing[ref_idx].first));
        continue;
      }
    }
    if (differs_in_second) {
      if (base[i].second.space() == ref_pairing[i].second.space()) {
        replacements.insert(
            std::make_pair(base[i].second, ref_pairing[ref_idx].second));
        continue;
      }
    }
  }

  ExprPtr expr = base_exprs.at(idx)->clone();

  if (!replacements.empty()) {
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
  assert(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.has_label() == result_exprs.front().has_label() &&
               (!expr.has_label() ||
                expr.label() == result_exprs.front().label());
      }));
  assert(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.symmetry() == result_exprs.front().symmetry();
      }));
  assert(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.braket_symmetry() == result_exprs.front().braket_symmetry();
      }));
  assert(std::all_of(result_exprs.begin(), result_exprs.end(),
                     [&](const ResultExpr& expr) {
                       return expr.particle_symmetry() ==
                              result_exprs.front().particle_symmetry();
                     }));
  assert(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.bra().size() == result_exprs.front().bra().size() &&
               std::is_permutation(expr.bra().begin(), expr.bra().end(),
                                   result_exprs.front().bra().begin());
      }));
  assert(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.ket().size() == result_exprs.front().ket().size() &&
               std::is_permutation(expr.ket().begin(), expr.ket().end(),
                                   result_exprs.front().ket().begin());
      }));
  assert(std::all_of(
      result_exprs.begin(), result_exprs.end(), [&](const ResultExpr& expr) {
        return expr.aux().size() == result_exprs.front().aux().size() &&
               std::is_permutation(expr.aux().begin(), expr.aux().end(),
                                   result_exprs.front().aux().begin());
      }));
  assert(std::all_of(result_exprs.begin(), result_exprs.end(),
                     [](const ResultExpr& res) {
                       return res.particle_symmetry() == ParticleSymmetry::symm;
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

  Eigen::MatrixXd coefficients = compute_biorth_coeffs(n_particles, threshold);

  auto num_perms = factorial(n_particles);
  assert(num_perms == coefficients.rows());
  assert(num_perms == coefficients.cols());

  auto original_exprs = result_exprs |
                        ranges::views::transform([](const ResultExpr& res) {
                          return res.expression();
                        }) |
                        ranges::to<container::svector<ExprPtr>>();

  for (std::size_t i = 0; i < result_exprs.size(); ++i) {
    result_exprs.at(i).expression() = ex<Constant>(0);
    perm::Permutation reference = perm::unrank(ranks.at(i), n_particles);
    reference->invert();

    for (std::size_t rank = 0; rank < num_perms; ++rank) {
      perm::Permutation perm = perm::unrank(rank, n_particles);
      perm->postMultiply(reference);

      result_exprs.at(i).expression() +=
          ex<Constant>(
              to_rational(coefficients(ranks.at(i), rank), threshold)) *
          create_expr_for(externals.at(i), perm, externals, original_exprs);
    }

    simplify(result_exprs.at(i).expression());
  }
}

ExprPtr biorthogonal_transform(
    const sequant::ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups,
    const double threshold) {
  ResultExpr res(
      bra(ext_index_groups | ranges::views::transform([](const auto& pair) {
            return pair.at(0);
          }) |
          ranges::to<container::svector<Index>>()),
      ket(ext_index_groups | ranges::views::transform([](const auto& pair) {
            return pair.at(1);
          }) |
          ranges::to<container::svector<Index>>()),
      aux(IndexList{}), Symmetry::nonsymm, BraKetSymmetry::nonsymm,
      ParticleSymmetry::symm, {}, expr);

  biorthogonal_transform(res, threshold);

  return res.expression();
}

}  // namespace sequant
