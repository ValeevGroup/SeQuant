#ifndef SEQUANT_CORE_OPTIMIZE_GA_GA_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_GA_HPP

// The search: an Engels-Putzka/Hanrath-style population loop whose mutation
// is k nearest-neighbour-interchange moves with k annealed over the
// generations, seeded at the per-term optimum (a GA from random genomes was
// measured 37,000x worse -- do not remove the seeding), with whole-block
// crossover and an NNI hill-climb polish.

#include <SeQuant/core/optimize/ga/fitness.hpp>

namespace sequant::opt::ga {

struct GAOptions {
  std::size_t population = 160;
  std::size_t generations = 150;
  std::size_t reproduction = 2;  ///< children per slot per generation
  double mixing = 0.20;          ///< non-elite fraction retained
  double kick0 = 10.0;           ///< mean NNI moves per mutation at gen 0
  double l1_move_prob = 0.85;    ///< NNI moves applied to L1 (vs L2) trees
  std::size_t restarts = 2;
  std::size_t hill_climb_sweeps = 40;
  std::uint64_t seed = 0;
};

/// Per-term optimal binarizations (an exact subset DP over the key table's
/// per-subset index masks -- see `seed_tree` in ga.cpp) plus flat L2 trees;
/// the mandatory starting point of the search.
///
/// Deliberately volatility-BLIND: the seed is the per-term flop optimum under
/// the unweighted cost, whatever weighting the `Fitness` the search then runs
/// on applies. That is the behavior the reference numbers were taken with.
Genome seed_genome(KeyTable const& kt);

/// Deterministic full-neighbourhood NNI descent on both layers; returns the
/// reached cost.
double hill_climb(Fitness const& fitness, Genome& genome,
                  std::size_t max_sweeps = 40);

/// Seed -> hill-climb -> `restarts` annealed GA runs -> polish.
/// Returns the best (cost, genome) found.
std::pair<double, Genome> run_ga(Fitness const& fitness, Genome seed,
                                 GAOptions const& opts = {});

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_GA_HPP
