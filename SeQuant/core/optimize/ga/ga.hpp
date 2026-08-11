#ifndef SEQUANT_CORE_OPTIMIZE_GA_GA_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_GA_HPP

// The search: an Engels-Putzka/Hanrath-style population loop whose mutation
// is k nearest-neighbour-interchange moves with k annealed over the
// generations, seeded at the per-term optimum (a GA from random genomes was
// measured 37,000x worse -- do not remove the seeding), with whole-block
// crossover and an NNI hill-climb polish.

#include <SeQuant/core/optimize/ga/fitness.hpp>
#include <SeQuant/core/utility/aggregate.hpp>

namespace sequant::opt::ga {

struct GAOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;
  std::size_t population = 160;
  std::size_t generations = 150;
  std::size_t reproduction = 2;  ///< children per slot per generation
  double mixing = 0.20;          ///< non-elite fraction retained
  double kick0 = 10.0;           ///< mean NNI moves per mutation at gen 0
  double l1_move_prob = 0.85;    ///< NNI moves applied to L1 (vs L2) trees
  std::size_t restarts = 2;
  std::size_t hill_climb_sweeps = 40;
  std::uint64_t seed = 0;
  /// Decoded-tree memo capacity per evaluation workspace, in entries.
  std::size_t memo_capacity = TreeMemo::default_capacity;
  /// Print the search's per-phase wall time (hill_climb / breed / eval /
  /// select) to stderr.
  bool report_phases = false;
  /// Seed the L1 trees with the volatility-WEIGHTED per-term DP
  /// (`seed_genome(kt, cost)`) instead of the blind one. Off by default: the
  /// blind seed is what the recorded reference numbers assume, and the flag
  /// changes every vw > 1 schedule by construction.
  ///
  /// Why it exists: the blind seed is the vw = 1 optimum, and the population
  /// never crosses out of its basin -- MEASURED, its winners keep ~2x the
  /// per-iteration volatile flops of the weighted seed's at every vw > 1, and
  /// the whole GA-vs-GA gap at vw > 1 is that basin. With this on and nothing
  /// else changed, the search reaches 0.845-0.911 x the per-term schedule's
  /// own (cached-runtime) cost at vw = 20, where the blind seed reaches 1.72x.
  ///
  /// Inert unless the cost is volatility-aware; `optimize_ga` additionally
  /// gates it on a volatile-leaf predicate being supplied at all, so blind
  /// universes stay bit-identical whatever this says.
  bool volatility_aware_seed = false;
};

/// Per-term optimal binarizations (an exact subset DP over the key table's
/// per-subset index masks -- see `seed_tree` in ga.cpp) plus flat L2 trees;
/// the mandatory starting point of the search. Deliberately volatility-BLIND:
/// the per-term flop optimum under the unweighted cost, whatever weighting the
/// `Fitness` applies. That is the behavior the reference numbers assume, and
/// it is what `GAOptions::volatility_aware_seed == false` (the default) runs.
Genome seed_genome(KeyTable const& kt);

/// The same DP, same subset order, same bipartition order and same `<=`
/// tie-break -- one shared body -- with each merge charged \p cost instead of
/// the unweighted flop count: `cost.merge(T, a, b, is_volatile(a|b))`, the
/// PER-TERM replay rule (only amplitude-dependent work is charged
/// `volatile_weight` times; persistent work is assumed built once, which is
/// what a `CostModel::amortize_persistent` emission delivers). Reached only
/// under `GAOptions::volatility_aware_seed`. Reproduces `seed_genome(kt)` bit
/// for bit whenever the charge coincides -- in particular for
/// `CostModel::native()` on a volatility-blind table, and for any volatile
/// mask at `volatile_weight == 1`. L2 genes stay flat zeros either way.
Genome seed_genome(KeyTable const& kt, CostModel const& cost);

/// Deterministic full-neighbourhood NNI descent on both layers; returns the
/// reached cost.
double hill_climb(Fitness const& fitness, Genome& genome,
                  std::size_t max_sweeps = 40);

/// What the search saw on the way to its winner, for reporting. Every figure
/// is in the `Fitness` objective's own units, i.e. directly comparable with
/// the cost `run_ga` returns. Filled only when a trace is requested, since
/// `seed_cost` costs one extra fitness evaluation.
struct SearchTrace {
  /// Cost of the SEED genome, evaluated before any move is made: the per-term
  /// optimum's cost under this objective. This is the per-term-equivalent
  /// baseline -- the number the GA's winner has to be compared against for
  /// its advantage to mean anything.
  double seed_cost = 0;
  /// Cost after the pre-search hill-climb polish of the seed.
  double hill_climbed_cost = 0;
  /// Best cost reached by each restart, in restart order. Its spread is the
  /// only visible evidence that the search landscape is multi-basin (and that
  /// a given restart count is or is not enough).
  container::svector<double> restart_costs;
};

/// Seed -> hill-climb -> `restarts` annealed GA runs -> polish.
/// Returns the best (cost, genome) found. \p trace, when non-null, receives
/// the seed / per-restart costs.
std::pair<double, Genome> run_ga(Fitness const& fitness, Genome seed,
                                 GAOptions const& opts = {},
                                 SearchTrace* trace = nullptr);

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_GA_HPP
