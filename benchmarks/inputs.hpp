#ifndef SEQUANT_BENCHMARKS_INPUTS_HPP
#define SEQUANT_BENCHMARKS_INPUTS_HPP

#include <SeQuant/core/expr_fwd.hpp>

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstddef>

// Generators for realistic, large benchmark inputs. Everything is derived
// in-place so that no equation files need to be checked in. Results are
// memoized: the derivation cost is paid once, outside of any timed loop.

#ifdef SEQUANT_BENCH_MAX_CC_RANK
#define SEQUANT_BENCH_MAX_CC_RANK_DEFAULT SEQUANT_BENCH_MAX_CC_RANK
#else
#define SEQUANT_BENCH_MAX_CC_RANK_DEFAULT 10
#endif

#ifdef SEQUANT_BENCH_MAX_DOWNSTREAM_CC_RANK
#define SEQUANT_BENCH_MAX_DOWNSTREAM_CC_RANK_DEFAULT \
  SEQUANT_BENCH_MAX_DOWNSTREAM_CC_RANK
#else
#define SEQUANT_BENCH_MAX_DOWNSTREAM_CC_RANK_DEFAULT 3
#endif

namespace sequant::bench {

/// Highest CC rank for benchmarks of the derivation itself
inline constexpr std::size_t max_derivation_cc_rank =
    SEQUANT_BENCH_MAX_CC_RANK_DEFAULT;

/// Highest CC rank for benchmarks of processing steps that operate on the
/// derived residuals (simplification, spintracing, ...)
inline constexpr std::size_t max_cc_rank = std::min<std::size_t>(
    max_derivation_cc_rank, SEQUANT_BENCH_MAX_DOWNSTREAM_CC_RANK_DEFAULT);

/// @returns the spin-orbital CC T-amplitude residual of excitation rank
/// @p projection for a CC model of rank @p cc_rank, as it comes out of the
/// derivation (i.e. antisymmetrized tensors, simplified)
/// @note The returned expression is shared; clone before mutating it
const ExprPtr& cc_residual(std::size_t cc_rank, std::size_t projection);

/// @returns the expression that closed_shell_spintrace hands to simplify for
/// the residual (antisymmetrizer written in terms of the symmetrizer and all
/// antisymmetric tensors expanded), i.e. a large expression that has not been
/// simplified yet
/// @note The returned expression is shared; clone before mutating it
const ExprPtr& unsimplified_cc_residual(std::size_t cc_rank,
                                        std::size_t projection);

/// Registers one benchmark argument pair (cc_rank, projection) for every
/// residual of CC models of rank 1 to max_cc_rank. For use with
/// `benchmark::Benchmark::Apply`.
void cc_residual_args(benchmark::Benchmark* benchmark);

}  // namespace sequant::bench

#endif  // SEQUANT_BENCHMARKS_INPUTS_HPP
