#ifndef SEQUANT_CORE_EVAL_FLOPS_HPP
#define SEQUANT_CORE_EVAL_FLOPS_HPP

/// \file flops.hpp
/// Process-wide counter of the floating-point work an evaluation backend
/// ACTUALLY performs.
///
/// Motivation: the optimizer's cost model prices a contraction as a product of
/// nominal index extents. For block-sparse and tensor-of-tensor (CSV/PNO)
/// arrays that is an upper bound and nothing more -- per-pair inner dimensions
/// vary, and screened blocks are never formed. This counter is the measured
/// side of that comparison: every figure it records is derived from the real
/// range of a real tile (or, for a tensor of tensors, of a real inner tensor),
/// never from an index's nominal size. Comparing the two answers the question
/// the cost model cannot: does a schedule that is symbolically cheaper execute
/// fewer real flops?
///
/// Off by default. When disabled the only cost is one load of a
/// statically-addressed \c bool and a predictable branch at each backend
/// operation -- no atomics are touched and no tile is iterated.
///
/// Threading: evaluation is multithreaded, so accumulation is into a
/// thread-local bin (no synchronization on the hot path) and the bins are
/// reduced only when the totals are read.
///
/// MPI: every figure is LOCAL to the calling rank -- this file performs no
/// collective, so it can be called from a path that only some ranks reach.
/// Reduce with the caller's own collective at a point all ranks reach.

#include <array>
#include <cstddef>
#include <map>
#include <string>

namespace sequant::eval {

/// Kind of work a recorded backend operation performed.
enum class FlopCategory : std::size_t {
  Contraction = 0,  ///< binary product with at least one contracted index
  Addition,         ///< `a + b`, `a += b`
  Permute,          ///< permutation / copy / adjoint: data movement, 0 flops
  Scale,            ///< multiplication by a scalar
  Reduction,        ///< full dot product (rank-0 result)
};

inline constexpr std::size_t flop_category_count = 5;

/// Per-category tally.
struct FlopTally {
  double flops = 0;     ///< real floating-point operations performed
  double elements = 0;  ///< result elements actually written
  double ops = 0;       ///< number of backend operations recorded
};

/// A snapshot of the counter. All figures are local to the calling rank.
struct FlopReport {
  std::array<FlopTally, flop_category_count> by_category{};

  /// Number of elementary contraction kernel invocations. For a tensor of
  /// tensors this is one per (outer element, inner contraction), i.e. one per
  /// GEMM-shaped piece of work, NOT one per backend operation.
  double kernels = 0;
  /// Kernel-weighted sums of the three GEMM extents, so that `m_sum / kernels`
  /// is the mean external-left extent, and so on.
  double m_sum = 0, n_sum = 0, k_sum = 0;

  /// Operations whose real shape could not be determined from the data at
  /// hand (see \c real_flops_unsourced_ops). Their flops are NOT included
  /// above; a nonzero value here is a coverage hole and must be reported as
  /// such rather than ignored.
  double unsourced_ops = 0;

  [[nodiscard]] double flops() const noexcept;
  [[nodiscard]] double elements() const noexcept;
  [[nodiscard]] double ops() const noexcept;
  [[nodiscard]] FlopTally const& operator[](FlopCategory c) const noexcept {
    return by_category[static_cast<std::size_t>(c)];
  }
};

/// Process-wide real-flop counter.
///
/// \note \c enable() must be called before evaluation starts (it is not
///       synchronized against concurrent recording), and \c read() /
///       \c reset() at a point where no evaluation task is running -- e.g.
///       between CC iterations, after a fence.
class FlopCounter {
 public:
  /// Turn counting on or off. Off by default. Call once, before any
  /// evaluation.
  static void enable(bool on) noexcept { enabled_ = on; }

  /// The gate. This is the whole cost of the disabled path: one load of a
  /// statically-addressed bool.
  [[nodiscard]] static bool enabled() noexcept { return enabled_; }

  /// Zero every thread's bin.
  static void reset();

  /// Reduce the thread bins and return the totals for this rank.
  [[nodiscard]] static FlopReport read();

  /// Record one backend operation. Only called when \c enabled().
  static void record(FlopCategory category, double flops, double elements,
                     double ops = 1) noexcept;

  /// Record the GEMM shape of \p kernels elementary contractions. Separate
  /// from \c record() because one backend contraction over a tensor of
  /// tensors is many differently-shaped kernels.
  static void record_kernels(double kernels, double m_sum, double n_sum,
                             double k_sum) noexcept;

  /// Record an operation whose real shape could not be determined, tagged
  /// with why and with the annotation triple that produced it. Off the hot
  /// path by construction (it only runs when a shape could NOT be read), so
  /// it may take a lock and build a string.
  static void record_unsourced(std::string reason);

  /// Per-reason census of the unsourced operations, for reporting the
  /// coverage hole rather than hiding it. Local to this rank.
  [[nodiscard]] static std::map<std::string, double> unsourced_detail();

 private:
  // Plain bool, not an atomic: read on the hot path, written once during
  // setup. Making it atomic would put a (relaxed, but still) synchronized
  // load in front of every backend operation for no benefit.
  inline static bool enabled_ = false;
};

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_FLOPS_HPP
