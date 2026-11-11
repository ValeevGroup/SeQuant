#ifndef SEQUANT_CORE_EVAL_DISK_SPILL_HPP
#define SEQUANT_CORE_EVAL_DISK_SPILL_HPP

#include <cstddef>

namespace sequant::eval {

/// \brief Switch and accounting for charging large idle values of an ordered
///        schedule as kept on disk.
///
/// Wired onto the root cache with \c CacheManager::set_disk_spill. While it
/// is set, the ordered executor treats a value as *on disk* when the cell
/// holding it is one the schedule can keep there (\c disk_eligible_cells,
/// cell_table.hpp) and the value is at least \c min_bytes:
///   - it is not counted in the live working set, nor listed in a peak's
///     live set, while it is held;
///   - an operation that reads it is charged the result of the read (the
///     whole value, or the slice the read declares) as an operand, for as
///     long as the operation runs.
/// That is the cost of a value materialized from disk on every read. Nothing
/// is evaluated differently and no value is moved anywhere: this is
/// bookkeeping, for the dry-run meter and for a backend that does keep such
/// values on disk to be measured against.
///
/// Sizes are those the results report, so the accounting is meaningful only
/// where sizes are measured (see \c Result::size_in_bytes ).
struct DiskSpill {
  /// A value smaller than this is never on disk.
  std::size_t min_bytes = 1'000'000'000;

  /// \name accounting, filled by the executor; cumulative over evaluations
  /// until reset by the owner
  /// @{

  /// Bytes on disk when last sampled, and their high-water mark. Sampled
  /// with the working set, i.e. only while the eval trace is on.
  std::size_t bytes_on_disk = 0;
  std::size_t peak_bytes_on_disk = 0;
  /// Stores of an on-disk value: how many, and their bytes.
  std::size_t n_writes = 0;
  std::size_t bytes_written = 0;
  /// Reads of an on-disk value: how many, and the bytes of the values read.
  /// Each read is charged the whole value, an upper bound on what a read of
  /// a slice needs.
  std::size_t n_reads = 0;
  std::size_t bytes_read = 0;
  /// @}
};

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_DISK_SPILL_HPP
