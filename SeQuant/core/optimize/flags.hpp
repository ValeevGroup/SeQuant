#ifndef SEQUANT_CORE_OPTIMIZE_FLAGS_HPP
#define SEQUANT_CORE_OPTIMIZE_FLAGS_HPP

namespace sequant {

/// Cost metric to optimize for in single-term and top-level optimize routines.
enum class OptFor { Flops, Memsize };

/// Whether to reorder summands so terms with shared intermediates appear
/// closer to each other.
enum class ReorderSum { Reorder, NoReorder };

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIMIZE_FLAGS_HPP
