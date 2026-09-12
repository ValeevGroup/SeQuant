//
// Created by Eduard Valeyev on 9/5/25.
//

#ifndef SEQUANT_CORE_TENSOR_NETWORK_TYPEDEFS_HPP
#define SEQUANT_CORE_TENSOR_NETWORK_TYPEDEFS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>

#include <boost/container/flat_set.hpp>

#include <cstddef>

namespace sequant {

namespace tensor_network {
using NamedIndexSet = IndexSet;
using VertexColor = std::uint32_t;

/// Per named (sliced) index loop-color: maps an Index to a small integer
/// identifying the DAG-scope loop that slices it. Passed to
/// TensorNetwork::canonicalize_slots / create_graph so that two same-space
/// named indices bound to DIFFERENT loops receive DIFFERENT graph colors (no
/// longer interchangeable), while same-loop indices stay interchangeable. An
/// empty (or null) map leaves canonicalization byte-identical to space-only
/// named coloring.
using NamedIndexColorMap =
    container::map<Index, std::size_t, Index::FullLabelCompare>;
}  // namespace tensor_network

}  // namespace sequant

#endif  // SEQUANT_CORE_TENSOR_NETWORK_TYPEDEFS_HPP
