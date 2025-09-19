//
// Created by Eduard Valeyev on 9/5/25.
//

#ifndef SEQUANT_CORE_TENSOR_NETWORK_TYPEDEFS_HPP
#define SEQUANT_CORE_TENSOR_NETWORK_TYPEDEFS_HPP

#include <SeQuant/core/index.hpp>

#include <boost/container/flat_set.hpp>

namespace sequant {

namespace tensor_network {
using NamedIndexSet = container::set<Index, Index::FullLabelCompare>;
using VertexColor = std::uint32_t;
}  // namespace tensor_network

}  // namespace sequant

#endif  // SEQUANT_CORE_TENSOR_NETWORK_TYPEDEFS_HPP
