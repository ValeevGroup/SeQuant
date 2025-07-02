#ifndef SEQUANT_CORE_EXPORT_EXPORT_NODE_HPP
#define SEQUANT_CORE_EXPORT_EXPORT_NODE_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/export/export_expr.hpp>

#include <type_traits>

namespace sequant {

template <typename T = ExportExpr>
  requires std::is_base_of_v<ExportExpr, T>
using ExportNode = FullBinaryNode<T>;

}

#endif  // SEQUANT_CORE_EXPORT_EXPORT_NODE_HPP
