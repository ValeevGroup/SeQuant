#ifndef SEQUANT_CORE_IRREPS_CONFIG_HPP
#define SEQUANT_CORE_IRREPS_CONFIG_HPP

#include <SeQuant/core/irreps/any.hpp>

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant {

/// Description of an IRREP context: the group (or rather its multiplication
/// rules) and the IRREPs that leaf tensors are assumed to transform as. See
/// make_irrep_context.
struct IrrepContextSpec {
  /// Names of IRREPs. The name "*" refers to any IRREP.
  using IrrepNames = std::vector<std::string>;

  /// An explicit multiplication table
  struct Table {
    /// Name of the totally symmetric IRREP
    std::string totally_symmetric;
    /// Products of the form "X*Y" together with the IRREPs they may contain
    std::vector<std::pair<std::string, IrrepNames>> products;
    /// Pairs of IRREPs and their conjugates
    std::vector<std::pair<std::string, std::string>> conjugates;
  };

  /// The name of the group: "symbolic" or one of the shipped point groups (C1,
  /// Ci, C2, Cs, C2v, C2h, D2, D2h). Unused if @c table is set.
  std::string group;
  /// If set, the group is given by this multiplication table
  std::optional<Table> table;
  /// Symbolic IRREPs. Only allowed for @c group (not for tables). If not set,
  /// symbols are only allowed for the "symbolic" group.
  std::optional<std::vector<std::string>> symbols;
  /// Relations of the form "lhs = rhs" between IRREPs. Only allowed for
  /// @c group (not for tables).
  std::optional<std::vector<std::string>> relations;
  /// IRREPs of leaf tensors. A key containing a '{' is a serialized tensor that
  /// is matched by block, all other keys are full-match regular expressions for
  /// the tensor label (an exact label takes precedence over patterns).
  std::vector<std::pair<std::string, IrrepNames>> tensors;
  /// IRREPs of leaf tensors that are not covered by @c tensors
  IrrepNames default_irreps = {"*"};
  /// Whether tensors that are not covered by @c tensors are an error (instead
  /// of having @c default_irreps)
  bool error_on_unspecified = false;
};

/// Creates an IRREP context from the given description
/// @throws Exception if the description is invalid
AnyIrrepContext make_irrep_context(const IrrepContextSpec &spec);

/// Removes the whitespace around the names in a textual IRREP specification
/// (a name, a product "X*Y" or a relation "X*Y = Z"), such that
/// "A * B = C" and "A*B=C" are equivalent. IRREP names must not contain
/// whitespace themselves.
/// @throws Exception if a name contains whitespace
std::string normalize_irrep_spec(std::string_view spec);

}  // namespace sequant

#endif  // SEQUANT_CORE_IRREPS_CONFIG_HPP
