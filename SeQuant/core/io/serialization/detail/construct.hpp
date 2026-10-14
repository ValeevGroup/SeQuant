#ifndef SEQUANT_CORE_IO_SERIALIZATION_DETAIL_CONSTRUCT_HPP
#define SEQUANT_CORE_IO_SERIALIZATION_DETAIL_CONSTRUCT_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>

#include <string>
#include <variant>

namespace sequant::io::serialization::detail {

/// The (bra/ket permutational, bra<->ket exchange, particle/column) symmetry
/// triple used while deserializing; the second element is either a concrete
/// BraKetSymmetry or an abstract Hermiticity (resolved against the tensor's
/// field at construction time).
using DefaultSymmetries =
    std::tuple<Symmetry, std::variant<BraKetSymmetry, Hermiticity>,
               ColumnSymmetry>;

/// Resolves the symmetry values to assume for input that does not spell them
/// out: the active Context's deserialization defaults, overridden by whatever
/// \p options provides. (Anti)symmetry in bra and ket implies particle
/// symmetry, so a non-Nonsymm permutation symmetry forces column symmetry.
DefaultSymmetries resolve_default_symmetries(
    const DeserializationOptions &options);

/// @return the vacuum a tilde-labelled NormalOperator refers to:
/// NormalOperator::label() is the tilde label for every non-Physical vacuum
Vacuum tilde_vacuum(Statistics s);

/// Builds a Tensor, applying the symmetry rules that are fixed by a reserved
/// label (the (anti)symmetrizer, the overlap/Kronecker metric and the reference
/// densities) so that a deserialized tensor hashes identically to one built by
/// the corresponding factory.
///
/// \p perm, \p braket and \p column are the symmetry values already resolved
/// against the defaults; the `*_specified` flags record whether the input
/// spelled each one out (only a spelled-out value may contradict a symmetry
/// fixed by a reserved label -- a defaulted one is replaced silently).
ExprPtr make_tensor(const std::wstring &name, container::vector<Index> bra,
                    container::vector<Index> ket, container::vector<Index> aux,
                    Symmetry perm,
                    std::variant<BraKetSymmetry, Hermiticity> braket,
                    ColumnSymmetry column, bool perm_specified,
                    bool braket_specified, bool column_specified);

}  // namespace sequant::io::serialization::detail

#endif  // SEQUANT_CORE_IO_SERIALIZATION_DETAIL_CONSTRUCT_HPP
