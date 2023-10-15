/// makes excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr T_(std::size_t Nbra,
           std::size_t Nket = std::numeric_limits<std::size_t>::max());

/// makes lambda deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr Λ_(std::size_t Nbra,
           std::size_t Nket = std::numeric_limits<std::size_t>::max());

/// makes generic excitation (right-hand eigenvector) operator of bra/ket ranks
/// @c Nbra/Nket
ExprPtr R_(std::size_t Nbra,
           std::size_t Nket = std::numeric_limits<std::size_t>::max());

/// makes generic deexcitation (left-hand eigenvector) of bra/ket ranks @c
/// Nbra/Nket
ExprPtr L_(std::size_t Nbra,
           std::size_t Nket = std::numeric_limits<std::size_t>::max());

/// makes excitation operator of all bra/ket ranks up to (and including) @c
/// Nbra/Nket
ExprPtr T(std::size_t Nbra,
          std::size_t Nket = std::numeric_limits<std::size_t>::max());

/// makes deexcitation operator of all bra/ket ranks up to (and including) @c
/// Nbra/Nket
ExprPtr Λ(std::size_t Nbra,
          std::size_t Nket = std::numeric_limits<std::size_t>::max());

/// makes generic bra/ket-antisymmetrizer
/// \param Kh if <0, annihilates this many holes, else creates this many
/// \param Kp if <0, annihilates this many particles, else creates this many
/// (default is to set \p Kp to \p Kh)
ExprPtr A(std::int64_t Kh,
          std::int64_t Kp = std::numeric_limits<std::int64_t>::max());

/// makes generic particle-symmetrizer
/// \param Kh if <0, annihilates this many holes, else creates this many
/// \param Kp if <0, annihilates this many particles, else creates this many
/// (default is to set \p Kp to \p Kh)
ExprPtr S(std::int64_t Kh,
          std::int64_t Kp = std::numeric_limits<std::int64_t>::max());

/// makes L deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr L(std::size_t Nbra,
          std::size_t Nket = std::numeric_limits<std::size_t>::max());

/// makes R excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr R(std::size_t Nbra,
          std::size_t Nket = std::numeric_limits<std::size_t>::max());

/// makes geminal excitation operator from @p geminal_generating_space for
/// ansatz @p ansatz
/// @param[in] geminal_generating_space the space from which the geminal
/// excitations originate from; default = IndexSpace::active_occupied
/// @param[in] ansatz 1 or 2
ExprPtr R12(
    IndexSpace::Type geminal_generating_space = IndexSpace::active_occupied,
    int ansatz = 2);
