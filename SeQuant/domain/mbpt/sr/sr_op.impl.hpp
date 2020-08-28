/// makes excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr T_(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false, bool antisymm = true);

/// makes lambda deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr Lambda_(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes excitation operator of all bra/ket ranks up to (and including) @c Nbra/Nket
ExprPtr T(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false, bool antisymm = true);

/// makes deexcitation operator of all bra/ket ranks up to (and including) @c Nbra/Nket
ExprPtr Lambda(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr A(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes L deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr L(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes R excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr R(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes geminal excitation operator from @p geminal_generating_space for ansatz @p ansatz
/// @param[in] geminal_generating_space the space from which the geminal excitations originate from; default = IndexSpace::active_occupied
/// @param[in] ansatz 1 or 2
ExprPtr R12(IndexSpace::Type geminal_generating_space = IndexSpace::active_occupied, int ansatz = 2);
