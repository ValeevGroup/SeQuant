/// makes excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr T_(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes lambda deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr Lambda_(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes excitation operator of all bra/ket ranks up to (and including) @c Nbra/Nket
ExprPtr T(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes deexcitation operator of all bra/ket ranks up to (and including) @c Nbra/Nket
ExprPtr Lambda(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr A(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes L deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr L(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);

/// makes R excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr R(std::size_t Nbra, std::size_t Nket = std::numeric_limits<std::size_t>::max(), bool complete_unoccupieds = false);
