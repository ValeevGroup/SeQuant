/// makes excitation operator of bra/ket ranks @c Nbra/Nket
template <std::size_t Nbra, std::size_t Nket = Nbra>
ExprPtr T_() {
  return Op<OpType::t, Nbra, Nket>();
}

/// makes lambda deexcitation operator of bra/ket ranks @c Nbra/Nket
template <std::size_t Nbra, std::size_t Nket = Nbra>
ExprPtr Lambda_() {
  return Op<OpType::l, Nbra, Nket>();
}

namespace detail {

template <OpType Op, size_t Nbra, size_t Nket = Nbra>
void op_impl(ExprPtr& result) {
  static_assert(Nbra > 0 || Nket > 0);
  if constexpr (Op == OpType::t)
    result = result ? result + T_<Nbra, Nket>() : T_<Nbra, Nket>();
  else if constexpr (Op == OpType::l)
    result = result ? result + Lambda_<Nbra, Nket>() : Lambda_<Nbra, Nket>();
  else
    assert(false && "unsupported Op value");

  if constexpr ((Nbra > 1 && Nket > 0) || (Nbra > 0 && Nket > 1)) {
    op_impl<Op, Nbra - 1, Nket - 1>(result);
  }
}
}  // namespace detail

/// makes excitation operator of all bra/ket ranks up to (and including) @c Nbra/Nket
template <std::size_t Nbra, std::size_t Nket = Nbra>
ExprPtr T() {
  static_assert(Nbra > 0 || Nket > 0);
  ExprPtr result;
  detail::op_impl<OpType::t, Nbra, Nket>(result);
  return result;
}

/// makes deexcitation operator of all bra/ket ranks up to (and including) @c Nbra/Nket
template <std::size_t Nbra, std::size_t Nket = Nbra>
ExprPtr Lambda() {
  static_assert(Nbra > 0 || Nket > 0);
  ExprPtr result;
  detail::op_impl<OpType::l, Nbra, Nket>(result);
  return result;
}

/// makes deexcitation operator of bra/ket ranks @c Nbra/Nket
template <std::size_t Nbra, std::size_t Nket = Nbra>
ExprPtr A(bool complete_unoccupieds = false) {
  return Op<OpType::A, Nbra, Nket>(complete_unoccupieds);
}

/// makes L deexcitation operator of bra/ket ranks @c Nbra/Nket
template <std::size_t Nbra, std::size_t Nket = Nbra>
ExprPtr L(bool complete_unoccupieds = false) {
  return Op<OpType::L, Nbra, Nket>(complete_unoccupieds);
}

/// makes R excitation operator of bra/ket ranks @c Nbra/Nket
template <std::size_t Nbra, std::size_t Nket = Nbra>
ExprPtr R(bool complete_unoccupieds = false) {
  return Op<OpType::R, Nbra, Nket>(complete_unoccupieds);
}
