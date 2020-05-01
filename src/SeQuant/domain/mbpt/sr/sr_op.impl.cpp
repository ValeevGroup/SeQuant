/// makes excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr T_(std::size_t Nbra, std::size_t Nket, bool complete_unoccupieds) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  assert(Nket > 0);
  const auto Nket_ = Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  return Op(OpType::t, Nbra, Nket_)(complete_unoccupieds);
}

/// makes lambda deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr Lambda_(std::size_t Nbra, std::size_t Nket, bool complete_unoccupieds) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  assert(Nket > 0);
  const auto Nket_ = Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  return Op(OpType::l, Nbra, Nket_)(complete_unoccupieds);
}

namespace detail {

/// constructs a sum of ops up to a given bra/ket rank
class op_impl {

  OpType op_;
  std::size_t nbra_, nket_;

 public:
  op_impl(OpType op, std::size_t nbra, std::size_t nket = std::numeric_limits<std::size_t>::max()) :
  op_(op), nbra_(nbra), nket_(nket == std::numeric_limits<std::size_t>::max() ? nbra : nket) {
    assert(nbra_ > 0 && nbra_ < std::numeric_limits<std::size_t>::max());
    assert(nket_ > 0);
  }

  void operator()(ExprPtr& result, bool complete_unoccupieds = false) {
    if (op_ == OpType::t)
      result = result ? result + T_(nbra_, nket_, complete_unoccupieds) : T_(nbra_, nket_, complete_unoccupieds);
    else if (op_ == OpType::l)
      result = result ? result + Lambda_(nbra_, nket_, complete_unoccupieds) : Lambda_(nbra_, nket_, complete_unoccupieds);
    else
      assert(false && "unsupported op value");

    if ((nbra_ > 1 && nket_ > 0) || (nbra_ > 0 && nket_ > 1)) {
      op_impl{op_, nbra_ - 1, nket_ - 1}(result, complete_unoccupieds);
    }
  }

};

}  // namespace detail

/// makes excitation operator of all bra/ket ranks up to (and including) @c Nbra/Nket
ExprPtr T(std::size_t Nbra, std::size_t Nket, bool complete_unoccupieds) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ = Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nket_ > 0);
  ExprPtr result;
  detail::op_impl{OpType::t, Nbra, Nket_}(result, complete_unoccupieds);
  return result;
}

/// makes deexcitation operator of all bra/ket ranks up to (and including) @c Nbra/Nket
ExprPtr Lambda(std::size_t Nbra, std::size_t Nket, bool complete_unoccupieds) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ = Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nket_ > 0);
  ExprPtr result;
  detail::op_impl{OpType::l, Nbra, Nket_}(result, complete_unoccupieds);
  return result;
}

/// makes geminal excitation operator for ansatz @p ansatz
ExprPtr R12(int ansatz) {
  assert(ansatz == 1 || ansatz == 2);
  return Op(OpType::R12, 2, 2)(ansatz == 2 ? IndexSpace::complete_unoccupied : IndexSpace::other_unoccupied);
}

/// makes deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr A(std::size_t Nbra, std::size_t Nket, bool complete_unoccupieds) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ = Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nket_ > 0);
  return Op(OpType::A, Nbra, Nket_)(complete_unoccupieds);
}

/// makes L deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr L(std::size_t Nbra, std::size_t Nket, bool complete_unoccupieds) {
  assert(Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ = Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nbra > 0 || Nket_ > 0);
  return Op(OpType::L, Nbra, Nket_)(complete_unoccupieds);
}

/// makes R excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr R(std::size_t Nbra, std::size_t Nket, bool complete_unoccupieds) {
  assert(Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ = Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nbra > 0 || Nket_ > 0);
  return Op(OpType::R, Nbra, Nket_)(complete_unoccupieds);
}
