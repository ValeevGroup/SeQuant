ExprPtr T_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return OpMaker(OpType::t, Nbra, Nket)();
}

ExprPtr Lambda_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return OpMaker(OpType::λ, Nbra, Nket)();
}

namespace detail {

/// constructs a sum of ops up to a given bra/ket rank
class op_impl {
  OpType op_;
  std::size_t nbra_, nket_;

 public:
  op_impl(OpType op, std::size_t nbra,
          std::size_t nket = std::numeric_limits<std::size_t>::max())
      : op_(op),
        nbra_(nbra),
        nket_(nket == std::numeric_limits<std::size_t>::max() ? nbra : nket) {
    assert(nbra_ > 0 && nbra_ < std::numeric_limits<std::size_t>::max());
    assert(nket_ > 0);
  }

  void operator()(ExprPtr& result) {
    if (op_ == OpType::t)
      result = result ? result + T_(nbra_, nket_) : T_(nbra_, nket_);
    else
      assert(false && "unsupported op value");

    if ((nbra_ > 1 && nket_ > 0) || (nbra_ > 0 && nket_ > 1)) {
      op_impl{op_, nbra_ - 1, nket_ - 1}(result);
    }
  }
};

}  // namespace detail

/// makes excitation operator of all bra/ket ranks up to (and including)
/// @c Nbra/Nket
ExprPtr T(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ =
      Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nket_ > 0);
  ExprPtr result;
  detail::op_impl{OpType::t, Nbra, Nket_}(result);
  return result;
}

/// makes geminal excitation operator for ansatz @p ansatz
ExprPtr R12(IndexSpace::Type gg_space, int ansatz) {
  assert(ansatz == 1 || ansatz == 2);
  if (ansatz == 2)
    return OpMaker(
        OpType::R12,
        {IndexSpace::complete_unoccupied, IndexSpace::complete_unoccupied},
        {gg_space, gg_space})();
  else
    return OpMaker(OpType::R12,
                   {IndexSpace::other_unoccupied, IndexSpace::other_unoccupied},
                   {gg_space, gg_space})();
}

/// makes deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr A(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return OpMaker(OpType::A, Nbra, Nket)();
}
