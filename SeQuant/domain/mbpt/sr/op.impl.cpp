/// makes excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr T_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return make_op(OpType::t, Nbra, Nket)();
}

/// makes lambda deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr Lambda_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return make_op(OpType::lambda, Nbra, Nket)();
}

/// makes R excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr R_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 || Nket > 0);
  return make_op(OpType::R, Nbra, Nket)();
}

/// makes L deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr L_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 || Nket > 0);
  return make_op(OpType::L, Nbra, Nket)();
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
    else if (op_ == OpType::lambda)
      result = result ? result + Lambda_(nbra_, nket_) : Lambda_(nbra_, nket_);
    else if (op_ == OpType::R)
      result = result ? result + R_(nbra_, nket_) : R_(nbra_, nket_);
    else if (op_ == OpType::L)
      result = result ? result + L_(nbra_, nket_) : L_(nbra_, nket_);
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

/// makes deexcitation operator of all bra/ket ranks up to (and including)
/// @c Nbra/Nket
ExprPtr Lambda(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ =
      Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nket_ > 0);
  ExprPtr result;
  detail::op_impl{OpType::lambda, Nbra, Nket_}(result);
  return result;
}

/// makes geminal excitation operator for ansatz @p ansatz
ExprPtr R12(IndexSpace::Type gg_space, int ansatz) {
  assert(ansatz == 1 || ansatz == 2);
  if (ansatz == 2)
    return make_op(
        OpType::R12,
        {IndexSpace::complete_unoccupied, IndexSpace::complete_unoccupied},
        {gg_space, gg_space})();
  else
    return make_op(OpType::R12,
                   {IndexSpace::other_unoccupied, IndexSpace::other_unoccupied},
                   {gg_space, gg_space})();
}

/// makes deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr A(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return make_op(OpType::A, Nbra, Nket)();
}

/// makes excitation operator of all bra/ket ranks up to (and including)
/// @c Nbra/Nket
ExprPtr R(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ =
      Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nket_ > 0);
  ExprPtr result;
  detail::op_impl{OpType::R, Nbra, Nket_}(result);
  return result;
}

/// makes deexcitation operator of all bra/ket ranks up to (and including)
/// @c Nbra/Nket
ExprPtr L(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ =
      Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nket_ > 0);
  ExprPtr result;
  detail::op_impl{OpType::L, Nbra, Nket_}(result);
  return result;
}
