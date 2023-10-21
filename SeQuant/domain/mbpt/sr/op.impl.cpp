/// makes excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr T_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return OpMaker(OpType::t, Nbra, Nket)();
}

/// makes lambda deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr Λ_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return OpMaker(OpType::λ, Nbra, Nket)();
}

/// makes R excitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr R_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 || Nket > 0);
  return OpMaker(OpType::R, Nbra, Nket)();
}

/// makes L deexcitation operator of bra/ket ranks @c Nbra/Nket
ExprPtr L_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 || Nket > 0);
  return OpMaker(OpType::L, Nbra, Nket)();
}

ExprPtr mu(std::size_t R) { return OpMaker(OpType::μ, R)(); }

ExprPtr T_pt_(std::size_t o, std::size_t Nbra, std::size_t Nket) {
  assert(o == 1 &&
         "sequant::sr::T_pt_(): only supports first order perturbation");
  assert(Nbra > 0);
  assert(Nket > 0);
  return OpMaker(OpType::t_1, Nbra, Nket)();
}

ExprPtr Λ_pt_(std::size_t o, std::size_t Nbra, std::size_t Nket) {
  assert(o == 1 &&
         "sequant::sr::Λ_pt_(): only supports first order perturbation");
  assert(Nbra > 0);
  assert(Nket > 0);
  return OpMaker(OpType::λ_1, Nbra, Nket)();
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
    else if (op_ == OpType::λ)
      result = result ? result + Λ_(nbra_, nket_) : Λ_(nbra_, nket_);
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
ExprPtr Λ(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0 && Nbra < std::numeric_limits<std::size_t>::max());
  const auto Nket_ =
      Nket == std::numeric_limits<std::size_t>::max() ? Nbra : Nket;
  assert(Nket_ > 0);
  ExprPtr result;
  detail::op_impl{OpType::λ, Nbra, Nket_}(result);
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

ExprPtr A(std::int64_t Kh, std::int64_t Kp) {
  assert(Kh != 0);
  if (Kp == std::numeric_limits<std::int64_t>::max()) Kp = Kh;
  assert(Kp != 0);

  // Kh and Kp should have same sign
  assert((Kh > 0 && Kp > 0) || (Kh < 0 && Kp < 0));

  container::svector<IndexSpace::Type> creators;
  container::svector<IndexSpace::Type> annihilators;
  if (Kh > 0)
    for (auto i : ranges::views::iota(0, Kh))
      annihilators.emplace_back(IndexSpace::active_occupied);
  else
    for (auto i : ranges::views::iota(0, -Kh))
      creators.emplace_back(IndexSpace::active_occupied);
  if (Kp > 0)
    for (auto i : ranges::views::iota(0, Kp))
      creators.emplace_back(IndexSpace::active_unoccupied);
  else
    for (auto i : ranges::views::iota(0, -Kp))
      annihilators.emplace_back(IndexSpace::active_unoccupied);

  std::optional<OpMaker::UseDepIdx> dep;
  if (get_default_formalism().csv() == mbpt::CSV::Yes)
    dep = Kh > 0 ? OpMaker::UseDepIdx::Bra : OpMaker::UseDepIdx::Ket;
  return OpMaker(OpType::A, creators, annihilators)(dep, {Symmetry::antisymm});
}

ExprPtr S(std::int64_t Kh, std::int64_t Kp) {
  assert(Kh != 0);
  if (Kp == std::numeric_limits<std::int64_t>::max()) Kp = Kh;
  assert(Kp != 0);

  // Kh and Kp should have same sign
  assert((Kh > 0 && Kp > 0) || (Kh < 0 && Kp < 0));

  container::svector<IndexSpace::Type> creators;
  container::svector<IndexSpace::Type> annihilators;
  if (Kh > 0)
    for (auto i : ranges::views::iota(0, Kh))
      annihilators.emplace_back(IndexSpace::active_occupied);
  else
    for (auto i : ranges::views::iota(0, -Kh))
      creators.emplace_back(IndexSpace::active_occupied);
  if (Kp > 0)
    for (auto i : ranges::views::iota(0, Kp))
      creators.emplace_back(IndexSpace::active_unoccupied);
  else
    for (auto i : ranges::views::iota(0, -Kp))
      annihilators.emplace_back(IndexSpace::active_unoccupied);

  std::optional<OpMaker::UseDepIdx> dep;
  if (get_default_formalism().csv() == mbpt::CSV::Yes)
    dep = Kh > 0 ? OpMaker::UseDepIdx::Bra : OpMaker::UseDepIdx::Ket;
  return OpMaker(OpType::S, creators, annihilators)(dep, {Symmetry::nonsymm});
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
