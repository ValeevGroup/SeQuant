/// makes excitation operator of rank @c N
template <std::size_t N>
ExprPtr T_() {
  return Op<N, OpType::t>();
}

/// makes lambda deexcitation operator of rank @c N
template <std::size_t N>
ExprPtr L_() {
  return Op<N, OpType::l>();
}

namespace detail {

template <size_t N, OpType Op>
void op_impl(ExprPtr& result) {
  static_assert(N > 0);
  if constexpr (Op == OpType::t)
    result = result ? result + T_<N>() : T_<N>();
  else if constexpr (Op == OpType::l)
    result = result ? result + L_<N>() : L_<N>();
  else
    assert(false && "unsupported Op value");

  if constexpr (N > 1) {
    op_impl<N - 1, Op>(result);
  }
}
}  // namespace detail

/// makes excitation operator of all ranks up to (and including) @c N
template <std::size_t N>
ExprPtr T() {
  static_assert(N > 0);
  ExprPtr result;
  detail::op_impl<N, OpType::t>(result);
  return result;
}

/// makes excitation operator of all ranks up to (and including) @c N
template <std::size_t N>
ExprPtr L() {
  static_assert(N > 0);
  ExprPtr result;
  detail::op_impl<N, OpType::l>(result);
  return result;
}

/// makes A deexcitation operator of rank @c N
template <std::size_t N>
ExprPtr A() {
  return Op<N, OpType::A>();
}
