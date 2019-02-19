/// makes excitation operator of rank @c N
template <std::size_t N>
ExprPtr T_() {
  return Op<N, OpType::t>();
}

/// makes excitation operator of all ranks up to (and including) @c N
template <std::size_t N>
ExprPtr T() {
  assert(N > 0);
  ExprPtr result = T_<1>();
  if (N > 1) result = result + T_<2>();
  if (N > 2) result = result + T_<3>();
  if (N > 3) result = result + T_<4>();
  if (N > 4)
    throw std::invalid_argument(
        "sequant2::mbpt::sr::so::T<N>: N>4 is not supported");
  return result;
}

/// makes lambda deexcitation operator of rank @c N
template <std::size_t N>
ExprPtr L_() {
  return Op<N, OpType::l>();
}

/// makes excitation operator of all ranks up to (and including) @c N
template <std::size_t N>
ExprPtr L() {
  assert(N > 0);
  ExprPtr result = L_<1>();
  if (N > 1) result = result + L_<2>();
  if (N > 2) result = result + L_<3>();
  if (N > 3) result = result + L_<4>();
  if (N > 4)
    throw std::invalid_argument(
        "sequant2::mbpt::sr::so::L<N>: N>4 is not supported");
  return result;
}

/// makes A deexcitation operator of rank @c N
template <std::size_t N>
ExprPtr A() {
  return Op<N, OpType::A>();
}
