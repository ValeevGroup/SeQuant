/// \brief Provides one row of the NNS projector matrix,
/// hardcoded from Mathematica to avoid numerical precision loss.
///
/// The Myrvold-Ruskey unrank1 algorithm (doi.org/10.1016/S0020-0190(01)00141-7)
/// was used to order permutations, then computing the permutational overlap
/// matrix (M) where elements are (-2)^{cycles} × (-1)^n.
/// The NNS projector weights are obtained from the normalized pseudoinverse
/// of M: first compute M_pinv (the pseudoinverse), then normalize it by
/// multiplying by (rank / rows) where rank is the matrix rank and rows is n!.
/// Finally, NNS projector = normalized_M_pinv · M
///
/// \param n_particles The rank of external index pairs
///
/// \return Vector of NNS projector weights representing the last row
///
/// \throw std::runtime_error if n_particles is not in the range [1,5]
template <typename T>
  requires(std::floating_point<T> || meta::is_complex_v<T>)
std::vector<T> hardcoded_nns_projector(std::size_t n_particles) {
  switch (n_particles) {
    case 1:
      return std::vector<T>{T(1) / T(1)};

    case 2:
      return std::vector<T>{T(1) / T(1), T(1) / T(1)};

    case 3:
      return std::vector<T>{T(-1) / T(5), T(-1) / T(5), T(-1) / T(5),
                            T(-1) / T(5), T(-1) / T(5), T(1) / T(1)};

    case 4:
      return std::vector<T>{
          T(1) / T(7),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(14),
          T(1) / T(7),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(1) / T(7),   T(-2) / T(7),
          T(-1) / T(14), T(1) / T(7),   T(-1) / T(14), T(-2) / T(7),
          T(1) / T(7),   T(-1) / T(14), T(-1) / T(14), T(-2) / T(7),
          T(-2) / T(7),  T(-2) / T(7),  T(-2) / T(7),  T(1) / T(1)};

    case 5:
      return std::vector<T>{
          T(-1) / T(14), T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(2) / T(21),  T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(2) / T(21),  T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(0) / T(1),
          T(-1) / T(14), T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(2) / T(21),  T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(2) / T(21),  T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(0) / T(1),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(2) / T(21),
          T(0) / T(1),   T(2) / T(21),  T(2) / T(21),  T(-1) / T(21),
          T(2) / T(21),  T(0) / T(1),   T(-1) / T(21), T(-1) / T(21),
          T(-1) / T(21), T(-1) / T(21), T(1) / T(7),   T(0) / T(1),
          T(0) / T(1),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(3),
          T(2) / T(21),  T(-1) / T(21), T(2) / T(21),  T(2) / T(21),
          T(0) / T(1),   T(-1) / T(21), T(-1) / T(21), T(-1) / T(21),
          T(-1) / T(21), T(1) / T(7),   T(2) / T(21),  T(-1) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(0) / T(1),   T(0) / T(1),
          T(1) / T(7),   T(0) / T(1),   T(1) / T(7),   T(-1) / T(3),
          T(-1) / T(21), T(-1) / T(21), T(-1) / T(21), T(-1) / T(21),
          T(1) / T(7),   T(-1) / T(21), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(0) / T(1),   T(-1) / T(21), T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(0) / T(1),   T(1) / T(7),
          T(0) / T(1),   T(0) / T(1),   T(1) / T(7),   T(-1) / T(3),
          T(0) / T(1),   T(1) / T(7),   T(1) / T(7),   T(0) / T(1),
          T(-1) / T(3),  T(1) / T(7),   T(0) / T(1),   T(1) / T(7),
          T(0) / T(1),   T(-1) / T(3),  T(1) / T(7),   T(1) / T(7),
          T(0) / T(1),   T(0) / T(1),   T(-1) / T(3),  T(-1) / T(3),
          T(-1) / T(3),  T(-1) / T(3),  T(-1) / T(3),  T(1) / T(1)};

    default:
      throw std::runtime_error(
          "hardcoded NNS weights only available for ranks 1-5, "
          "requested rank is : " +
          std::to_string(n_particles));
  }
}
