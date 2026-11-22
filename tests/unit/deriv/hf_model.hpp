#ifndef SEQUANT_TESTS_UNIT_DERIV_HF_MODEL_HPP
#define SEQUANT_TESTS_UNIT_DERIV_HF_MODEL_HPP

#include <Eigen/Dense>

#include <complex>
#include <cstddef>
#include <vector>

#include "evaluator.hpp"

/// A numerical model of the energy of mbpt::hf::energy() and
/// mbpt::hf::orthonormalized(), built from the evaluator's own element values:
/// the energy of one determinant, evaluated exactly at given perturbation
/// strengths x,
///
///   E(x) = Σ_i h̃_ii + 1/2 Σ_ij g̃_{ij,ij},  i over the occupied orbitals,
///
/// with h̃ = C† h(x) C and g̃ the same on each of the four legs of g(x). The
/// orbital coefficients are C = X · exp(K): orthonormalized, then rotated,
/// X = S(x)^(−1/2) (S Hermitian, by diagonalization) and exp(K) by
/// diagonalizing the Hermitian iK. K is the anti-Hermitian matrix of κ(x):
/// K_ai = κ_ai, K_ia = −conj(κ_ai), i occupied, a unoccupied, every other
/// element 0. Every quantity q is a series in x,
/// q(x) = q(0) + Σ_z x_z q¹_z + 1/2 Σ_zw x_z x_w q²_zw, with S(0) = 1 and
/// κ(0) = 0; for a scalar strength there is one z and no auxiliary index.
/// No Hartree-Fock equation is solved: the data is random.
namespace sequant::deriv::test::hf_model {

using Number = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
/// g{p,q;r,s} at ((p n + q) n + r) n + s, n orbitals
using Array4 = std::vector<Number>;

/// the orbitals: the values an index of the complete space runs over in the
/// evaluator, and the positions of the reference-occupied ones and of the
/// unoccupied ones among them
struct Orbitals {
  container::svector<int> values;
  container::svector<std::size_t> occupied, unoccupied;
};

/// a quantity as a series in the strengths: its value at zero, and the
/// coefficients of first order, one per z, and of second order, one per (z,w);
/// a missing coefficient is zero
template <typename T>
struct Series {
  T zeroth;
  std::vector<T> first;
  std::vector<std::vector<T>> second;
};

/// the data of the model; κ is held as an n × n matrix whose only nonzero
/// elements are κ_ai
struct Model {
  Orbitals orbitals;
  Series<Matrix> h, s, kappa;
  Series<Array4> g;
};

/// how the orbital coefficients are made from X and exp(K)
enum class Coefficients {
  /// C = X · exp(K), the model
  orthonormalized_then_rotated,
  /// C = exp(K) · X: agrees with the model at first order only
  rotated_then_orthonormalized
};

/// the model over the orbitals of the default context, read from the
/// evaluator with @p data ; the orbitals are the values an index of the
/// complete space runs over in a sum of the evaluator, spin cases distinct if
/// @p quantum_numbers (as element_of() makes them where the indices of its
/// input differ in quantum numbers). Read: h, g and their coefficients h¹, h²,
/// g¹, g² with the symmetries of hf::energy(), the overlap coefficients s¹, s²
/// spelled as hf::orthonormalized() writes them, and κ̂¹, κ̂² over (unoccupied;
/// occupied). @p z holds the values of the auxiliary index of a strength
/// tensor, which the coefficients carry after their slots (z_1 for the first
/// order, z_1,z_2 for the second); empty: a scalar strength, no auxiliary
/// index
Model read(const container::svector<int>& z, const Data& data,
           bool quantum_numbers);

/// E(x) of @p model at the strengths @p x (one per z), the orbital
/// coefficients made as @p coefficients says
Number energy(
    const Model& model, const std::vector<double>& x,
    Coefficients coefficients = Coefficients::orthonormalized_then_rotated);

}  // namespace sequant::deriv::test::hf_model

#endif  // SEQUANT_TESTS_UNIT_DERIV_HF_MODEL_HPP
