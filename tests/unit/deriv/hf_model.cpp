#include "hf_model.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <algorithm>
#include <functional>
#include <string>

namespace sequant::deriv::test::hf_model {

namespace {

using Positions = container::svector<std::size_t>;

/// the values an index of @p space runs over in a sum of the evaluator
container::svector<int> values_of(const IndexSpace& space, bool quantum_numbers,
                                  const Data& data) {
  Domains domains;
  add_domain(Index(std::wstring(space.base_key()) + L"_1"), {}, domains,
             quantum_numbers, data.occupied);
  return domains.front().second;
}

Orbitals orbitals_of(bool quantum_numbers, const Data& data) {
  const auto& registry = *get_default_context().index_basis_registry();
  const auto values = [&](const IndexSpace& s) {
    return values_of(s, quantum_numbers, data);
  };
  Orbitals result;
  result.values = values(registry.complete_space(mbpt::Spin::any));
  const auto occupied =
      values(registry.reference_occupied_space(mbpt::Spin::any));
  const auto unoccupied =
      values(registry.vacuum_unoccupied_space(mbpt::Spin::any));
  for (std::size_t k = 0; k != result.values.size(); ++k) {
    const auto v = result.values[k];
    if (std::ranges::find(occupied, v) != occupied.end())
      result.occupied.push_back(k);
    else if (std::ranges::find(unoccupied, v) != unoccupied.end())
      result.unoccupied.push_back(k);
    else
      throw Exception("hf_model: an orbital neither occupied nor unoccupied");
  }
  return result;
}

/// the matrix of the one-leaf expression @p text over the labels @p row and
/// @p column at the positions @p rows and @p columns , zero elsewhere
Matrix matrix(const Orbitals& o, const std::wstring& text,
              const std::wstring& row, const Positions& rows,
              const std::wstring& column, const Positions& columns,
              const Values& fixed, const Data& data) {
  const auto leaf = deserialize(text);
  const auto n = o.values.size();
  Matrix result = Matrix::Zero(n, n);
  auto values = fixed;
  for (const auto r : rows)
    for (const auto c : columns) {
      values[row] = o.values[r];
      values[column] = o.values[c];
      result(r, c) = element_of(leaf, values, data);
    }
  return result;
}

/// the array of the one-leaf expression @p text over p_1, p_2, p_3, p_4
Array4 array4(const Orbitals& o, const std::wstring& text, const Values& fixed,
              const Data& data) {
  const auto leaf = deserialize(text);
  const auto n = o.values.size();
  Array4 result(n * n * n * n);
  auto values = fixed;
  std::size_t at = 0;
  for (std::size_t a = 0; a != n; ++a)
    for (std::size_t b = 0; b != n; ++b)
      for (std::size_t c = 0; c != n; ++c)
        for (std::size_t d = 0; d != n; ++d) {
          values[L"p_1"] = o.values[a];
          values[L"p_2"] = o.values[b];
          values[L"p_3"] = o.values[c];
          values[L"p_4"] = o.values[d];
          result[at++] = element_of(leaf, values, data);
        }
  return result;
}

/// the coefficients of first and second order of a series, read by @p read
/// from the auxiliary slots of a coefficient and their values, the labels
/// @p first and @p second given
template <typename T>
void read_coefficients(
    Series<T>& series, const container::svector<int>& z,
    const std::wstring& first, const std::wstring& second,
    const std::function<T(const std::wstring&, const std::wstring&,
                          const Values&)>& read) {
  if (z.empty()) {
    series.first = {read(first, L"", {})};
    series.second = {{read(second, L"", {})}};
    return;
  }
  for (const auto zj : z) {
    series.first.push_back(read(first, L";z_1", {{L"z_1", zj}}));
    series.second.emplace_back();
    for (const auto zk : z)
      series.second.back().push_back(
          read(second, L";z_1,z_2", {{L"z_1", zj}, {L"z_2", zk}}));
  }
}

void add(Matrix& to, double c, const Matrix& q) { to += c * q; }

void add(Array4& to, double c, const Array4& q) {
  for (std::size_t k = 0; k != to.size(); ++k) to[k] += c * q[k];
}

/// q(x) of the series @p q
template <typename T>
T at(const Series<T>& q, const std::vector<double>& x) {
  T result = q.zeroth;
  for (std::size_t z = 0; z != q.first.size(); ++z)
    add(result, x[z], q.first[z]);
  for (std::size_t z = 0; z != q.second.size(); ++z)
    for (std::size_t w = 0; w != q.second[z].size(); ++w)
      add(result, x[z] * x[w] / 2, q.second[z][w]);
  return result;
}

/// S^(−1/2) of the Hermitian positive definite @p S
Matrix inverse_sqrt(const Matrix& S) {
  const Eigen::SelfAdjointEigenSolver<Matrix> solver(S);
  if (solver.eigenvalues().minCoeff() <= 0)
    throw Exception("hf_model: the overlap is not positive definite");
  const Eigen::VectorXd d = solver.eigenvalues().array().rsqrt();
  return solver.eigenvectors() * d.asDiagonal() *
         solver.eigenvectors().adjoint();
}

/// exp(K) of the anti-Hermitian @p K : K = −iH with H = iK Hermitian
Matrix exp_antihermitian(const Matrix& K) {
  const Eigen::SelfAdjointEigenSolver<Matrix> solver(Number(0, 1) * K);
  const Eigen::VectorXcd d =
      (Number(0, -1) * solver.eigenvalues().cast<Number>()).array().exp();
  return solver.eigenvectors() * d.asDiagonal() *
         solver.eigenvectors().adjoint();
}

/// throws "hf_model: @p what" unless every coefficient of @p q satisfies
/// @p holds
template <typename T, typename P>
void require(const Series<T>& q, const P& holds, const std::string& what) {
  bool ok = holds(q.zeroth);
  for (const auto& c : q.first) ok = ok && holds(c);
  for (const auto& row : q.second)
    for (const auto& c : row) ok = ok && holds(c);
  if (!ok) throw Exception("hf_model: " + what);
}

}  // namespace

Model read(const container::svector<int>& z, const Data& data,
           bool quantum_numbers) {
  const auto o = orbitals_of(quantum_numbers, data);
  const auto n = o.values.size();
  Positions all(n);
  for (std::size_t k = 0; k != n; ++k) all[k] = k;
  const auto& registry = *get_default_context().index_basis_registry();
  const auto label = [&](const IndexSpace& s) {
    return std::wstring(s.base_key()) + L"_1";
  };
  const auto occupied =
      label(registry.reference_occupied_space(mbpt::Spin::any));
  const auto unoccupied =
      label(registry.vacuum_unoccupied_space(mbpt::Spin::any));

  // over p_1 and p_2 with the symmetries `flags`
  const auto one_body = [&](const std::wstring& flags) {
    return [&o, &all, &data, flags](const std::wstring& l,
                                    const std::wstring& aux, const Values& v) {
      return matrix(o, l + L"{p_1;p_2" + aux + L"}:" + flags, L"p_1", all,
                    L"p_2", all, v, data);
    };
  };
  const auto two_body = [&](const std::wstring& l, const std::wstring& aux,
                            const Values& v) {
    return array4(o, l + L"{p_1,p_2;p_3,p_4" + aux + L"}:A-C-S", v, data);
  };
  const auto kappa = [&](const std::wstring& l, const std::wstring& aux,
                         const Values& v) {
    return matrix(o, l + L"{" + unoccupied + L";" + occupied + aux + L"}:A-N-S",
                  unoccupied, o.unoccupied, occupied, o.occupied, v, data);
  };

  Model result;
  result.orbitals = o;
  result.h.zeroth = one_body(L"A-C-S")(L"h", L"", {});
  read_coefficients<Matrix>(result.h, z, L"h¹", L"h²", one_body(L"A-C-S"));
  result.s.zeroth = Matrix::Identity(n, n);
  read_coefficients<Matrix>(result.s, z, L"s¹", L"s²", one_body(L"N-C-S"));
  result.kappa.zeroth = Matrix::Zero(n, n);
  read_coefficients<Matrix>(result.kappa, z, L"κₒ¹", L"κₒ²", kappa);
  result.g.zeroth = two_body(L"g", L"", {});
  read_coefficients<Array4>(result.g, z, L"g¹", L"g²", two_body);

  const auto hermitian = [](const Matrix& q) {
    return q.isApprox(q.adjoint());
  };
  require(result.h, hermitian, "h is not Hermitian");
  require(result.s, hermitian, "S is not Hermitian");
  require(
      result.g,
      [n](const Array4& g) {
        const auto at = [n](std::size_t p, std::size_t q, std::size_t r,
                            std::size_t s) {
          return ((p * n + q) * n + r) * n + s;
        };
        const auto close = [](Number x, Number y) {
          return std::abs(x - y) <= 1e-12 * std::max(1., std::abs(x));
        };
        for (std::size_t p = 0; p != n; ++p)
          for (std::size_t q = 0; q != n; ++q)
            for (std::size_t r = 0; r != n; ++r)
              for (std::size_t s = 0; s != n; ++s) {
                const auto x = g[at(p, q, r, s)];
                if (!close(x, -g[at(q, p, r, s)]) ||
                    !close(x, -g[at(p, q, s, r)]) ||
                    !close(x, std::conj(g[at(r, s, p, q)])))
                  return false;
              }
        return true;
      },
      "g is not antisymmetric and Hermitian");
  return result;
}

Number energy(const Model& model, const std::vector<double>& x,
              Coefficients coefficients) {
  const auto& o = model.orbitals;
  const auto n = o.values.size();
  const Matrix X = inverse_sqrt(at(model.s, x));
  const Matrix kappa = at(model.kappa, x);
  const Matrix U = exp_antihermitian(kappa - kappa.adjoint());
  const Matrix C = coefficients == Coefficients::orthonormalized_then_rotated
                       ? X * U
                       : U * X;
  // ρ_ba = Σ_i C_bi conj(C_ai): Σ_i h̃_ii = Σ_ab h_ab ρ_ba and
  // Σ_ij g̃_{ij,ij} = Σ_abcd g_abcd ρ_ca ρ_db
  Matrix occupied(n, o.occupied.size());
  for (std::size_t i = 0; i != o.occupied.size(); ++i)
    occupied.col(i) = C.col(o.occupied[i]);
  const Matrix rho = occupied * occupied.adjoint();
  const Matrix h = at(model.h, x);
  const auto g = at(model.g, x);
  Number two_body = 0;
  std::size_t k = 0;
  for (std::size_t a = 0; a != n; ++a)
    for (std::size_t b = 0; b != n; ++b)
      for (std::size_t c = 0; c != n; ++c)
        for (std::size_t d = 0; d != n; ++d)
          two_body += g[k++] * rho(c, a) * rho(d, b);
  return (h * rho).trace() + two_body / 2.;
}

}  // namespace sequant::deriv::test::hf_model
