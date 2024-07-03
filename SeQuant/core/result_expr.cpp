#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/indices.hpp>

namespace sequant {

ResultExpr::ResultExpr(const Tensor &tensor, ExprPtr expression)
    : m_expr(std::move(expression)),
      m_symm(tensor.symmetry()),
      m_bksymm(tensor.braket_symmetry()),
      m_psymm(tensor.particle_symmetry()),
      m_braIndices(tensor.bra().begin(), tensor.bra().end()),
      m_ketIndices(tensor.ket().begin(), tensor.ket().end()),
      m_label(tensor.label()) {}

ResultExpr::ResultExpr(const Variable &variable, ExprPtr expression)
    : m_expr(std::move(expression)), m_label(variable.label()) {}

ResultExpr &ResultExpr::operator=(ExprPtr expression) {
  m_expr = std::move(expression);

  return *this;
}

bool ResultExpr::has_label() const { return m_label.has_value(); }

const std::wstring &ResultExpr::label() const { return m_label.value(); }

Symmetry ResultExpr::symmetry() const { return m_symm; }

void ResultExpr::set_symmetry(Symmetry symm) { m_symm = symm; }

BraKetSymmetry ResultExpr::braket_symmetry() const { return m_bksymm; }

void ResultExpr::set_braket_symmetry(BraKetSymmetry symm) { m_bksymm = symm; }

ParticleSymmetry ResultExpr::particle_symmetry() const { return m_psymm; }

void ResultExpr::set_particle_symmetry(ParticleSymmetry symm) {
  m_psymm = symm;
}

const ResultExpr::IndexContainer &ResultExpr::bra() const {
  return m_braIndices;
}

const ResultExpr::IndexContainer &ResultExpr::ket() const {
  return m_ketIndices;
}

const ExprPtr &ResultExpr::expression() const { return m_expr; }

ExprPtr &ResultExpr::expression() { return m_expr; }

}  // namespace sequant
