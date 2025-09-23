#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/utility/indices.hpp>

namespace sequant {

bool ResultExpr::ResultCmp::operator()(const ResultExpr &lhs,
                                       const ResultExpr &rhs) const {
  if (lhs.has_label() != rhs.has_label()) {
    return false;
  }

  if (lhs.has_label() && lhs.label() != rhs.label()) {
    return false;
  }

  return lhs.symmetry() == rhs.symmetry() &&
         lhs.braket_symmetry() == rhs.braket_symmetry() &&
         lhs.column_symmetry() == rhs.column_symmetry() &&
         lhs.bra() == rhs.bra() && lhs.ket() == rhs.ket() &&
         lhs.aux() == rhs.aux();
}

ResultExpr::ResultExpr(const Tensor &tensor, ExprPtr expression)
    : m_expr(std::move(expression)),
      m_symm(tensor.symmetry()),
      m_bksymm(tensor.braket_symmetry()),
      m_csymm(tensor.column_symmetry()),
      m_braIndices(tensor.bra().begin(), tensor.bra().end()),
      m_ketIndices(tensor.ket().begin(), tensor.ket().end()),
      m_auxIndices(tensor.aux().begin(), tensor.aux().end()),
      m_label(tensor.label()) {}

ResultExpr::ResultExpr(const Variable &variable, ExprPtr expression)
    : m_expr(std::move(expression)), m_label(variable.label()) {}

ResultExpr::ResultExpr(IndexContainer bra, IndexContainer ket,
                       IndexContainer aux, Symmetry symm,
                       BraKetSymmetry braket_symm, ColumnSymmetry column_symm,
                       std::optional<std::wstring> label, ExprPtr expression)
    : m_expr(std::move(expression)),
      m_symm(symm),
      m_bksymm(braket_symm),
      m_csymm(column_symm),
      m_braIndices(std::move(bra)),
      m_ketIndices(std::move(ket)),
      m_auxIndices(std::move(aux)),
      m_label(std::move(label)) {}

ResultExpr &ResultExpr::operator=(ExprPtr expression) {
  m_expr = std::move(expression);

  return *this;
}

bool ResultExpr::has_label() const { return m_label.has_value(); }

const std::wstring &ResultExpr::label() const { return m_label.value(); }

void ResultExpr::set_label(std::wstring label) { m_label = std::move(label); }

Symmetry ResultExpr::symmetry() const { return m_symm; }

void ResultExpr::set_symmetry(Symmetry symm) { m_symm = symm; }

BraKetSymmetry ResultExpr::braket_symmetry() const { return m_bksymm; }

void ResultExpr::set_braket_symmetry(BraKetSymmetry symm) { m_bksymm = symm; }

ColumnSymmetry ResultExpr::column_symmetry() const { return m_csymm; }

void ResultExpr::set_column_symmetry(ColumnSymmetry symm) { m_csymm = symm; }

const ResultExpr::IndexContainer &ResultExpr::bra() const {
  return m_braIndices;
}

const ResultExpr::IndexContainer &ResultExpr::ket() const {
  return m_ketIndices;
}

const ResultExpr::IndexContainer &ResultExpr::aux() const {
  return m_auxIndices;
}

const ExprPtr &ResultExpr::expression() const { return m_expr; }

ExprPtr &ResultExpr::expression() { return m_expr; }

ResultExpr ResultExpr::clone() const {
  return ResultExpr(m_braIndices, m_ketIndices, m_auxIndices, m_symm, m_bksymm,
                    m_csymm, m_label, m_expr->clone());
}

ResultExpr &canonicalize(ResultExpr &expr) {
  expr.expression() = canonicalize(expr.expression());

  return expr;
}

ResultExpr &simplify(ResultExpr &expr) {
  expr.expression() = simplify(expr.expression());

  return expr;
}

ResultExpr &rapid_simplify(ResultExpr &expr) {
  expr.expression() = rapid_simplify(expr.expression());

  return expr;
}

ResultExpr &expand(ResultExpr &expr) {
  expr.expression() = expand(expr.expression());

  return expr;
}

ResultExpr &optimize(ResultExpr &expr) {
  expr.expression() = optimize(expr.expression());

  return expr;
}

ResultExpr &canonicalize(ResultExpr &&expr) { return canonicalize(expr); }

ResultExpr &simplify(ResultExpr &&expr) { return simplify(expr); }

ResultExpr &rapid_simplify(ResultExpr &&expr) { return rapid_simplify(expr); }

ResultExpr &expand(ResultExpr &&expr) { return expand(expr); }

ResultExpr &optimize(ResultExpr &&expr) { return optimize(expr); }

}  // namespace sequant
