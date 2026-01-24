#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/export/compute_selection.hpp>
#include <SeQuant/core/export/export_expr.hpp>
#include <SeQuant/core/utility/atomic.hpp>

#include <atomic>
#include <cstddef>

namespace sequant {

std::size_t detail::next_export_expr_id() {
  static std::atomic<std::size_t> next_id = 0;

  return fetch_and_increment(next_id);
}

ExportExpr::ExportExpr(const EvalExpr &other) : EvalExpr(other) {}

std::size_t ExportExpr::id() const { return m_id; }

void ExportExpr::set_id(std::size_t id) { m_id = id; }

ComputeSelection ExportExpr::compute_selection() const { return m_selection; }

void ExportExpr::set_compute_selection(ComputeSelection selection) {
  m_selection = selection;
}

void ExportExpr::select_left() { m_selection |= ComputeSelection::Left; }

void ExportExpr::select_right() { m_selection |= ComputeSelection::Right; }

void ExportExpr::deselect_left() { m_selection &= ~ComputeSelection::Left; }

void ExportExpr::deselect_right() { m_selection &= ~ComputeSelection::Right; }

void ExportExpr::set_expr(ExprPtr expr) { expr_ = std::move(expr); }

bool ExportExpr::operator==(const ExportExpr &other) const {
  return m_id == other.m_id;
}

}  // namespace sequant
