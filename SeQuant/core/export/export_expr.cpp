#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/export/compute_selection.hpp>
#include <SeQuant/core/export/export_expr.hpp>

#include <atomic>
#include <cstddef>

namespace sequant {

std::size_t detail::next_export_expr_id() {
  static std::atomic<std::size_t> next_id = 0;

  std::size_t id = next_id.load();

  // This ensures that we are updating next_id with the next id while
  // also ensuring that no other thread is currently producing the same
  // id that this one is doing.
  while (!next_id.compare_exchange_weak(id, id + 1)) {
  }

  return id;
}

ExportExpr::ExportExpr(const EvalExpr &other) : EvalExpr(other) {}

std::size_t ExportExpr::id() const { return m_id; }

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
