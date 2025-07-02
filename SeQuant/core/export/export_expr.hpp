#ifndef SEQUANT_CORE_EXPORT_EXPORT_EXPR_HPP
#define SEQUANT_CORE_EXPORT_EXPORT_EXPR_HPP

#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/export/compute_selection.hpp>

#include <cstddef>

namespace sequant {

namespace detail {
std::size_t next_export_expr_id();
}

/// Dedicated Node data class to be used for exporting expressions (i.e. code
/// generation)
class ExportExpr : public EvalExpr {
 public:
  using EvalExpr::EvalExpr;

  ExportExpr(const EvalExpr &other);

  /// @returns The ID of this object. IDs uniquely determine the object identity
  [[nodiscard]] std::size_t id() const;
  /// Sets the ID of this object
  void set_id(std::size_t id);

  /// @returns The ComputeSelection for this object
  [[nodiscard]] ComputeSelection compute_selection() const;
  /// Sets the ComputeSelection for this object
  void set_compute_selection(ComputeSelection selection);
  /// Modifies the ComputeSelection to have the left subtree selected
  void select_left();
  /// Modifies the ComputeSelection to have the right subtree selected
  void deselect_left();
  /// Modifies the ComputeSelection to deselect the left subtree
  void select_right();
  /// Modifies the ComputeSelection to deselect the right subtree
  void deselect_right();

  /// Sets the expression stored by this object
  void set_expr(ExprPtr expr);

  bool operator==(const ExportExpr &other) const;

 private:
  ComputeSelection m_selection = ComputeSelection::Both;
  std::size_t m_id = detail::next_export_expr_id();
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_EXPORT_EXPR_HPP
