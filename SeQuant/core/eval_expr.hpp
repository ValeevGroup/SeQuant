#ifndef SEQUANT_EVAL_EXPR_HPP
#define SEQUANT_EVAL_EXPR_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/tensor.hpp>

namespace sequant {

enum class EvalOp { Id, Sum, Prod };

enum class EvalResult { Tensor, Constant };

class EvalExpr {
 private:
  EvalOp op_type_;
  EvalResult result_type_;
  size_t hash_value_;
  size_t id_;
  ExprPtr expr_;
  bool tot_;

  static size_t global_id_;

 public:
  explicit EvalExpr(Tensor const& tnsr);

  explicit EvalExpr(Constant const& c);

  EvalExpr(EvalExpr const& left, EvalExpr const& right, EvalOp op);

  [[nodiscard]] EvalOp op_type() const noexcept;

  [[nodiscard]] EvalResult result_type() const noexcept;

  [[nodiscard]] size_t hash_value() const noexcept;

  [[nodiscard]] size_t id() const noexcept;

  [[nodiscard]] ExprPtr expr() const noexcept;

  [[nodiscard]] bool tot() const noexcept;

  [[nodiscard]] std::wstring to_latex() const noexcept;

  [[nodiscard]] sequant::Tensor const& as_tensor() const noexcept(false);

  [[nodiscard]] sequant::Constant const& as_constant() const noexcept(false);
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_EXPR_HPP
