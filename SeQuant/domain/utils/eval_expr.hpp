#ifndef SEQUANT_UTILS_EVAL_EXPR_HPP
#define SEQUANT_UTILS_EVAL_EXPR_HPP

#include <SeQuant/core/tensor.hpp>

namespace sequant::utils {

/**
 * An object to be stored in the binary evaluation tree nodes.
 *
 * @author Bimal Gaudel
 * @version 03 Dec, 2020
 */
class eval_expr final {
 public:
  /**
   * Type of evaluation resulting into eval_expr.
   */
  enum class eval_op {
    /** Identity type evaluation. */
    Id,
    /** Sum type evaluation. */
    Sum,
    /** Product type evaluation. */
    Prod,
  };  // eval_op

  /**
   * Construct eval_expr that goes into the leaf nodes
   * of the binary evaluation tree.
   */
  explicit eval_expr(Tensor);

  /**
   * Construct eval_expr that goes into the internal nodes
   * of the binary evaluation tree.
   */
  eval_expr(const eval_expr&, const eval_expr&);

  /**
   * Operation type of the object.
   */
  [[nodiscard]] eval_op op() const;

  /**
   * Hash value of the object.
   */
  [[nodiscard]] size_t hash() const;

  /**
   * Tensor expression stored by the object.
   */
  [[nodiscard]] const Tensor& tensor() const;

  /** Factor to scale tensor by. */
  [[nodiscard]] const Constant& scalar() const;

  /** Phase of the tensor resulting from canonicalization. */
  [[nodiscard]] const Constant& phase() const;

  /** Set the scalar. **/
  void scale(const Constant&);

  friend inline bool operator==(const eval_expr& lhs, const eval_expr& rhs) {
    return lhs.hash() == rhs.hash() &&
           (lhs.tensor().to_latex() == rhs.tensor().to_latex());
  }

 private:
  using index_container_type = container::svector<Index>;
  using braket_type = std::pair<index_container_type, index_container_type>;

  eval_op op_;

  size_t hash_;

  Tensor tensor_;

  Constant scalar_{1};

  Constant phase_{1};

  /**
   * Make an intermediate tensor from two expressions.
   */
  static Tensor make_imed_expr(const eval_expr& expr1, const eval_expr& expr2,
                                eval_op op);

  /**
   * Figure out the type of evaluation between a pair of expressions.
   */
  static eval_op infer_eval_op(const Tensor&, const Tensor&);

  /**
   * Infer the symmetry of the resulting tensor after summing two tensors.
   */
  static Symmetry infer_tensor_symmetry_sum(const Tensor&, const Tensor&);

  /**
   * Infer the symmetry of the resulting tensor from product of two tensors.
   */
  static Symmetry infer_tensor_symmetry_prod(const Tensor&, const Tensor&);

  /**
   * Infers the particle symmetry of a tensor.
   *
   * @param s The tensor symmetry of the tensor.
   */
  static ParticleSymmetry infer_particle_symmetry(Symmetry s);

  /**
   * Infers the braket symmetry of a tensor based on the default sequant
   * context.
   */
  static BraKetSymmetry infer_braket_symmetry();

  /**
   * Combined hash of the index spaces of the indices in a braket.
   */
  static size_t hash_braket(
      const decltype(std::declval<Tensor>().const_braket())&);

  /**
   * Hash a tensor based on its label and index space of braket indices.
   */
  static size_t hash_terminal_tensor(const Tensor&);

  /**
   * Hash an intermediate tensor based on the topology of braket connectivity.
   *
   * @param expr1 A sequant expression.
   * @param expr2 A sequant expression.
   *
   * @param hash1 Hash value of expr1.
   * @param hash2 Hash value of expr2.
   *
   * @param op Type of evaluation operation between @c expr1 and @c expr2.
   */
  static size_t hash_imed(const eval_expr& expr1, const eval_expr& expr2,
                          eval_op op);

  /**
   * Hash tensor pair(T1, T2) based on the tensor labels and the topology of
   * braket connection between them.
   *
   * @note result of pair(T1, T2) may or may not equal that of pair(T2, T1).
   */
  static size_t hash_tensor_pair(const Tensor&, const Tensor&);

  /**
   * Hash the topology of braket connectivity between a pair of tensor. The
   * positions of the connected indices in the corresponding brakets are hashed.
   */
  static size_t hash_tensor_pair_topology(const Tensor&, const Tensor&);

  /**
   * Figure the target braket of resulting tensor from sum between a pair of
   * tensors.
   *
   * @return Pair of bra and ket index containers where
   * the first is the target bra.
   */
  static braket_type target_braket_sum(const Tensor&, const Tensor&);

  /**
   * Figure the target braket of resulting tensor from product of a pair of
   * tensors.
   *
   * @return Pair of bra and ket index containers where
   * the first is the target bra.
   */
  static braket_type target_braket_prod(const Tensor&, const Tensor&);
};  // eval_expr

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_EVAL_EXPR_HPP
