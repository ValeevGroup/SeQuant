#ifndef SEQUANT_EVAL_EXPR_HPP
#define SEQUANT_EVAL_EXPR_HPP

#include <SeQuant/core/tensor.hpp>

namespace sequant {

/**
 * Type of evaluation resulting into EvalExpr.
 */
enum class EvalOp {
  /** Identity type evaluation. */
  Id,
  /** Sum type evaluation. */
  Sum,
  /** Product type evaluation. */
  Prod,
  /** Symmetrization type evaluation. */
  Symm,
  /** Antisymmetrization type evaluation. */
  Antisymm,
};  // EvalOp

/**
 * An object to be stored in the binary evaluation tree nodes.
 *
 * @author Bimal Gaudel
 * @version 03 Dec, 2020
 */
class EvalExpr {
 public:
  using hash_t = size_t;

  EvalExpr() = delete;

  /**
   * Construct EvalExpr that goes into the leaf nodes
   * of the binary evaluation tree.
   */
  explicit EvalExpr(const Tensor&) noexcept;

  /**
   * Construct EvalExpr that goes into the internal nodes
   * of the binary evaluation tree.
   */
  EvalExpr(const EvalExpr&, const EvalExpr&, EvalOp op) noexcept;

#ifdef DEBUG_EVAL_EXPR
  EvalExpr(EvalExpr&&) noexcept;

  EvalExpr& operator=(EvalExpr&&) noexcept;

  EvalExpr(EvalExpr const&) noexcept;

  EvalExpr& operator=(EvalExpr const&) noexcept;

  ~EvalExpr() noexcept;
#endif

  /**
   * The unique id of this object.
   */
  [[nodiscard]] size_t id() const noexcept;

  /**
   * The store Tensor's label plus the id.
   */
  [[nodiscard]] std::string label(std::string_view sep = "_") const noexcept;

  /**
   * Operation type of the object.
   */
  [[nodiscard]] EvalOp op() const noexcept;

  /**
   * Hash value of the object.
   */
  [[nodiscard]] hash_t hash_value() const noexcept;

  /**
   * Tensor expression stored by the object.
   */
  [[nodiscard]] const Tensor& tensor() const noexcept;

  /** Factor to scale tensor by. */
  [[nodiscard]] const Constant& scalar() const noexcept;

  /** Scale the scalar prefactor by @c fac. */
  template <typename T = std::complex<double>>
  EvalExpr& operator*=(T fac) noexcept {
    scalar_ *= Constant{std::move(fac)};
    return *this;
  }

  /** Set the scalar prefactor to @c fac. */
  template <typename T = std::complex<double>>
  void scale(T fac) noexcept {
    scalar_ = Constant{std::move(fac)};
  }

  friend inline bool operator==(const EvalExpr& lhs,
                                const EvalExpr& rhs) noexcept {
    return lhs.hash_value() == rhs.hash_value() &&
           (lhs.tensor().to_latex() == rhs.tensor().to_latex());
  }

 private:
  using index_container_type = container::svector<Index>;
  using braket_type = std::pair<index_container_type, index_container_type>;

  size_t id_;

  EvalOp op_;

  hash_t hash_;

  Tensor tensor_;

  Constant scalar_{1};

  static size_t global_id;

  /**
   * Infer the symmetry of the resulting tensor after summing two tensors.
   */
  static Symmetry infer_tensor_symmetry_sum(EvalExpr const&,
                                            EvalExpr const&) noexcept;

  /**
   * Infer the symmetry of the resulting tensor from product of two tensors.
   */
  static Symmetry infer_tensor_symmetry_prod(EvalExpr const&,
                                             EvalExpr const&) noexcept;

  /**
   * Infers the particle symmetry of a tensor.
   *
   * @param s The tensor symmetry of the tensor.
   */
  static ParticleSymmetry infer_particle_symmetry(Symmetry s) noexcept;

  /**
   * Infers the braket symmetry of a tensor based on the default sequant
   * context.
   */
  static BraKetSymmetry infer_braket_symmetry() noexcept;

  /**
   * Combined hash of the index spaces of the indices in a braket.
   */
  static hash_t hash_braket(
      const decltype(std::declval<Tensor>().const_braket())&) noexcept;

  /**
   * Hash a tensor based on its label and index space of braket indices.
   */
  static hash_t hash_terminal_tensor(const Tensor&) noexcept;

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
  static hash_t hash_imed(const EvalExpr& expr1, const EvalExpr& expr2,
                          EvalOp op) noexcept;

  /**
   * Hash the topology of braket connectivity between a pair of tensor. The
   * positions of the connected indices in the corresponding brakets are hashed.
   */
  static hash_t hash_tensor_pair_topology(const Tensor&,
                                          const Tensor&) noexcept;

  /**
   * Figure the target braket of resulting tensor from product of a pair of
   * tensors.
   *
   * @return Pair of bra and ket index containers where
   * the first is the target bra.
   */
  static braket_type target_braket_prod(const Tensor&, const Tensor&) noexcept;
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_EXPR_HPP
