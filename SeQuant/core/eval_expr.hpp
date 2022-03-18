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
class EvalExpr final {
 public:
  using hash_t = size_t;

  /**
   * Construct EvalExpr that goes into the leaf nodes
   * of the binary evaluation tree.
   */
  explicit EvalExpr(const Tensor&);

  /**
   * Construct EvalExpr that goes into the internal nodes
   * of the binary evaluation tree.
   */
  EvalExpr(const EvalExpr&, const EvalExpr&, EvalOp op);

  /**
   * Operation type of the object.
   */
  [[nodiscard]] EvalOp op() const;

  /**
   * Hash value of the object.
   */
  [[nodiscard]] hash_t hash() const;

  /**
   * Tensor expression stored by the object.
   */
  [[nodiscard]] const Tensor& tensor() const;

  /** Factor to scale tensor by. */
  [[nodiscard]] const Constant& scalar() const;

  ///
  /// annotation for TiledArray
  ///
  std::string annot() const {return annot_; }

  template <typename T = std::complex<double>>
  EvalExpr& operator*=(T fac) {
    scalar_ *= Constant{std::move(fac)};
    return *this;
  }

  template <typename T = std::complex<double>>
  void scale(T fac) {
    scalar_ = Constant{std::move(fac)};
  }

  friend inline bool operator==(const EvalExpr& lhs, const EvalExpr& rhs) {
    return lhs.hash() == rhs.hash() &&
           (lhs.tensor().to_latex() == rhs.tensor().to_latex());
  }

 private:
  using index_container_type = container::svector<Index>;
  using braket_type = std::pair<index_container_type, index_container_type>;

  EvalOp op_;

  hash_t hash_;

  Tensor tensor_;

  Constant scalar_{1};

  std::string annot_;

  /**
   * Infer the symmetry of the resulting tensor after summing two tensors.
   */
  static Symmetry infer_tensor_symmetry_sum(EvalExpr const&, EvalExpr const&);

  /**
   * Infer the symmetry of the resulting tensor from product of two tensors.
   */
  static Symmetry infer_tensor_symmetry_prod(EvalExpr const&, EvalExpr const&);

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
  static hash_t hash_braket(
      const decltype(std::declval<Tensor>().const_braket())&);

  /**
   * Hash a tensor based on its label and index space of braket indices.
   */
  static hash_t hash_terminal_tensor(const Tensor&);

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
                          EvalOp op);

  /**
   * Hash the topology of braket connectivity between a pair of tensor. The
   * positions of the connected indices in the corresponding brakets are hashed.
   */
  static hash_t hash_tensor_pair_topology(const Tensor&, const Tensor&);

  /**
   * Figure the target braket of resulting tensor from product of a pair of
   * tensors.
   *
   * @return Pair of bra and ket index containers where
   * the first is the target bra.
   */
  static braket_type target_braket_prod(const Tensor&, const Tensor&);

  ///
  /// Given an iterable of Index objects, generate a string annotation
  /// that can be used for TiledArray tensor expressions.
  /// Tensor-of-tensors also supported.
  template <typename Indices_t>
  static std::string braket_to_annot(Indices_t const& indices) {
    using ranges::views::transform;

    // make a comma-separated string out of an iterable of strings
    auto add_commas = [](auto const& strs) -> std::string {
      auto result = std::string{ranges::front(strs)};
      for (auto&& s: ranges::views::tail(strs))
        result += "," + s;
      return result;
    };

    container::vector<std::string> outer_labels{}, inner_labels{};
    for (auto&& idx: indices) {
      if (idx.has_proto_indices()) {
        inner_labels.emplace_back(idx.ascii_label());
        for (auto&& pidx : idx.proto_indices())
          outer_labels.emplace_back(pidx.ascii_label());
      } else {
        outer_labels.emplace_back(idx.ascii_label());
      }
    }

    if (inner_labels.empty())
      return add_commas(outer_labels);

    // support CSV methods
    ranges::sort(inner_labels);
    ranges::sort(outer_labels);
    ranges::actions::unique(outer_labels);
    return add_commas(outer_labels) + ";" + add_commas(inner_labels);
  }

};

}  // namespace sequant

#endif  // SEQUANT_EVAL_EXPR_HPP
