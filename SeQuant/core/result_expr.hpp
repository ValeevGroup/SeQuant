#ifndef SEQUANT_RESULT_EXPR_HPP
#define SEQUANT_RESULT_EXPR_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <cassert>
#include <initializer_list>
#include <optional>
#include <string>

namespace sequant {

class Tensor;
class Variable;

class ResultExpr {
 public:
  using IndexContainer = container::svector<Index>;

  ResultExpr(const Tensor &tensor, ExprPtr expression);
  ResultExpr(const Variable &variable, ExprPtr expression);

  ResultExpr(const ResultExpr &other) = default;
  ResultExpr(ResultExpr &&other) = default;

  ResultExpr &operator=(const ResultExpr &other) = default;
  ResultExpr &operator=(ResultExpr &&other) = default;

  /// Assigns a new expression to this result
  ResultExpr &operator=(ExprPtr expression);

  bool has_label() const;
  const std::wstring &label() const;

  Symmetry symmetry() const;
  void set_symmetry(Symmetry symm);

  BraKetSymmetry braket_symmetry() const;
  void set_braket_symmetry(BraKetSymmetry symm);

  ParticleSymmetry particle_symmetry() const;
  void set_particle_symmetry(ParticleSymmetry symm);

  const IndexContainer &bra() const;
  const IndexContainer &ket() const;

  const ExprPtr &expression() const;
  ExprPtr &expression();

  template <typename Group>
  container::svector<Group> index_particle_grouping() const {
    container::svector<Group> groups;

    assert(m_braIndices.size() == m_ketIndices.size() &&
           "Not yet generalized to particle non-conserving results");

    groups.reserve(m_braIndices.size());

    // Note that the assumption is that indices are sorted
    // based on the particle they belong to and that bra and
    // ket indices are assigned to the same set of particles.
    for (std::size_t i = 0; i < m_braIndices.size(); ++i) {
      groups.emplace_back(
          std::initializer_list<Index>{m_braIndices.at(i), m_ketIndices.at(i)});
    }

    return groups;
  }

 private:
  ExprPtr m_expr;

  Symmetry m_symm = Symmetry::nonsymm;
  BraKetSymmetry m_bksymm = BraKetSymmetry::nonsymm;
  ParticleSymmetry m_psymm = ParticleSymmetry::nonsymm;
  IndexContainer m_braIndices;
  IndexContainer m_ketIndices;
  std::optional<std::wstring> m_label;
};

}  // namespace sequant

#endif  // SEQUANT_RESULT_EXPR_HPP
