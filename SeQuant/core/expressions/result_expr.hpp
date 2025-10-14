#ifndef SEQUANT_EXPRESSIONS_RESULT_EXPR_HPP
#define SEQUANT_EXPRESSIONS_RESULT_EXPR_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <initializer_list>
#include <optional>
#include <string>
#include <type_traits>

namespace sequant {

/// Represents an expression containing a left-hand-side of the form
/// `<lhs> = <expression>`
/// Thus, objects of this class are used to keep track of the result
/// of a given expression as well as the expression itself.
///
/// This is particularly important when it comes to tracking the exact
/// pairing (and symmetry) of external indices in cases the given
/// expression computes a tensorial property.
///
/// Importantly, it is possible to track the result without giving
/// an explicit name (label) to it.
class ResultExpr {
 public:
  using IndexContainer = container::svector<Index>;

  struct ResultCmp {
    bool operator()(const ResultExpr &lhs, const ResultExpr &rhs) const;
  };

  ResultExpr(const Tensor &tensor, ExprPtr expression);
  ResultExpr(const Variable &variable, ExprPtr expression);

  template <typename Range1, typename Range2, typename Range3>
  ResultExpr(bra<Range1> &&bra, ket<Range2> &&ket, aux<Range3> &&aux,
             Symmetry symm, BraKetSymmetry braket_symm,
             ColumnSymmetry column_symm, std::optional<std::wstring> label,
             ExprPtr expr)
      : ResultExpr(std::forward<Range1>(bra), std::forward<Range2>(ket),
                   std::forward<Range3>(aux), symm, braket_symm, column_symm,
                   std::move(label), std::move(expr)) {}

  ResultExpr(const ResultExpr &other) = default;
  ResultExpr(ResultExpr &&other) = default;

  ResultExpr &operator=(const ResultExpr &other) = default;
  ResultExpr &operator=(ResultExpr &&other) = default;

  /// Assigns a new expression to this result
  ResultExpr &operator=(ExprPtr expression);

  bool has_label() const;
  const std::wstring &label() const;

  void set_label(std::wstring label);

  Symmetry symmetry() const;
  void set_symmetry(Symmetry symm);

  BraKetSymmetry braket_symmetry() const;
  void set_braket_symmetry(BraKetSymmetry symm);

  ColumnSymmetry column_symmetry() const;
  void set_column_symmetry(ColumnSymmetry symm);

  const IndexContainer &bra() const;
  const IndexContainer &ket() const;
  const IndexContainer &aux() const;

  const ExprPtr &expression() const;
  ExprPtr &expression();

  ResultExpr clone() const;

  /// Obtains the exact grouping (pairing) of indices in the result. Typically,
  /// this represents particle-assignments (i.e. indices in the same group are
  /// associated with the same particle in the underlying theory).
  ///
  /// @tparam Group The type of the object to represent an index group. Must be
  /// constructible from an initializer_list<Index> or from a set of two
  /// indices. This allows for pair/tuple and vector-like types.
  template <typename Group>
  container::svector<Group> index_particle_grouping() const {
    container::svector<Group> groups;

    SEQUANT_ASSERT(m_braIndices.size() == m_ketIndices.size() &&
                   "Not yet generalized to particle non-conserving results");
    SEQUANT_ASSERT(m_auxIndices.empty() &&
                   "Not yet clear how auxiliary indices should be handled");

    groups.reserve(m_braIndices.size());

    // Note that the assumption is that indices are sorted
    // based on the particle they belong to and that bra and
    // ket indices are assigned to the same set of particles.
    for (std::size_t i = 0; i < m_braIndices.size(); ++i) {
      if constexpr (std::is_constructible_v<Group,
                                            std::initializer_list<Index>>) {
        groups.emplace_back(std::initializer_list<Index>{m_braIndices.at(i),
                                                         m_ketIndices.at(i)});
      } else {
        static_assert(
            std::is_constructible_v<Group, Index, Index>,
            "Group is expected to be constructible from two indices or from an "
            "initializer_list of indices");
        groups.emplace_back(m_braIndices.at(i), m_ketIndices.at(i));
      }
    }

    return groups;
  }

  Tensor result_as_tensor(std::wstring default_label = L"Unnamed") const {
    return Tensor(m_label.has_value() ? m_label.value() : default_label,
                  sequant::bra(m_braIndices), sequant::ket(m_ketIndices),
                  sequant::aux(m_auxIndices), m_symm, m_bksymm, m_csymm);
  }

  Variable result_as_variable(std::wstring default_label = L"Unnamed") const {
    SEQUANT_ASSERT(!produces_tensor());
    return Variable(m_label.has_value() ? m_label.value() : default_label);
  }

  bool produces_tensor() const {
    return !m_braIndices.empty() || !m_ketIndices.empty() ||
           !m_auxIndices.empty();
  }

  bool operator==(const ResultExpr &other) const = default;

 private:
  ExprPtr m_expr;

  template <typename Range1, typename Range2, typename Range3>
    requires(!std::is_constructible_v<IndexContainer, Range1> ||
             !std::is_constructible_v<IndexContainer, Range2> ||
             !std::is_constructible_v<IndexContainer, Range3>)
  ResultExpr(Range1 &&bra, Range2 &&ket, Range3 &&aux, Symmetry symm,
             BraKetSymmetry braket_symm, ColumnSymmetry column_symm,
             std::optional<std::wstring> label, ExprPtr expression)
      : ResultExpr(IndexContainer(bra.begin(), bra.end()),
                   IndexContainer(ket.begin(), ket.end()),
                   IndexContainer(aux.begin(), aux.end()), symm, braket_symm,
                   column_symm, std::move(label), std::move(expression)) {}

  ResultExpr(IndexContainer bra, IndexContainer ket, IndexContainer aux,
             Symmetry symm, BraKetSymmetry braket_symm,
             ColumnSymmetry column_symm, std::optional<std::wstring> label,
             ExprPtr expression);

  Symmetry m_symm = Symmetry::Nonsymm;
  BraKetSymmetry m_bksymm = BraKetSymmetry::Nonsymm;
  ColumnSymmetry m_csymm = ColumnSymmetry::Nonsymm;
  IndexContainer m_braIndices;
  IndexContainer m_ketIndices;
  IndexContainer m_auxIndices;
  std::optional<std::wstring> m_label;
};

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_RESULT_EXPR_HPP
