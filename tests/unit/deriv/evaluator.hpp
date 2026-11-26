#ifndef SEQUANT_TESTS_UNIT_DERIV_EVALUATOR_HPP
#define SEQUANT_TESTS_UNIT_DERIV_EVALUATOR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <utility>

namespace sequant::deriv::test {

using Expressions = container::svector<ExprPtr>;

/// the summands of @p e , nested sums flattened
Expressions summands_of(const ExprPtr& e);

/// the factors of @p term , nested products flattened, and its scalar
Expressions factors_of(const ExprPtr& term, Constant::scalar_type& scalar);

/// whether @p e is a projector, Â or Ŝ
bool is_projector(const Expr& e);

/// the indices of a product, by core's get_used_indices_with_counts
struct ProductIndices {
  /// those outside exactly two non-proto slots (a projector's count like any
  /// other), and the protoindices of those
  container::set<Index> free;
  /// those a projector carries
  container::set<Index> projected;
};

ProductIndices indices_of_product(const ExprPtr& product);

/// the products of @p term : itself, or the terms of its expansion where a
/// factor is a Sum
Expressions products_of(const ExprPtr& term);

/// the external indices of @p f : those free in its first product and those a
/// projector of it carries
container::set<Index> externals_of(const ExprPtr& f);

/// the slots of the leaf @p x and their protoindices
container::set<Index> indices_of_leaf(const ExprPtr& x);

/// the symbol of a leaf, apart from its indices: a Variable's label and
/// conjugation, a Tensor's label and ranks
std::wstring symbol_of(const Expr& leaf);

/// a new name for each of @p indices that @p name gives one, those without
/// protoindices first; an index with protoindices also takes the new names of
/// its protoindices
container::map<Index, Index> names_for(
    const container::set<Index>& indices,
    const std::function<std::optional<Index>(const Index&)>& name);

/// index labels, each with the value it takes
using Values = container::map<std::wstring, int>;

/// index labels, each with the values it runs over
using Domains =
    container::svector<std::pair<std::wstring, container::svector<int>>>;

/// adds the label of @p i to @p domains , unless @p values or @p domains hold
/// it: @p occupied values per occupied and three per virtual elementary space
/// its space spans, an elementary space being one type bit and, if @p
/// quantum_numbers , one quantum-number (spin) bit, so that a subspace runs
/// over fewer values than its superspace and spin cases differ
void add_domain(const Index& i, const Values& values, Domains& domains,
                bool quantum_numbers = true, int occupied = 2);

/// what the data of the evaluator knows beyond each leaf
struct Data {
  /// the labels of the real symbols (a tensor's without the adjoint marker):
  /// their data is real
  container::set<std::wstring> real = {};
  /// Â and Ŝ are the normalized projectors they stand for, applied to the
  /// rest of their term over their labels, which must be external labels;
  /// otherwise they are data like any other tensor
  bool projectors = false;
  /// the values an occupied elementary space runs over: three give a triple
  /// of distinct values
  int occupied = 2;
};

/// the element of @p e at @p values of some index labels, summed over the
/// values of the others (a factor that is a Sum is expanded first). Spin cases
/// are distinguished where the indices of the input differ in spin. An index
/// with protoindices is addressed by its value and those of its protoindices,
/// sorted if symmetric. A δ, and an overlap s between indices of one
/// orthonormal basis whose protoindices take equal values, is 1 where its two
/// indices agree, else 0; an overlap between two bases is data. The data of a
/// tensor carries exactly the symmetry it declares: (anti)symmetric bra and
/// ket, symmetric columns, bra-ket symmetric (real data, as core's adjoint
/// assumes) or Hermitian (`Conjugate`: complex data, the exchanged element its
/// conjugate); t⁺ is the conjugate of t at the exchanged slots, x^* of x. A
/// Variable and the data of every other tensor are complex unless
/// @p data says the symbol is real. A Power raises the element of its base,
/// a TensorPower element-wise
std::complex<double> element_of(const ExprPtr& e, const Values& values,
                                const Data& data = {});

/// which elements a comparison reads: every one, or, where there are more than
/// @c count , elements drawn at random from a generator seeded with @c seed
/// until @c count of them are nonzero (an element that vanishes by symmetry
/// does not count) or @c draws are drawn. A drawn label repeats,
/// half the time, the value of an earlier label over its range in another of
/// @c groups , else takes one no earlier label of its own group has
struct Sampling {
  /// 0: every element
  std::size_t count = 0;
  std::uint64_t seed = 0;
  /// the most draws; 0: eight times count
  std::size_t draws = 0;
  /// the external indices by the tensor they belong to (the function and each
  /// variable); a label of none is in the first
  container::svector<container::set<Index>> groups = {};
};

/// the comparison of two arrays over the labels of some indices
struct Comparison {
  /// how many elements were compared, and how many of them are nonzero in
  /// either array
  std::size_t elements = 0, nonzero = 0;
  /// the first element that differs; empty if none does
  std::string difference;
  /// whether every element was read, not a sample
  bool every_element = false;
};

/// @p a and @p b compared as arrays of the labels of @p free and their
/// protoindices, every other label summed over; elements equal to 1e-12;
/// the elements @p sampling reads
Comparison compare_arrays(const ExprPtr& a, const ExprPtr& b,
                          const container::set<Index>& free,
                          const Data& data = {}, const Sampling& sampling = {});

/// df/dx against the evaluator: for each value of the free labels of @p f
/// and each element J of @p x , @p r at J, the projectors applied, times the
/// number of distinct arrangements J' the operations of x make of J equals
/// the number of operations of x (its rearrangements_of(), real as @p data
/// says) times the derivative of @p f , its projectors applied, when the
/// element J is raised by ε, its partners with their signs (exact, over dual
/// numbers; Wirtinger's: the conjugate spelling of x, x^* or t⁺, and for a
/// Hermitian x the elements at the other exchange do not count, while for a
/// real x the conjugate spelling is x itself and counts); that @p r has the
/// symmetry of x is checked, r(J') =
/// sign(J') r(J), at each J' some values of the labels of x read (a slot
/// carries its protoindices along, which a label need not follow); the values
/// of the labels of f and x that @p sampling reads
Comparison compare_tangents(const ExprPtr& f, const ExprPtr& r,
                            const ExprPtr& x, const Data& data = {},
                            const Sampling& sampling = {});

/// @p e serialized, or as LaTeX where a leaf has no serialization
std::string text_of(const ExprPtr& e);

/// an order of the labels of a tensor that its declared symmetry reads as
/// the same element
struct Rearrangement {
  /// the slots, bra, ket and aux in order
  container::svector<Index> slots;
  /// the protoindices of the slots that are no slots, each read in place of
  /// the one at its position in the first rearrangement
  container::svector<Index> protoindices;
  int sign = 1;
};

/// the slots of @p t , bra, ket and aux in order, in every order its declared
/// symmetry reads as the same element, with the sign of each: the bra and the
/// ket permuted apart for an (anti)symmetric tensor, the paired columns and
/// the unpaired slots of a column-symmetric one; where bra and ket ranks are
/// equal and nonzero, bra and ket exchanged too for a bra-ket symmetric one,
/// and for a Hermitian one that is @p real , whose exchanged element is the
/// same real number. Each with every order of the protoindices that are no
/// slots of @p t that keeps the protoindex list of every slot, as a set where
/// it is symmetric. The slots and protoindices as they stand come first
container::svector<Rearrangement> rearrangements_of(const AbstractTensor& t,
                                                    bool real);

}  // namespace sequant::deriv::test

#endif  // SEQUANT_TESTS_UNIT_DERIV_EVALUATOR_HPP
