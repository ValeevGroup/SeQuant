#ifndef SEQUANT_TESTS_UNIT_DERIV_IDENTITIES_HPP
#define SEQUANT_TESTS_UNIT_DERIV_IDENTITIES_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include "evaluator.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::deriv::test {

// the family

/// a model of the family
struct Model {
  bool csv = false;
  /// closed-shell spin-traced
  bool traced = false;
  /// open-shell spin-traced
  bool open_shell = false;
};

std::string name_of(const Model& model);

/// CC{2}: canonical and CSV, spin-orbital and closed-shell traced, and
/// canonical open-shell; under SEQUANT_SKIP_LONG_TESTS the spin-orbital ones
/// only. CSV open-shell is left out until core is fixed: core's
/// open_shell_CC_spintrace under CSV gives a protoindex the spin of the slot
/// that carries it, so the index of a mixed-spin pair holds two spins within a
/// term
container::svector<Model> family_models();

/// the equations of a model and the variables to differentiate them by
struct Family {
  Model model;
  std::vector<std::pair<std::wstring, ExprPtr>> equations;
  /// every leaves_as_variables() spelling and the generated spellings
  Expressions variables;
};

/// calls @p f with the family of @p model , made once per process and used
/// under its context
void with_family(const Model& model,
                 const std::function<void(const Family&)>& f);

/// the equation of @p family named @p name ; it must exist
ExprPtr equation_named(const Family& family, std::wstring_view name);

// the inputs

/// the elements J the evaluator reads per (equation, variable) pair
inline constexpr std::size_t samples_per_pair = 1;

/// the draws I1 takes per pair on its second try
inline constexpr std::size_t retry_draws = 64;

/// the variables other than amplitudes a model is checked by
inline constexpr std::size_t sampled_variables = 4;

/// a hash of @p text that is the same on every platform
std::uint64_t stable_hash(std::string_view text);

/// whether @p x is an amplitude: t or λ, perturbed or conjugated
bool is_amplitude(const ExprPtr& x);

/// the variables the identities run by, and those diff rejects
struct CheckedVariables {
  /// every amplitude, and the sampled_variables other variables (integrals,
  /// overlaps) of least stable_hash()
  Expressions variables;
  /// the variables diff rejects with some equation of the family, by the text
  /// of the first rejection
  container::map<std::string, Expressions> rejected;
};

/// the variables of @p family diff accepts with every equation, and those it
/// rejects; anything else diff throws propagates
CheckedVariables checked_variables(const Family& family);

/// @p x with the label @p label , its symmetries and indices kept; the
/// conjugate spelling stays conjugate
ExprPtr relabeled(const ExprPtr& x, std::wstring_view label);

/// @p x with each index and protoindex renamed to the ordinal @p by higher
ExprPtr renamed_apart(const ExprPtr& x, std::size_t by);

/// @p e without the factors Â, Ŝ whose slots are all among @p indices
ExprPtr without_projectors_over(const ExprPtr& e,
                                const container::set<Index>& indices);

// the identities

/// I1, definition: compare_tangents() of @p r = diff(f, x), at
/// samples_per_pair nonzero elements drawn with @p seed ; where r is not 0
/// and no draw finds a nonzero element, in retry_draws draws over three
/// occupied values, then at every element
Comparison definition(const ExprPtr& f, const ExprPtr& x, const ExprPtr& r,
                      std::uint64_t seed);

/// I3, Taylor about x = 0: f with every leaf of the block of x replaced by h
/// (x⁺ by h⁺) equals Σ (1/(n! m! |B|^(n+m))) ⟨diff(f, x…x x⁺…x⁺) at x = 0,
/// h…h h⁺…h⁺⟩ over n derivatives by x and m by x⁺ up to f's degree in each,
/// exactly; h is x relabelled ζ, |B| is @p order , 0: block_order(x). Where
/// x⁺ holds the elements of x, m is 0. The arrays are compared at
/// samples_per_pair elements drawn with @p seed
Comparison taylor(const ExprPtr& f, const ExprPtr& x, std::uint64_t seed,
                  std::size_t order = 0);

/// I4, mixed partials: diff(diff(f, {x}), {z}) and
/// @p scale diff(diff(f, {z}), {x}) are the same array, z being @p y renamed
/// apart from x; @p scale other than 1 is a wrong result
Comparison schwarz(const ExprPtr& f, const ExprPtr& x, const ExprPtr& y,
                   std::uint64_t seed, Constant::scalar_type scale = 1);

/// the first amplitude after @p x in @p amplitudes , cyclically, by which
/// @p r = diff(f, {x}) has a derivative not identically zero, spelled as
/// schwarz() spells it; the one right after x if there is none
const ExprPtr& schwarz_partner(const ExprPtr& r, const ExprPtr& x,
                               const Expressions& amplitudes);

/// one comparison of I2: an equation SeQuant generates against the derivative
/// that gives it
struct OwnEquation {
  std::wstring name;
  /// the generated equation, as compared
  ExprPtr generated;
  /// the derivative, as compared
  ExprPtr expected;
  container::set<Index> free;
};

/// I2, SeQuant's own equations of @p family : t()[k] against diff(ℒ, λ_k),
/// λ()[k] against diff(ℒ, t_k) (closed-shell: after SeQuant's biorthogonal
/// transform), and tʼ[k], λʼ[k] without their h¹ and ω terms against the
/// ε-derivative of t[k], λ[k] at t + ε t¹ (λ + ε λ¹), at ε = 0
container::svector<OwnEquation> own_equations(const Family& family);

}  // namespace sequant::deriv::test

#endif  // SEQUANT_TESTS_UNIT_DERIV_IDENTITIES_HPP
