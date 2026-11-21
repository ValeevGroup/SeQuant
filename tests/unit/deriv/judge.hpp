#ifndef SEQUANT_TESTS_UNIT_DERIV_JUDGE_HPP
#define SEQUANT_TESTS_UNIT_DERIV_JUDGE_HPP

#include "evaluator.hpp"

#include <SeQuant/core/expr_fwd.hpp>

#include <cstddef>
#include <string>
#include <string_view>

namespace sequant::deriv::test {

/// the number of operations of @p x : the size of the evaluator's enumeration
/// of the rearrangements of x (rearrangements_of()), the one
/// compare_tangents() counts; 1 for a Variable
std::size_t group_order(const ExprPtr& x, const Data& data = {});

/// the number of those that keep the block of @p x : each slot maps to a slot
/// of the same basis
std::size_t block_order(const ExprPtr& x, const Data& data = {});

/// the first term of @p r that core's simplify rejects or whose free indices
/// are not those of @p f , with the labels its projector carries, plus those
/// of @p x , and how; empty if none
std::string validity(const ExprPtr& f, const ExprPtr& x, const ExprPtr& r);

/// requires r == diff(f, x) to be @p expected (EquivalentTo; null: the judge
/// alone decides), to agree with compare_tangents() at every element (256
/// drawn where there are more) and to pass validity()
void check_rule(const ExprPtr& f, const ExprPtr& x, const ExprPtr& expected,
                const Data& data = {});
/// check_rule() of the deserialized texts; an empty @p expected is null
void check_rule(std::wstring_view f, std::wstring_view x,
                std::wstring_view expected, const Data& data = {});

}  // namespace sequant::deriv::test

#endif  // SEQUANT_TESTS_UNIT_DERIV_JUDGE_HPP
