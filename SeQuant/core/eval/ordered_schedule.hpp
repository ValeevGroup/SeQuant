#ifndef SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
#define SEQUANT_EVAL_ORDERED_SCHEDULE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <utility>
#include <variant>

namespace sequant::eval {

///
/// \brief Task 2 of the ordered-scope batched-eval design (SP2): the \c
/// OrderedSchedule IR -- an ORDERED tree of loop blocks and build steps.
/// Purely data + a structural \c well_formed check; no sequencer/executor
/// here (that is a later task).
///
/// \details Unlike \c ScopeSchedule (\c scope_schedule.hpp), whose \c
/// ScopeNode groups all values homed at a node into one unordered \c
/// homed_values bag, \c OrderedSchedule threads builds and child loop blocks
/// through a single ORDERED sequence (\c ScopeBlock::steps): a value built
/// AFTER a child block in that sequence may read the child block's
/// accumulated output, so relative order among steps is load-bearing, not
/// incidental.
///

///
/// \brief What happens to a value at the close of its home \c ScopeBlock.
///
enum class OutputKind {
  Transient,          //!< the value's home IS this block; nothing carries out
  AccumulateSum,      //!< reduction: summed into an outer accumulator
  AccumulateScatter,  //!< loop-carried: scattered into a disjoint outer slice
};

///
/// \brief Build this value here: one contraction of already-cached operands.
///
struct BuildStep {
  std::size_t value_id;
};

struct ScopeBlock;  // fwd; ScopeBlock::steps (below) holds container::vector
                    // <Step>, so Step must be at least forward-declared
                    // before ScopeBlock -- see Step's own doc comment (right
                    // after ScopeBlock) for why Step's FULL definition has to
                    // come AFTER ScopeBlock's, not before it.
struct Step;        // fwd; see immediately below.

///
/// \brief One loop block: a batch loop over \c axis (sentinel/default on the
/// root block, which sits outside every loop), containing an ORDERED
/// sequence of build steps and nested child blocks, plus the set of values
/// that leave this block (and how) when it closes.
///
struct ScopeBlock {
  Index axis{};     //!< the loop axis; default (sentinel) on the root block.
  int ordinal = 0;  //!< disambiguates recurring sibling blocks realizing the
                    //!< same axis (e.g. two disjoint passes over the same
                    //!< axis TYPE at one nesting level).
  BatchModeType kind =
      BatchModeType::Contracted;    //!< Contracted (accumulate on block exit)
                                    //!< or External (scatter on block exit);
                                    //!< meaningless on the root.
  container::vector<Step> steps{};  //!< ORDERED: build-or-child-block,
                                    //!< interleaved. \c container::vector
                                    //!< (not \c svector) because \c Step is
                                    //!< only forward-declared here -- see \c
                                    //!< Step's doc comment.
  container::svector<std::pair<std::size_t, OutputKind>>
      outputs{};  //!< value_id -> how it leaves this block on close.
};

///
/// \brief One ORDERED step of a \c ScopeBlock: build a value at this scope,
/// or enter a nested child loop block.
///
/// \details The design brief's target shape is a bare alias, \c using Step =
/// std::variant<BuildStep, ScopeBlock>. That is not directly expressible:
/// \c ScopeBlock::steps must hold a sequence of \c Step (so a build
/// interleaves with child blocks in one ORDERED list, per the class doc
/// above), which means \c Step has to be at least forward-declared before \c
/// ScopeBlock; but \c std::variant requires every alternative type COMPLETE
/// at the point the variant specialization is instantiated (unlike \c
/// container::vector / \c std::vector, which the C++17 library tolerates
/// holding an incomplete element type up to first use -- the same relaxation
/// \c ScopeNode::children in \c scope_schedule.hpp relies on for ITS
/// self-reference). A bare \c using Step = std::variant<BuildStep,
/// ScopeBlock> declared before \c ScopeBlock therefore cannot compile
/// (\c ScopeBlock incomplete there), and a type alias cannot be
/// forward-declared separately from its definition (no "using Step;"
/// forward declaration exists in C++) to defer it to after \c ScopeBlock
/// either.
///
/// The fix: make \c Step a real (forward-declarable) class wrapping the
/// variant, not a bare alias -- \c ScopeBlock::steps holds \c
/// container::vector<Step> while \c Step is still only forward-declared
/// (legal, per the \c std::vector incomplete-type allowance above), and \c
/// Step's own definition (with the variant member) follows \c ScopeBlock,
/// where \c ScopeBlock is by then complete. This preserves the single
/// ORDERED sequence of interleaved build/child-block steps the design
/// requires, and preserves real \c std::variant semantics (\c
/// std::holds_alternative / \c std::get_if / \c std::visit all work on \c
/// Step::value) -- the only change from the brief's literal shape is the one
/// extra wrapping layer forced by the forward-declaration ordering.
///
struct Step {
  std::variant<BuildStep, ScopeBlock> value;

  Step(BuildStep b) : value(std::move(b)) {}   // NOLINT(*-explicit-*)
  Step(ScopeBlock b) : value(std::move(b)) {}  // NOLINT(*-explicit-*)
};

///
/// \brief The whole ordered schedule: the root block plus the total value
/// count (every \c BuildStep::value_id and every \c ScopeBlock::outputs
/// value_id is expected to be < \c num_values; see \c well_formed).
///
struct OrderedSchedule {
  ScopeBlock root{};
  std::size_t num_values = 0;
};

namespace detail {

///
/// \brief \c well_formed's recursive worker: checks \p block's own steps and
/// outputs, recurses into every child \c ScopeBlock step, and checks ordinal
/// uniqueness among \p block's own same-axis (\c IndexSpace::base_key())
/// sibling child blocks.
///
inline bool ordered_schedule_block_well_formed(ScopeBlock const& block,
                                               std::size_t num_values) {
  for (Step const& step : block.steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      if (build->value_id >= num_values) return false;
    } else {
      auto const& child = std::get<ScopeBlock>(step.value);
      if (!ordered_schedule_block_well_formed(child, num_values)) return false;
    }
  }

  // Ordinal uniqueness among same-axis sibling child blocks (direct children
  // of THIS block only; deeper levels are checked by the recursion above).
  for (std::size_t i = 0; i < block.steps.size(); ++i) {
    auto const* ci = std::get_if<ScopeBlock>(&block.steps[i].value);
    if (!ci) continue;
    for (std::size_t j = i + 1; j < block.steps.size(); ++j) {
      auto const* cj = std::get_if<ScopeBlock>(&block.steps[j].value);
      if (!cj) continue;
      if (ci->axis.space().base_key() == cj->axis.space().base_key() &&
          ci->ordinal == cj->ordinal)
        return false;
    }
  }

  for (auto const& [value_id, kind] : block.outputs) {
    (void)kind;
    if (value_id >= num_values) return false;
  }

  return true;
}

}  // namespace detail

///
/// \brief Structural sanity check on \p sched (no sequencer logic): every \c
/// BuildStep::value_id is < \c sched.num_values; ordinals are unique among
/// same-axis (\c IndexSpace::base_key()) sibling blocks within a parent; and
/// every \c ScopeBlock::outputs value_id is < \c sched.num_values.
///
[[nodiscard]] inline bool well_formed(OrderedSchedule const& sched) {
  return detail::ordered_schedule_block_well_formed(sched.root,
                                                    sched.num_values);
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
