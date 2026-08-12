// Task 2 of the ordered-scope batched-eval design (SP2): pins the
// OrderedSchedule IR (SeQuant/core/eval/ordered_schedule.hpp) -- an ORDERED
// tree of loop blocks and build steps -- plus its well_formed structural
// sanity check. No sequencer/executor here.

#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/index.hpp>

#include <catch2/catch_test_macros.hpp>

#include <utility>

using sequant::Index;
using sequant::eval::BuildStep;
using sequant::eval::OrderedSchedule;
using sequant::eval::OutputKind;
using sequant::eval::ScopeBlock;
using sequant::eval::Step;
using sequant::eval::well_formed;

TEST_CASE(
    "well_formed accepts a root block with one build step and one child "
    "loop block carrying an AccumulateSum output",
    "[ordered-schedule]") {
  Index const i1{L"i_1"};

  ScopeBlock child;
  child.axis = i1;
  child.ordinal = 0;
  child.steps.push_back(Step{BuildStep{1}});
  child.outputs.push_back({1, OutputKind::AccumulateSum});

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{BuildStep{0}});
  sched.root.steps.push_back(Step{std::move(child)});
  sched.num_values = 2;

  CHECK(well_formed(sched));
}

TEST_CASE("well_formed rejects an out-of-range BuildStep::value_id",
          "[ordered-schedule]") {
  OrderedSchedule sched;
  sched.root.steps.push_back(Step{BuildStep{5}});
  sched.num_values = 1;  // value_id 5 is out of range

  CHECK_FALSE(well_formed(sched));
}

TEST_CASE(
    "well_formed rejects duplicate ordinals among same-axis sibling blocks",
    "[ordered-schedule]") {
  Index const i1{L"i_1"};
  Index const i2{L"i_2"};  // same TYPE ("i") as i1, different physical label

  ScopeBlock child_a;
  child_a.axis = i1;
  child_a.ordinal = 0;
  child_a.steps.push_back(Step{BuildStep{0}});

  ScopeBlock child_b;
  child_b.axis = i2;
  child_b.ordinal = 0;  // duplicate ordinal at the same axis TYPE as child_a
  child_b.steps.push_back(Step{BuildStep{1}});

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{std::move(child_a)});
  sched.root.steps.push_back(Step{std::move(child_b)});
  sched.num_values = 2;

  CHECK_FALSE(well_formed(sched));
}

TEST_CASE("well_formed rejects an out-of-range output value_id",
          "[ordered-schedule]") {
  ScopeBlock child;
  child.axis = Index{L"i_1"};
  child.outputs.push_back({7, OutputKind::AccumulateScatter});

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{std::move(child)});
  sched.num_values = 1;  // output value_id 7 is out of range

  CHECK_FALSE(well_formed(sched));
}
