#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "processing_step.hpp"

#include "test_extint.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <cstddef>
#include <string>
#include <string_view>

using namespace sequant;
using namespace sequant::util::extint;

namespace {

/// Removes all expressions from its input. Optionally reports the input as
/// unchanged, in which case it is aliased by the base class.
class ClearingStep final
    : public OneByOneProcessingStep<ExpressionData, true,
                                    InputAccess::WorkingCopy> {
 public:
  explicit ClearingStep(bool report_unchanged = false,
                        bool keep_original = false)
      : report_unchanged_(report_unchanged), keep_original_(keep_original) {}

  std::string kind() const override { return "clearing"; }
  bool accepts_options() const override { return false; }
  bool requires_options() const override { return false; }
  void set_options(const nlohmann::json &) override {}

  /// Whether the input handed to the last process call was a copy
  bool last_was_copy = false;

 protected:
  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      WorkingCopy<ExpressionData> input) override {
    last_was_copy = input.is_copy();

    if (report_unchanged_) {
      return 0;
    }

    input.data.expressions.clear();
    ctx.set_data(id_prefix, id_start, std::move(input.data));
    return 1;
  }

  bool alias_unchanged_inputs() const override { return report_unchanged_; }
  bool needs_original_input() const override { return keep_original_; }

 private:
  bool report_unchanged_;
  bool keep_original_;
};

ExecutionContext seeded_context() {
  ExecutionContext ctx;
  ctx.set_data("in", 0,
               ProcessingData{
                   test_support::expr_data({L"R1{a1;i1} = + f{a1;i1}:A-H-S"})});
  ctx.add_data_alias("in.0", "in.energy");
  return ctx;
}

std::size_t n_expressions(ExecutionContext &ctx, std::string_view id) {
  return convert_data<ExpressionData>(ctx.get_data(id).front().data.get())
      .expressions.size();
}

}  // namespace

TEST_CASE("OneByOneProcessingStep working copies", "[extint]") {
  ExecutionContext ctx = seeded_context();

  SECTION("input needed later is copied") {
    ctx.set_needed_later_predicate([](std::string_view) { return true; });

    ClearingStep step;
    REQUIRE(step.run("out", ctx, {"in.0"}) == 1);
    REQUIRE(step.last_was_copy);
    REQUIRE(n_expressions(ctx, "in.0") == 1);
    REQUIRE(n_expressions(ctx, "out.0") == 0);
    REQUIRE(ctx.has_data("out.energy"));
  }

  SECTION("input not needed later is consumed") {
    ctx.set_needed_later_predicate([](std::string_view) { return false; });

    ClearingStep step;
    REQUIRE(step.run("out", ctx, {"in.0"}) == 1);
    REQUIRE_FALSE(step.last_was_copy);
    REQUIRE_THROWS_AS(ctx.get_data("in.0"), Exception);
    REQUIRE(n_expressions(ctx, "out.0") == 0);
    REQUIRE(ctx.has_data("out.energy"));
  }

  SECTION("the original can be required") {
    ctx.set_needed_later_predicate([](std::string_view) { return false; });

    ClearingStep step(false, /* keep_original = */ true);
    REQUIRE(step.run("out", ctx, {"in.0"}) == 1);
    REQUIRE(step.last_was_copy);
  }

  SECTION("unchanged inputs are aliased") {
    ClearingStep step(/* report_unchanged = */ true);
    REQUIRE(step.run("out", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.ids_are_equivalent("out.0", "in.0"));
    REQUIRE(ctx.has_data("out.energy"));
  }

  SECTION("a consumed input can't be reported as unchanged") {
    ctx.set_needed_later_predicate([](std::string_view) { return false; });

    ClearingStep step(/* report_unchanged = */ true);
    REQUIRE_THROWS_WITH(step.run("out", ctx, {"in.0"}),
                        Catch::Matchers::ContainsSubstring("consumed"));
  }

  SECTION("errors name the input") {
    ctx.set_data("tree", 0,
                 ProcessingData{test_support::tree_data(
                     {L"X{a1;i1} = + f{a1;i1}:A-H-S"})});

    ClearingStep step;
    REQUIRE_THROWS_WITH(step.run("out", ctx, {"tree.0"}),
                        Catch::Matchers::ContainsSubstring("tree.0"));
  }
}
