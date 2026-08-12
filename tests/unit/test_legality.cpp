#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>

TEST_CASE(
    "build_site_of collects the batchable contracted axis of a binary "
    "contraction and excludes a non-batchable one",
    "[legality]") {
  using sequant::eval::dryrun::EvalExprDryRun;

  // g{i_1;a_3} * t{a_3;i_1}: contracts BOTH i_1 (occupied) and a_3
  // (virtual), fully -- so canon_indices() on the product node is empty and
  // contracted_indices() reports both i_1 and a_3. Mark only a_3 batchable
  // so the assertion distinguishes "in the build-site" from "in the
  // subtree at all".
  auto expr = sequant::deserialize<sequant::ExprPtr>("g{i_1;a_3} * t{a_3;i_1}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"a";
  };

  auto const site = sequant::eval::build_site_of(node, policy);

  auto const contains = [&](std::wstring const& key) {
    return std::any_of(site.begin(), site.end(), [&](sequant::Index const& ix) {
      return ix.space().base_key() == key;
    });
  };
  CHECK(contains(L"a"));
  CHECK_FALSE(contains(L"i"));
}
