// Unit tests for the shared slicing-signature helper
// (SeQuant/core/eval/slicing_signature.hpp): index_position, slicing_signature,
// and signatures_consistent -- the criterion that decides whether a CSE-folded
// value's occurrences may share one sliced materialization (consistent
// signature) or must be SPLIT (a relabeled mode diverges). See
// doc/dev/specs/2026-08-07-remat-cse-aware-split-design.md.

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <optional>

namespace {

using sequant::bra;
using sequant::ColumnSymmetry;
using sequant::ex;
using sequant::ExprPtr;
using sequant::Index;
using sequant::ket;
using sequant::Symmetry;
using sequant::Tensor;
using sequant::container::svector;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;

// Distinctly named (Unity-build safe): a single-tensor leaf eval node.
EvalNodeDryRun sig_test_leaf(ExprPtr const& t) {
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(t);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  return node;
}

ExprPtr sig_test_g(svector<Index> const& idxs) {
  return ex<Tensor>(L"g", bra(idxs), ket{}, Symmetry::Nonsymm, std::nullopt,
                    ColumnSymmetry::Nonsymm);
}

}  // namespace

TEST_CASE("slicing_signature: relabeled mode diverges, shared-label agrees",
          "[slicing-signature]") {
  Index i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"}, i4{L"i_4"};
  // A carries i_3 on a result slot; B carries i_4 there (the g.C legs' shape).
  auto A = sig_test_leaf(sig_test_g({i1, i2, i3}));
  auto B = sig_test_leaf(sig_test_g({i1, i2, i4}));

  // index_position: a mode present on the result has a position; absent ->
  // none.
  CHECK(sequant::index_position(A, i3).has_value());
  CHECK_FALSE(sequant::index_position(A, i4).has_value());
  CHECK(sequant::index_position(B, i4).has_value());
  CHECK_FALSE(sequant::index_position(B, i3).has_value());

  svector<EvalNodeDryRun> occs{A, B};

  // A RELABELED mode (present in A, absent in B) -> signatures inconsistent ->
  // the value cannot be shared sliced and must be SPLIT.
  CHECK_FALSE(sequant::signatures_consistent(occs, svector<Index>{i3}));
  CHECK_FALSE(sequant::signatures_consistent(occs, svector<Index>{i4}));

  // A SHARED-LABEL mode (i_1 occupies the same canonical slot in both, since A
  // and B are canonically equivalent up to the i_3<->i_4 relabel) -> consistent
  // -> shareable while folded.
  CHECK(sequant::signatures_consistent(occs, svector<Index>{i1}));

  // A single occurrence is trivially consistent.
  svector<EvalNodeDryRun> one{A};
  CHECK(sequant::signatures_consistent(one, svector<Index>{i3}));
}
