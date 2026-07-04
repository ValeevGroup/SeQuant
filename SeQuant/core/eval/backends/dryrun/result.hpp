#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <any>
#include <array>
#include <cstddef>
#include <memory>
#include <utility>

namespace sequant::eval::dryrun {

///
/// \brief Annotation type DryRun's Result ops decode from the eval engine's
///        std::any [l,r,res] / [pre,post] triples/pairs.
///
/// Unlike \c EvalExprTAPP (opaque \c int64_t index-label hashes -- see
/// \c backends/tapp/eval_expr.hpp), DryRun's annotation IS the plain literal
/// (canon-order) index list itself: \c Result::prod/sum/permute need each
/// index's actual space/extent (via \c CostModel), not just an opaque
/// identity, to compute a modeled size.
///
using annot_t = container::svector<Index>;

class ResultDryRun;
class ResultDryRunNested;

///
/// \brief Builds whichever concrete DryRun Result type matches \p idx's
///        content: a nested \c ResultDryRunNested if any index in \p idx is
///        proto-indexed (a CSV/PNO composite leg, e.g. a CSV amplitude's PNO
///        domain leg `a_1<i_1,i_2>`), otherwise a flat \c ResultDryRun.
///
/// Dispatch is by CONTENT of the decoded result annotation, not by either
/// operand's concrete type -- exactly mirroring how the real eval engine
/// itself decides tensor-of-tensor-ness (\c EvalExpr::tot(), from the same
/// proto-indexed-leg criterion). This is what lets \c prod()/sum() freely
/// combine a flat operand (e.g. a bare 3-center DF integral) with a nested
/// one (e.g. a CSV/PNO coefficient), exactly as real CSV-CCSD terms do,
/// without either side needing to know the other's concrete type.
///
[[nodiscard]] inline ResultPtr make_dryrun_result(
    container::svector<Index> idx, std::shared_ptr<CostModel const> cm,
    ExtentOverrides overrides = {});

namespace detail {

[[nodiscard]] inline bool has_proto(container::svector<Index> const& idx) {
  return std::any_of(idx.begin(), idx.end(),
                     [](Index const& ix) { return ix.has_proto_indices(); });
}

[[nodiscard]] inline ExtentOverrides merge_overrides(ExtentOverrides const& a,
                                                     ExtentOverrides const& b) {
  ExtentOverrides out = a;
  for (auto const& [ix, n] : b) out[ix] = n;
  return out;
}

// Uniform read access to a DryRun Result's (index list, overrides, cost
// model) regardless of which concrete DryRun type `r` is. `is<T>()`/`as<T>()`
// are public Result methods, so no friendship is needed; declared here (and
// defined below, after both concrete classes) purely because their bodies
// need the concrete classes' definitions.
[[nodiscard]] container::svector<Index> indices_of(Result const& r);
[[nodiscard]] ExtentOverrides overrides_of(Result const& r);

///
/// \brief Shared op bodies for the two DryRun Result concrete types.
///
/// Both \c ResultDryRun and \c ResultDryRunNested carry exactly an (index
/// list, ExtentOverrides, CostModel) triple and differ only in what they
/// additionally expose (\c ResultDryRunNested splits its index list into
/// outer()/inner() views for CSV-composite-aware inspection/testing).
/// Implemented once here so the two classes' prod/sum/permute/slice_mode/
/// mode_batches bodies are one-line forwards, not near-duplicated logic.
///
struct DryRunOps {
  [[nodiscard]] static ResultPtr sum(container::svector<Index> const& idx,
                                     ExtentOverrides const& ov,
                                     std::shared_ptr<CostModel const> const& cm,
                                     Result const& other,
                                     std::array<std::any, 3> const& annot) {
    auto const a = Annot<annot_t>{annot};
    auto merged = merge_overrides(ov, overrides_of(other));
    return make_dryrun_result(
        container::svector<Index>(a.this_annot.begin(), a.this_annot.end()), cm,
        std::move(merged));
  }

  [[nodiscard]] static ResultPtr prod(
      container::svector<Index> const& idx, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm, Result const& other,
      std::array<std::any, 3> const& annot) {
    if (other.is<ResultScalar<double>>()) {
      // Scalar * tensor: shape (and any accumulated slicing) unchanged.
      return make_dryrun_result(idx, cm, ov);
    }
    auto const a = Annot<annot_t>{annot};
    auto merged = merge_overrides(ov, overrides_of(other));
    if (a.this_annot.empty()) {
      // Full contraction -> scalar. No real numeric value is ever tracked by
      // this zero-data backend (only sizes/costs), so the placeholder 0.0
      // is never meant to be read as a physical result.
      return eval_result<ResultScalar<double>>(0.0);
    }
    return make_dryrun_result(
        container::svector<Index>(a.this_annot.begin(), a.this_annot.end()), cm,
        std::move(merged));
  }

  [[nodiscard]] static ResultPtr permute(
      container::svector<Index> const& /*idx*/, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm,
      std::array<std::any, 2> const& ann) {
    auto const post = std::any_cast<annot_t>(ann[1]);
    return make_dryrun_result(
        container::svector<Index>(post.begin(), post.end()), cm, ov);
  }

  [[nodiscard]] static ResultPtr slice_mode(
      container::svector<Index> const& idx, ExtentOverrides const& ov,
      std::shared_ptr<CostModel const> const& cm, std::size_t mode,
      std::size_t elem_lo, std::size_t elem_hi) {
    SEQUANT_ASSERT(mode < idx.size());
    auto merged = ov;
    merged[idx[mode]] = elem_hi - elem_lo;
    return make_dryrun_result(idx, cm, std::move(merged));
  }

  [[nodiscard]] static container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(container::svector<Index> const& idx, ExtentOverrides const& ov,
               std::shared_ptr<CostModel const> const& cm, std::size_t mode,
               std::size_t target_batch_size) {
    SEQUANT_ASSERT(mode < idx.size());
    Index const& ix = idx[mode];
    std::size_t extent;
    if (auto it = ov.find(ix); it != ov.end())
      extent = it->second;
    else
      extent = cm->regime().extent(ix);

    container::svector<std::pair<std::size_t, std::size_t>> out;
    if (target_batch_size == 0 || extent == 0) {
      out.push_back({0, extent});
      return out;
    }
    for (std::size_t lo = 0; lo < extent; lo += target_batch_size)
      out.push_back({lo, std::min(extent, lo + target_batch_size)});
    return out;
  }
};

}  // namespace detail

///
/// \brief Flat (non-CSV) zero-data tensor token.
///
/// Carries only its own literal outer index list (canon order -- the same
/// order \c EvalExpr::canon_indices()/annot() use, so \c slice_mode()/
/// \c mode_batches()'s positional `mode` argument indexes it correctly), an
/// \c ExtentOverrides table recording any runtime \c slice_mode()/
/// \c mode_batches() narrowing (keyed by Index so it survives reshaping
/// across prod/sum/permute), and a shared \c CostModel. No tensor data is
/// ever allocated or copied; every op is index-set bookkeeping plus a
/// CostModel query. Mirrors \c ResultTensorTAPP's structure
/// (backends/tapp/result.hpp) with every real-tensor line replaced by that
/// bookkeeping.
///
class ResultDryRun final : public Result {
 public:
  using Result::id_t;

  ResultDryRun(container::svector<Index> idxset,
               std::shared_ptr<CostModel const> cm,
               ExtentOverrides overrides = {})
      : Result{Payload{}},
        indices_{std::move(idxset)},
        cm_{std::move(cm)},
        overrides_{std::move(overrides)} {}

  [[nodiscard]] container::svector<Index> const& indices() const noexcept {
    return indices_;
  }
  [[nodiscard]] ExtentOverrides const& overrides() const noexcept {
    return overrides_;
  }

 private:
  struct Payload {};

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<ResultDryRun>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    return detail::DryRunOps::sum(indices_, overrides_, cm_, other, annot);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    return detail::DryRunOps::prod(indices_, overrides_, cm_, other, annot);
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann);
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann);
  }

  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    return detail::DryRunOps::slice_mode(indices_, overrides_, cm_, mode,
                                         elem_lo, elem_hi);
  }

  [[nodiscard]] container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(std::size_t mode, std::size_t target_batch_size) const override {
    return detail::DryRunOps::mode_batches(indices_, overrides_, cm_, mode,
                                           target_batch_size);
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<ResultDryRun>() || other.is<ResultDryRunNested>());
    overrides_ =
        detail::merge_overrides(overrides_, detail::overrides_of(other));
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_);
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_);
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t /*factor*/) const override {
    return eval_result<ResultDryRun>(indices_, cm_, overrides_);
  }

  [[nodiscard]] std::size_t size_in_bytes() const final {
    return cm_->memsize(indices_, overrides_);
  }

  container::svector<Index> indices_;
  std::shared_ptr<CostModel const> cm_;
  ExtentOverrides overrides_;
};

///
/// \brief CSV/PNO tensor-of-tensor zero-data token.
///
/// Like \c ResultDryRun, but additionally exposes an outer()/inner() split of
/// its (canon-order) index list -- inner = the proto-indexed (composite)
/// legs, e.g. a CSV amplitude's PNO domain leg `a_1<i_1,i_2>`; outer = every
/// other (plain) leg, e.g. the PAO index `mu~_1`. The split is purely an
/// observability/testing convenience: \c size_in_bytes()'s arithmetic is
/// IDENTICAL to \c ResultDryRun's (\c CostModel::memsize already routes any
/// index list containing a proto-indexed entry through the moment-aware
/// `inner_pow` path internally, via \c tot_indices/inner_aware_volume --
/// content-driven, not type-driven), so tests that want to confirm "this used
/// the k-th moment, not extent^k" can inspect inner() directly.
///
/// Position semantics for \c slice_mode()/\c mode_batches(): the `mode`
/// argument the runtime passes is always resolved against the FULL
/// canon-order list (an optional trailing constructor argument, defaulting to
/// `outer ++ inner` when the caller does not need position accuracy, e.g. a
/// hand-built test instance); the \c DryRunLeafEvaluator (eval_expr.hpp)
/// always supplies the leaf's true \c canon_indices() order there, since only
/// LEAF-constructed instances are ever sliced by the runtime (\c slice_mode()
/// is invoked only inside the batched evaluator's leaf-wrapping closure, never
/// on a prod()/sum()-produced intermediate).
///
class ResultDryRunNested final : public Result {
 public:
  using Result::id_t;

  ResultDryRunNested(container::svector<Index> outer,
                     container::svector<Index> inner,
                     std::shared_ptr<CostModel const> cm,
                     ExtentOverrides overrides = {},
                     container::svector<Index> canon_order = {})
      : Result{Payload{}},
        outer_{std::move(outer)},
        inner_{std::move(inner)},
        indices_{canon_order.empty()
                     ? [this] {
                         container::svector<Index> c = outer_;
                         c.insert(c.end(), inner_.begin(), inner_.end());
                         return c;
                       }()
                     : std::move(canon_order)},
        cm_{std::move(cm)},
        overrides_{std::move(overrides)} {}

  [[nodiscard]] container::svector<Index> const& outer() const noexcept {
    return outer_;
  }
  [[nodiscard]] container::svector<Index> const& inner() const noexcept {
    return inner_;
  }
  [[nodiscard]] container::svector<Index> const& indices() const noexcept {
    return indices_;
  }
  [[nodiscard]] ExtentOverrides const& overrides() const noexcept {
    return overrides_;
  }

 private:
  struct Payload {};

  [[nodiscard]] id_t type_id() const noexcept override {
    return id_for_type<ResultDryRunNested>();
  }

  [[nodiscard]] ResultPtr sum(
      Result const& other,
      std::array<std::any, 3> const& annot) const override {
    return detail::DryRunOps::sum(indices_, overrides_, cm_, other, annot);
  }

  [[nodiscard]] ResultPtr prod(Result const& other,
                               std::array<std::any, 3> const& annot,
                               DeNest /*DeNestFlag*/) const override {
    return detail::DryRunOps::prod(indices_, overrides_, cm_, other, annot);
  }

  [[nodiscard]] ResultPtr permute(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann);
  }

  [[nodiscard]] ResultPtr adjoint(
      std::array<std::any, 2> const& ann) const override {
    return detail::DryRunOps::permute(indices_, overrides_, cm_, ann);
  }

  [[nodiscard]] ResultPtr slice_mode(std::size_t mode, std::size_t elem_lo,
                                     std::size_t elem_hi) const override {
    return detail::DryRunOps::slice_mode(indices_, overrides_, cm_, mode,
                                         elem_lo, elem_hi);
  }

  [[nodiscard]] container::svector<std::pair<std::size_t, std::size_t>>
  mode_batches(std::size_t mode, std::size_t target_batch_size) const override {
    return detail::DryRunOps::mode_batches(indices_, overrides_, cm_, mode,
                                           target_batch_size);
  }

  void add_inplace(Result const& other) override {
    SEQUANT_ASSERT(other.is<ResultDryRun>() || other.is<ResultDryRunNested>());
    overrides_ =
        detail::merge_overrides(overrides_, detail::overrides_of(other));
  }

  [[nodiscard]] ResultPtr symmetrize() const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_);
  }

  [[nodiscard]] ResultPtr antisymmetrize(size_t /*bra_rank*/) const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_);
  }

  [[nodiscard]] ResultPtr mult_by_phase(std::int8_t /*factor*/) const override {
    return eval_result<ResultDryRunNested>(outer_, inner_, cm_, overrides_,
                                           indices_);
  }

  [[nodiscard]] std::size_t size_in_bytes() const final {
    return cm_->memsize(indices_, overrides_);
  }

  container::svector<Index> outer_;
  container::svector<Index> inner_;
  container::svector<Index> indices_;  // canon order; outer_++inner_ content
  std::shared_ptr<CostModel const> cm_;
  ExtentOverrides overrides_;
};

[[nodiscard]] inline ResultPtr make_dryrun_result(
    container::svector<Index> idx, std::shared_ptr<CostModel const> cm,
    ExtentOverrides overrides) {
  if (!detail::has_proto(idx))
    return eval_result<ResultDryRun>(std::move(idx), std::move(cm),
                                     std::move(overrides));
  container::svector<Index> outer, inner;
  for (auto const& ix : idx)
    (ix.has_proto_indices() ? inner : outer).push_back(ix);
  return eval_result<ResultDryRunNested>(std::move(outer), std::move(inner),
                                         std::move(cm), std::move(overrides),
                                         std::move(idx));
}

namespace detail {

[[nodiscard]] inline container::svector<Index> indices_of(Result const& r) {
  if (r.is<ResultDryRun>()) return r.as<ResultDryRun>().indices();
  SEQUANT_ASSERT(r.is<ResultDryRunNested>());
  return r.as<ResultDryRunNested>().indices();
}

[[nodiscard]] inline ExtentOverrides overrides_of(Result const& r) {
  if (r.is<ResultDryRun>()) return r.as<ResultDryRun>().overrides();
  SEQUANT_ASSERT(r.is<ResultDryRunNested>());
  return r.as<ResultDryRunNested>().overrides();
}

}  // namespace detail

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_RESULT_HPP
