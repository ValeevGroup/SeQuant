#include "processing.hpp"
#include "format_support.hpp"
#include "utils.hpp"

#include <SeQuant/core/biorthogonalization.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/expr.hpp>

#include <SeQuant/domain/mbpt/rules/df.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <spdlog/spdlog.h>

using namespace sequant;

container::svector<ResultExpr> postProcess(ResultExpr result,
                                           const IndexSpaceMeta &spaceMeta,
                                           const ProcessingOptions &options) {
  if (result.expression().is<Constant>() ||
      result.expression().is<Variable>()) {
    return {result};
  }

  container::svector<container::svector<Index>> externals =
      result.index_particle_grouping<container::svector<Index>>();

  if (options.density_fitting) {
    IndexSpace aux_space =
        get_default_context().index_space_registry()->retrieve(L"F");
    result.expression() =
        mbpt::density_fit(result.expression(), aux_space, L"g", L"DF");
  }

  switch (options.spintrace) {
    case SpinTracing::None:
      break;
    case SpinTracing::ClosedShell:
    case SpinTracing::Rigorous:
      // Spintracing expects to be fed a sum
      if (!result.expression().is<Sum>()) {
        result.expression() =
            ex<Sum>(ExprPtrList{result.expression()->clone()});
      } else {
        result.expression() = result.expression()->clone();
      }
      break;
  }

  assert(result.expression());

  if (options.expand_symmetrizer) {
    spdlog::debug("Fully expanding all antisymmetrizers before spintracing...");
    result.expression() = simplify(expand_A_op(result.expression()));
    spdlog::debug("Expanded expression:\n{}", result);
  }

  container::svector<ResultExpr> processed;

  switch (options.spintrace) {
    case SpinTracing::None:
      processed = {result};
      break;
    case SpinTracing::ClosedShell:
      processed = closed_shell_spintrace(result);
      break;
    case SpinTracing::Rigorous:
      processed = spintrace(result);
      break;
  }

  for (ResultExpr &current : processed) {
    current.expression() = simplify(current.expression());
    if (options.spintrace != SpinTracing::None) {
      spdlog::debug("Expression after spintracing:\n{}", current);
    }
  }

  // Remove symmetrization operators from expressions
  container::svector<std::optional<ExprPtr>> symmetrizers;
  symmetrizers.reserve(processed.size());
  for (ResultExpr &current : processed) {
    std::optional<ExprPtr> symmetrizer = pop_tensor(current.expression(), L"S");
    if (!symmetrizer.has_value()) {
      symmetrizer = pop_tensor(current.expression(), L"A");
    }

    symmetrizers.push_back(std::move(symmetrizer));
  }

  switch (options.transform) {
    case ProjectionTransformation::None:
      break;
    case ProjectionTransformation::Biorthogonal:
      biorthogonal_transform(processed);

      for (ResultExpr &current : processed) {
        spdlog::debug("Expression after biorthogonal transformation:\n{}",
                      current);
      }
      break;
  }

  if (options.factorize_to_binary) {
    for (ResultExpr &current : processed) {
      current.expression() = optimize(current.expression());
    }
  }

  // Add symmetrizers back into the expressions
  for (std::size_t i = 0; i < processed.size(); ++i) {
    if (symmetrizers.at(i).has_value()) {
      processed.at(i).expression() = ex<Product>(
          ExprPtrList{symmetrizers.at(i).value(), processed.at(i).expression()},
          Product::Flatten::No);
    }
  }

  return processed;
}
