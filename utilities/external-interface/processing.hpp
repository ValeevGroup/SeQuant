#ifndef SEQUANT_EXTERNAL_INTERFACE_PROCESSING_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PROCESSING_HPP

#include "utils.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

enum class SpinTracing {
  None,
  ClosedShell,
  Rigorous,
};

enum class ProjectionTransformation {
  None,
  Biorthogonal,
};

enum class IndexBatching {
  None,
  Fastest,
  Slowest,
};

struct ProcessingOptions {
  ProcessingOptions() = default;

  bool density_fitting = false;
  SpinTracing spintrace = SpinTracing::Rigorous;
  ProjectionTransformation transform = ProjectionTransformation::None;
  bool factorize_to_binary = true;
  bool expand_symmetrizer = false;
  bool term_by_term = false;
  bool subexpression_elimination = true;
  std::size_t min_cse_usage = 2;
  std::variant<IndexBatching, std::vector<sequant::Index>> batching =
      IndexBatching::Slowest;
  std::size_t min_unbatched_indices = 2;
};

sequant::container::svector<sequant::ResultExpr> postProcess(
    sequant::ResultExpr expression, const ProcessingOptions &options = {});

#endif
