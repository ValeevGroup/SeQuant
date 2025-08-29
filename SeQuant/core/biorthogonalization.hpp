#ifndef SEQUANT_BIORTHOGONALIZE_HPP
#define SEQUANT_BIORTHOGONALIZE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

namespace sequant {

namespace {
static constexpr double default_biorth_threshold = 1e-12;
}

[[nodiscard]] ResultExpr biorthogonal_transform_copy(
    const ResultExpr& expr, double threshold = default_biorth_threshold);

[[nodiscard]] container::svector<ResultExpr> biorthogonal_transform_copy(
    const container::svector<ResultExpr>& exprs,
    double threshold = default_biorth_threshold);

void biorthogonal_transform(ResultExpr& expr,
                            double threshold = default_biorth_threshold);

void biorthogonal_transform(container::svector<ResultExpr>& exprs,
                            double threshold = default_biorth_threshold);

/// performs symbolic biorthogonal transform of CC-like equation using
///(for rank-3 and higher
/// Wang-Knizia biorthogonalization (https://arxiv.org/abs/1805.00565) is used
[[nodiscard]] ExprPtr biorthogonal_transform(
    const ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups = {},
    double threshold = default_biorth_threshold);

}  // namespace sequant

#endif
