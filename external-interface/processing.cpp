#include "processing.hpp"
#include "format_support.hpp"
#include "utils.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <SeQuant/domain/mbpt/spin.hpp>

#include <spdlog/spdlog.h>

using namespace sequant;

ExprPtr custom_biorthogonalize(ExprPtr expr, const container::svector< container::svector< Index > > &externals) {
	assert(externals.size() <= 2);

	if (externals.size() < 2) {
		return biorthogonal_transform(expr, externals);
	}

	const bool braSpacesSame = externals[0][0].space() == externals[1][0].space();
	const bool ketSpacesSame = externals[0][1].space() == externals[1][1].space();
	if (braSpacesSame && ketSpacesSame) {
		// P0, P2, I2
		return biorthogonal_transform(expr, externals);
	}

	if (!braSpacesSame && !ketSpacesSame) {
		spdlog::error("Biorthogonalization for S1 not yet supported");
		return expr;
	}

	container::map< Index, Index > permutation;
	if (braSpacesSame) {
		permutation = { { { externals[0][0], externals[1][0] }, { externals[1][0], externals[0][0] } } };
	} else {
		assert(ketSpacesSame);
		permutation = { { { externals[0][1], externals[1][1] }, { externals[1][1], externals[0][1] } } };
	}

	ExprPtr permuted = expr.clone();
	permuted         = transform_expr(permuted, permutation);

	spdlog::warn("Adding additional factor of 1/2 to P1/S2 expression to remain consistent with GeCCo");

	return ex< Constant >(ratio{ 1, 2 }) * ex< Constant >(ratio{ 1, 3 })
		   * ex< Sum >(ex< Constant >(2) * expr + permuted);
}

container::svector< ResultExpr > postProcess(ResultExpr result, const IndexSpaceMeta &spaceMeta,
											 const ProcessingOptions &options) {
	if (result.expression().is< Constant >() || result.expression().is< Variable >()) {
		return { result };
	}

	container::svector< container::svector< Index > > externals =
		result.index_particle_grouping< container::svector< Index > >();


	if (options.density_fitting) {
		throw std::runtime_error("DF insertion not yet implemented");
	}

	switch (options.spintrace) {
		case SpinTracing::None:
			break;
		case SpinTracing::ClosedShell:
		case SpinTracing::Rigorous:
			// Spintracing expects to be fed a sum
			if (!result.expression().is< Sum >()) {
				result.expression() = ex< Sum >(ExprPtrList{ result.expression()->clone() });
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

	container::svector< ResultExpr > processed;

	switch (options.spintrace) {
		case SpinTracing::None:
			break;
		case SpinTracing::ClosedShell:
			processed = closed_shell_spintrace(result);
			break;
		case SpinTracing::Rigorous:
			processed.push_back(spintrace(result));
			break;
	}

	for (ResultExpr &current : processed) {
		externals = current.index_particle_grouping< container::svector< Index > >();

		current.expression() = simplify(current.expression());
		if (options.spintrace != SpinTracing::None) {
			spdlog::debug("Expression after spintracing:\n{}", current);
		}

		switch (options.transform) {
			case ProjectionTransformation::None:
				break;
			case ProjectionTransformation::Biorthogonal:
				// TODO: pop S tensor first?
				std::optional< ExprPtr > symmetrizer = popTensor(current.expression(), L"S");
				current.expression()                 = custom_biorthogonalize(current.expression(), externals);
				spdlog::debug("Biorthogonalized without simplification:\n{}", current);
				current.expression() = simplify(current.expression());
				if (symmetrizer) {
					current.expression() = simplify(
						ex< Product >(ExprPtrList{ symmetrizer.value(), current.expression() }, Product::Flatten::No));
				}
				spdlog::debug("Expression after biorthogonal transformation:\n{}", current);
				break;
		}

		if (options.factorize_to_binary) {
			std::optional< ExprPtr > symmetrizer = popTensor(current.expression(), L"S");
			if (!symmetrizer.has_value()) {
				symmetrizer = popTensor(current.expression(), L"A");
			}

			current.expression() = optimize(current.expression());

			if (symmetrizer.has_value()) {
				current.expression() = ex< Product >(1, ExprPtrList{ symmetrizer.value(), current.expression() });
			}
		}
	}

	return processed;
}
