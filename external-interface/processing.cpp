#include "processing.hpp"
#include "format_support.hpp"
#include "utils.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <SeQuant/domain/mbpt/spin.hpp>

#include <spdlog/spdlog.h>

using namespace sequant;

ExprPtr custom_biorthogonalize(ExprPtr expr, const container::svector<container::svector<Index>> &externals) {
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

	container::map<Index, Index> permutation;
	if (braSpacesSame) {
		permutation = {
			{{externals[0][0], externals[1][0]}, {externals[1][0], externals[0][0]}}
		};
	} else {
		assert(ketSpacesSame);
		permutation = {
			{{externals[0][1], externals[1][1]}, {externals[1][1], externals[0][1]}}
		};
	}

	ExprPtr permuted = expr.clone();
	permuted = transform_expr(permuted, permutation);

	spdlog::warn("Adding additional factor of 1/2 to P1/S2 expression to remain consistent with GeCCo");

	return ex<Constant>(ratio{1,2}) * ex<Constant>(ratio{1,3}) * ex<Sum>(ex<Constant>(2) * expr + permuted);
}

ResultExpr postProcess(ResultExpr result, const IndexSpaceMeta &spaceMeta, const ProcessingOptions &options) {
	if (result.expression().is< Constant >() || result.expression().is< Variable >()) {
		return result;
	}

	ExprPtr processed;
	const container::svector< container::svector< Index > > externals =
		result.index_particle_grouping< container::svector< Index > >();


	if (options.density_fitting) {
		throw std::runtime_error("DF insertion not yet implemented");
	}

	switch (options.spintrace) {
		case SpinTracing::None:
			processed = result.expression()->clone();
			break;
		case SpinTracing::ClosedShell:
		case SpinTracing::Rigorous:
			// Spintracing expects to be fed a sum
			if (!result.expression().is< Sum >()) {
				processed = ex< Sum >(ExprPtrList{ result.expression()->clone() });
			} else {
				processed = result.expression()->clone();
			}
			break;
	}

	assert(processed);

	if (options.expand_symmetrizer) {
		spdlog::debug("Fully expanding all antisymmetrizers before spintracing...");
		processed = simplify(expand_A_op(processed));
		spdlog::debug("Expanded expression:\n{}", processed);
	}

	switch (options.spintrace) {
		case SpinTracing::None:
			break;
		case SpinTracing::ClosedShell:
			processed = closed_shell_spintrace(processed, externals);
			break;
		case SpinTracing::Rigorous:
			processed = spintrace(processed, externals);
			break;
	}

	processed = simplify(processed);
	if (options.spintrace != SpinTracing::None) {
		spdlog::debug("Expression after spintracing:\n{}", processed);
	}

	switch (options.transform) {
		case ProjectionTransformation::None:
			break;
		case ProjectionTransformation::Biorthogonal:
			// TODO: pop S tensor first?
			std::optional< ExprPtr > symmetrizer = popTensor(processed, L"S");
			processed                            = custom_biorthogonalize(processed, externals);
			spdlog::debug("Biorthogonalized without simplification:\n{}", processed);
			processed = simplify(processed);
			if (symmetrizer) {
				processed =
					simplify(ex< Product >(ExprPtrList{ symmetrizer.value(), processed }, Product::Flatten::No));
			}
			spdlog::debug("Expression after biorthogonal transformation:\n{}", processed);
			break;
	}

	if (options.factorize_to_binary) {
		std::optional< ExprPtr > symmetrizer = popTensor(processed, L"S");
		if (!symmetrizer.has_value()) {
			symmetrizer = popTensor(processed, L"A");
		}

		processed = optimize(processed);

		if (symmetrizer.has_value()) {
			processed = ex< Product >(1, ExprPtrList{ symmetrizer.value(), processed });
		}
	}

	result.expression() = std::move(processed);

	return result;
}
