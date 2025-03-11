#include "processing.hpp"
#include "format_support.hpp"
#include "utils.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <SeQuant/domain/mbpt/spin.hpp>

#include <spdlog/spdlog.h>

using namespace sequant;

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
			processed                            = simplify(biorthogonal_transform(processed, externals));
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
