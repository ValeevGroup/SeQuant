#include "processing.hpp"
#include "utils.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/tensor.hpp>

#include <SeQuant/domain/mbpt/spin.hpp>

using namespace sequant;

ExprPtr postProcess(const ExprPtr &expression, const IndexSpaceMeta &spaceMeta, const ProcessingOptions &options) {
	if (expression.is< Constant >() || expression.is< Variable >()) {
		return expression;
	}

	ExprPtr processed;
	if (options.density_fitting) {
		throw std::runtime_error("DF insertion not yet implemented");
	}

	switch (options.spintrace) {
		case SpinTracing::None:
			processed = expression->clone();
			break;
		case SpinTracing::ClosedShell:
		case SpinTracing::Rigorous:
			// Spintracing expects to be fed a sum
			if (!expression.is< Sum >()) {
				processed = ex< Sum >(ExprPtrList{ expression->clone() });
			} else {
				processed = expression->clone();
			}
			break;
	}

	container::svector< container::svector< Index > > externals = getExternalIndexPairs(processed);

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

	switch (options.transform) {
		case ProjectionTransformation::None:
			break;
		case ProjectionTransformation::Biorthogonal:
			processed = simplify(biorthogonal_transform(processed, externals));
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

	return processed;
}
