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

void custom_biorthogonalize(container::svector< ResultExpr > &exprs) {
	container::svector< container::svector< Index > > externals =
		exprs.at(0).index_particle_grouping< container::svector< Index > >();

	if (externals.size() < 2) {
		// Singles or scalar
		for (ResultExpr &current : exprs) {
			std::optional< ExprPtr > symmetrizer = popTensor(current.expression(), L"S");
			current.expression()                 = biorthogonal_transform(current.expression(), externals);
			current.expression()                 = simplify(current.expression());
			if (symmetrizer) {
				current.expression() = simplify(
					ex< Product >(ExprPtrList{ symmetrizer.value(), current.expression() }, Product::Flatten::No));
			}
		}

		return;
	}

	assert(externals.size() == 2);

	const bool braSpacesSame = externals[0][0].space() == externals[1][0].space();
	const bool ketSpacesSame = externals[0][1].space() == externals[1][1].space();

	if (braSpacesSame && ketSpacesSame) {
		// P0, P2, I2
		for (ResultExpr &current : exprs) {
			std::optional< ExprPtr > symmetrizer = popTensor(current.expression(), L"S");
			current.expression()                 = biorthogonal_transform(current.expression(), externals);
			current.expression()                 = simplify(current.expression());
			if (symmetrizer) {
				current.expression() = simplify(
					ex< Product >(ExprPtrList{ symmetrizer.value(), current.expression() }, Product::Flatten::No));
			}
		}

		return;
	}

	if (!braSpacesSame && !ketSpacesSame) {
		// S1

		assert(exprs.size() == 2);

		ExprPtr clone = exprs[0].expression().clone();

		exprs[0].expression() =
			simplify(ex< Constant >(ratio{ 1, 6 }) * (ex< Constant >(2) * clone + exprs[1].expression()));
		exprs[1].expression() =
			simplify(ex< Constant >(ratio{ 1, 6 }) * (ex< Constant >(2) * exprs[1].expression() + clone));

		return;
	}

	// P1, I1

	for (ResultExpr &current : exprs) {
		externals = current.index_particle_grouping< container::svector< Index > >();

		container::map< Index, Index > permutation;
		if (braSpacesSame) {
			permutation = { { { externals[0][0], externals[1][0] }, { externals[1][0], externals[0][0] } } };
		} else {
			assert(ketSpacesSame);
			permutation = { { { externals[0][1], externals[1][1] }, { externals[1][1], externals[0][1] } } };
		}

		std::optional< ExprPtr > symmetrizer = popTensor(current.expression(), L"S");

		ExprPtr permuted = current.expression().clone();
		permuted         = transform_expr(permuted, permutation);

		spdlog::warn("Adding additional factor of 1/2 to P1/S2 expression to remain consistent with GeCCo");

		current.expression() = ex< Constant >(ratio{ 1, 2 }) * ex< Constant >(ratio{ 1, 3 })
							   * ex< Sum >(ex< Constant >(2) * current.expression() + permuted);
		current.expression() = simplify(current.expression());

		if (symmetrizer) {
			current.expression() =
				simplify(ex< Product >(ExprPtrList{ symmetrizer.value(), current.expression() }, Product::Flatten::No));
		}
	}
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
		current.expression() = simplify(current.expression());
		if (options.spintrace != SpinTracing::None) {
			spdlog::debug("Expression after spintracing:\n{}", current);
		}
	}

	switch (options.transform) {
		case ProjectionTransformation::None:
			break;
		case ProjectionTransformation::Biorthogonal:
			custom_biorthogonalize(processed);

			for (const ResultExpr &current : processed) {
				spdlog::debug("Expression after biorthogonal transformation:\n{}", current);
			}
			break;
	}

	for (ResultExpr &current : processed) {
		if (options.factorize_to_binary) {
			std::optional< ExprPtr > symmetrizer = popTensor(current.expression(), L"S");
			if (!symmetrizer.has_value()) {
				symmetrizer = popTensor(current.expression(), L"A");
			}

			current.expression() = optimize(current.expression());

			if (symmetrizer.has_value()) {
				current.expression() =
					ex< Product >(ExprPtrList{ symmetrizer.value(), current.expression() }, Product::Flatten::No);
			}
		}
	}

	return processed;
}
