#ifndef SEQUANT_EXTERNAL_INTERFACE_PROCESSING_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PROCESSING_HPP

#include "utils.hpp"

#include <SeQuant/core/expr_fwd.hpp>

enum class SpinTracing {
	None,
	ClosedShell,
	Rigorous,
};

enum class ProjectionTransformation {
	None,
	Biorthogonal,
};

struct ProcessingOptions {
	ProcessingOptions() = default;

	bool density_fitting               = false;
	SpinTracing spintrace              = SpinTracing::Rigorous;
	ProjectionTransformation transform = ProjectionTransformation::None;
	bool factorize_to_binary           = true;
	bool expand_symmetrizer            = false;
};

sequant::ExprPtr postProcess(const sequant::ExprPtr &expression, const IndexSpaceMeta &space_meta,
							 const ProcessingOptions &options = {});

#endif
