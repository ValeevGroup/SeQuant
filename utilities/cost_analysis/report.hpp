#ifndef SEQUANT_UTILITIES_COST_ANALYSIS_REPORT_HPP
#define SEQUANT_UTILITIES_COST_ANALYSIS_REPORT_HPP

#include "cost_analysis.hpp"

#include <ostream>
#include <string>
#include <utility>
#include <vector>

/// Serializes a binarized tree back to a parenthesized string (also used by the
/// pipeline's --dump_tree output).
std::string full_expr(const TreeNode& n);

/// Writes the full Markdown report (per-result tables + cache section) to \p
/// out.
void write_report(
    const Config& cfg,
    const std::vector<std::pair<std::string, EquationResult>>& results,
    const SimResult& sim, std::ostream& out);

#endif  // SEQUANT_UTILITIES_COST_ANALYSIS_REPORT_HPP
