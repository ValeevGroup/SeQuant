#ifndef SEQUANT_EXTERNAL_INTERFACE_TEMPLATE_EXPANSION_HPP
#define SEQUANT_EXTERNAL_INTERFACE_TEMPLATE_EXPANSION_HPP

#include <nlohmann/json_fwd.hpp>

namespace sequant::util::extint {

/// Expands all template uses contained in a driver's @p steps array into plain
/// steps, so that the result can be handed to the Executor unchanged.
///
/// A step is a template use if it carries a @c template key (instead of a
/// @c kind); it names a definition in @p templates and binds that definition's
/// @c parameters through its own @c parameters object. Every such use is
/// replaced, in place, by the definition's steps with all @c ${param}
/// placeholders substituted and all internal step IDs (and the references to
/// them) namespaced under the use's @c id. Steps that are not template uses are
/// passed through unchanged.
///
/// @param steps The top-level @c steps array (template uses allowed).
/// @param templates The top-level @c templates object (name -> definition);
///        may be @c null, which is treated as no templates.
/// @returns A new @c steps array free of template uses.
nlohmann::json expand_templates(const nlohmann::json &steps,
                                const nlohmann::json &templates);

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_TEMPLATE_EXPANSION_HPP
