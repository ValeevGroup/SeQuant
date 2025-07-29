#ifndef SEQUANT_CORE_EXPORT_GENERATOR_HPP
#define SEQUANT_CORE_EXPORT_GENERATOR_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>

#include <string>
#include <string_view>
#include <type_traits>

namespace sequant {

/// The scope at which declarations may happen
enum class DeclarationScope {
  Global,
  Section,
  Expression,
};

/// Abstract base class for all (code) generators. These work by implementing
/// callbacks for semantic actions and their job is to perform (piecewise)
/// translation of those semantic actions into the corresponding code.
template <typename C>
class Generator {
 public:
  using Context = C;
  static_assert(
      std::is_base_of_v<ExportContext, Context>,
      "Generator context class must inherit from sequant::ExportContext");

  virtual ~Generator() = default;

  /// @returns The name of the format this generator will produce (i.e. the name
  /// of the programming language and/or framework)
  virtual std::string get_format_name() const = 0;

  /// @returns Whether this generator supports named sections
  virtual bool supports_named_sections() const = 0;
  /// @returns Whether this generator requires names (that is, it can't deal
  /// with code outside of named sections)
  virtual bool requires_named_sections() const = 0;

  /// @returns The scope at which this generator would like declare indices
  virtual DeclarationScope index_declaration_scope() const = 0;
  /// @returns The scope at which this generator would like declare variables
  virtual DeclarationScope variable_declaration_scope() const = 0;
  /// @returns The scope at which this generator would like declare tensors
  virtual DeclarationScope tensor_declaration_scope() const = 0;

  /// @returns A backend-specific string representation of the given Index
  virtual std::string represent(const Index &idx, const Context &ctx) const = 0;
  /// @returns A backend-specific string representation of the given Tensor
  virtual std::string represent(const Tensor &tensor,
                                const Context &ctx) const = 0;
  /// @returns A backend-specific string representation of the given Variable
  virtual std::string represent(const Variable &variable,
                                const Context &ctx) const = 0;
  /// @returns A backend-specific string representation of the given Constant
  virtual std::string represent(const Constant &constant,
                                const Context &ctx) const = 0;

  /// Semantic callback for creating the given tensor. This is expected to make
  /// the created tensor available for further use.
  virtual void create(const Tensor &tensor, bool zero_init,
                      const Context &ctx) = 0;
  /// Semantic callback for loading the given tensor. Loading implies making the
  /// tensor available for subsequent use (e.g. by reading it from disk in order
  /// to bring it into RAM).
  virtual void load(const Tensor &tensor, bool set_to_zero,
                    const Context &ctx) = 0;
  /// Semantic callback for setting the given tensor to zero
  virtual void set_to_zero(const Tensor &tensor, const Context &ctx) = 0;
  /// Semantic callback for unloading the given tensor. This is expected to
  /// discard the tensor without saving its current value. However, the tensor
  /// must still be loadable afterwards.
  virtual void unload(const Tensor &tensor, const Context &ctx) = 0;
  /// Semantic callback for destroying the given tensor. Destroying means that
  /// this tensor will never be needed again and is no longer required to be
  /// loadable.
  virtual void destroy(const Tensor &tensor, const Context &ctx) = 0;
  /// Semantic callback for persisting the given tensor. Persisting means that
  /// the current value of the tensor is stored in some way that makes its value
  /// available at a later time.
  virtual void persist(const Tensor &tensor, const Context &ctx) = 0;

  /// Semantic callback for creating the given variable. This is expected to
  /// make the created variable available for further use.
  virtual void create(const Variable &variable, bool zero_init,
                      const Context &ctx) = 0;
  /// Semantic callback for loading the given variable. Loading implies making
  /// the variable available for subsequent use (e.g. by reading it from disk in
  /// order to bring it into RAM).
  virtual void load(const Variable &variable, bool set_to_zero,
                    const Context &ctx) = 0;
  /// Semantic callback for setting the given variable to zero
  virtual void set_to_zero(const Variable &variable, const Context &ctx) = 0;
  /// Semantic callback for unloading the given variable. This is expected to
  /// discard the variable without saving its current value. However, the
  /// variable must still be loadable afterwards.
  virtual void unload(const Variable &variable, const Context &ctx) = 0;
  /// Semantic callback for destroying the given variable. Destroying means that
  /// this variable will never be needed again and is no longer required to be
  /// loadable.
  virtual void destroy(const Variable &variable, const Context &ctx) = 0;
  /// Semantic callback for persisting the given variable. Persisting means that
  /// the current value of the variable is stored in some way that makes its
  /// value available at a later time.
  virtual void persist(const Variable &tensor, const Context &ctx) = 0;

  /// Semantic callback for encoding the computation of the given expression.
  /// The result is expected to be written into the provided Tensor.
  virtual void compute(const Expr &expression, const Tensor &result,
                       const Context &ctx) = 0;
  /// Semantic callback for encoding the computation of the given expression.
  /// The result is expected to be written into the provided Variable.
  virtual void compute(const Expr &expression, const Variable &result,
                       const Context &ctx) = 0;

  /// Semantic callback for declaring the provided Index.
  virtual void declare(const Index &idx, const Context &ctx) = 0;
  /// Semantic callback for declaring the provided Variable
  /// @param variable The Variable to declare
  /// @param usage A UsageSet describing how the Variable will be used in the
  /// code
  virtual void declare(const Variable &variable, UsageSet usage,
                       const Context &ctx) = 0;
  /// Semantic callback for declaring the provided Tensor
  /// @param variable The Tensor to declare
  /// @param usage A UsageSet describing how the Tensor will be used in the code
  virtual void declare(const Tensor &tensor, UsageSet usage,
                       const Context &ctx) = 0;

  /// Callback signalling that all indices have been declared
  /// @param amount The amount of indices that have been declared
  virtual void all_indices_declared(std::size_t amount, const Context &ctx) = 0;
  /// Callback signalling that all variables have been declared
  /// @param amount The amount of variables that have been declared
  virtual void all_variables_declared(std::size_t amount,
                                      const Context &ctx) = 0;
  /// Callback signalling that all tensors have been declared
  /// @param amount The amount of tensors that have been declared
  virtual void all_tensors_declared(std::size_t amount, const Context &ctx) = 0;

  /// Callback signalling that we now start declarations for the given scope
  virtual void begin_declarations(DeclarationScope scope,
                                  const Context &ctx) = 0;
  /// Callback signalling that all declarations for the given scope have been
  /// performed
  virtual void end_declarations(DeclarationScope scope, const Context &ctx) = 0;

  /// Semantic callback for inserting a comment at the current position in the
  /// generated code
  /// @param comment The comment to insert
  virtual void insert_comment(const std::string &comment,
                              const Context &ctx) = 0;

  /// Semantic callback that we are starting a new named section
  /// @param name The name of the named section
  virtual void begin_named_section(std::string_view name,
                                   const Context &ctx) = 0;
  /// Semantic callback to signal that we have reached the end of the current
  /// named section
  /// @param name The name of the named section that is ending now
  virtual void end_named_section(std::string_view name, const Context &ctx) = 0;

  /// Callback signalling that we now begin export of a new (high-level)
  /// expression
  virtual void begin_expression(const Context &ctx) = 0;
  /// Callback signalling that the export of the current (high-level) expression
  /// is finished
  virtual void end_expression(const Context &ctx) = 0;

  /// Callback signalling that we now begin with the export process
  virtual void begin_export(const Context &ctx) = 0;
  /// Callback signalling that the export process is finished
  virtual void end_export(const Context &ctx) = 0;

  /// @returns The generated code that corresponds to all semantic actions
  /// encountered during the current export
  virtual std::string get_generated_code() const = 0;

 protected:
  Generator() = default;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_GENERATOR_HPP
